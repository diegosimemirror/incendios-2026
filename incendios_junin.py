# ============================================================
# INCENDIOS EN JUNÍN · ¿A QUIÉN AFECTA?
# Pregunta: ¿cuántas personas y qué establecimientos de salud están a 5 km
#           o menos de los focos de calor detectados en los últimos 7 días?
# Camino D del profesor + valor agregado: establecimientos de salud (RENIPRESS)
#
# Cómo ejecutar:   python incendios_junin.py
# Resultados:      carpeta salidas/ (mapa HTML, mapa PNG y tablas CSV)
# ============================================================

# --- Paquetes ---
# Instalar (una sola vez por computadora):
#   uv pip install geopandas rasterio rasterstats folium scipy requests matplotlib
# Cargar (en cada ejecución): las líneas import de abajo.
import sys
import warnings
from datetime import datetime, timedelta, timezone
from pathlib import Path

import folium                        # folium: arma mapas interactivos (HTML) desde Python
import geopandas as gpd              # geopandas: tablas con geometría (puntos, polígonos); pandas + mapas
import matplotlib.pyplot as plt      # matplotlib: mapa estático (imagen PNG) para el README
import numpy as np                   # numpy: cálculos con arreglos de números (horas, distancias)
import pandas as pd                  # pandas: tablas de filas y columnas, como Excel pero con código
import requests                      # requests: descarga archivos de internet
from rasterstats import zonal_stats  # rasterstats: resume un raster dentro de polígonos (estadística zonal)
from scipy.sparse import coo_matrix  # scipy: guarda qué focos están conectados con cuáles
from scipy.sparse.csgraph import connected_components  # scipy: encuentra los grupos de focos conectados
from scipy.spatial import cKDTree    # scipy: busca rápido los focos vecinos a 1 km o menos

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8")  # para que las tildes se vean bien en la consola de Windows


# ============================================================
# PROCESO 0: Configuración
# Insumo:    nada (son decisiones del análisis)
# Resultado: rutas, direcciones web y parámetros en un solo lugar
# ============================================================

CARPETA = Path(__file__).resolve().parent
DATOS = CARPETA / "datos"          # datos descargados (no se suben a GitHub: pesan mucho)
SALIDAS = CARPETA / "salidas"      # resultados que sí se suben
DATOS.mkdir(exist_ok=True)
SALIDAS.mkdir(exist_ok=True)

# Límites distritales del INEI (Censo 2017), con UBIGEO
SHAPE_DISTRITOS = Path(
    r"C:\Users\Principal\Desktop\1\1_Doc\2_Cursos\2_Ciencia de datos\6_Python basico"
    r"\Trabajo final Python basico\shapefiles\3_Distrital\Distrito_INEI_2017.shp"
)
DEPARTAMENTO = "JUNIN"

# Parámetros del método (los mismos que usa el profesor en su monitor)
RADIO_M = 5_000        # zona de exposición: 5 km alrededor de cada foco
DIST_FOCOS_M = 1_000   # dos focos son del mismo incendio si están a <= 1 km ...
HORAS_FOCOS = 48       # ... y a <= 48 horas entre sí
LADO_PIXEL_M = 375     # tamaño del píxel VIIRS: cada foco representa un cuadrado de 375 m
MIN_FOCOS = 3          # total "estricto": solo incendios con >= 3 focos (descarta detecciones aisladas)
HORAS_RECIENTE = 72    # "actividad reciente": último foco hace 72 h o menos (activo u observación)
CRS_METROS = 32718     # UTM zona 18 Sur: sistema de coordenadas en METROS (Junín cae en esta zona)
NODATA_WORLDPOP = -99999

# Focos de calor: NASA FIRMS, sensor VIIRS (375 m), 3 satélites, últimos 7 días, sin clave
URL_FIRMS = {
    "S-NPP": "https://firms.modaps.eosdis.nasa.gov/data/active_fire/suomi-npp-viirs-c2/csv/SUOMI_VIIRS_C2_South_America_7d.csv",
    "NOAA-20": "https://firms.modaps.eosdis.nasa.gov/data/active_fire/noaa-20-viirs-c2/csv/J1_VIIRS_C2_South_America_7d.csv",
    "NOAA-21": "https://firms.modaps.eosdis.nasa.gov/data/active_fire/noaa-21-viirs-c2/csv/J2_VIIRS_C2_South_America_7d.csv",
}
# Población: WorldPop 2026, raster de 100 m (personas estimadas por celda)
URL_WORLDPOP = ("https://data.worldpop.org/GIS/Population/Global_2015_2030/R2024B/2026/PER/v1/"
                "100m/constrained/per_pop_2026_CN_100m_R2024B_v1.tif")
ARCHIVO_WORLDPOP = DATOS / "per_pop_2026_CN_100m_R2024B_v1.tif"
# Establecimientos de salud: RENIPRESS (SUSALUD), corte del 30-09-2026
URL_RENIPRESS = "https://www.datosabiertos.gob.pe/sites/default/files/RENIPRESS_30-09-2026.csv"
ARCHIVO_RENIPRESS = DATOS / "RENIPRESS_30-09-2026.csv"


def descargar(url, destino):
    """Descarga `url` en `destino` solo si el archivo todavía no existe."""
    if destino.exists():
        print(f"  ya estaba descargado: {destino.name}")
        return destino
    print(f"  descargando {destino.name} ...")
    respuesta = requests.get(url, timeout=900, headers={"User-Agent": "Mozilla/5.0"})
    respuesta.raise_for_status()           # si la descarga falla, se detiene con un error claro
    destino.write_bytes(respuesta.content)
    return destino


# Hora de corte = momento de la ejecución, en hora de Lima (UTC-5; el Perú no cambia de hora)
CORTE = pd.Timestamp(datetime.now(timezone.utc).replace(tzinfo=None) - timedelta(hours=5))
print(f"Corte: {CORTE:%Y-%m-%d %H:%M} (hora de Lima)")


# ============================================================
# PROCESO 1: Cargar los distritos de Junín
# Insumo:    shapefile INEI 2017 (1,874 distritos de todo el Perú)
# Resultado: distritos_junin (solo los distritos de Junín)
# ============================================================

# gpd.read_file(): convierte el shapefile en una tabla con geometría (GeoDataFrame)
#   Antes: archivo en el disco | Después: tabla con una fila por distrito y su polígono
distritos_peru = gpd.read_file(SHAPE_DISTRITOS)

# Filtro con condición + selección de columnas: conserva solo Junín y las columnas útiles
#   Antes: 1,874 filas × 15 columnas | Después: distritos de Junín × 4 columnas
distritos_junin = (distritos_peru.loc[distritos_peru["nombdep"] == DEPARTAMENTO,
                                      ["ubigeo", "nombprov", "nombdist", "geometry"]]
                   .reset_index(drop=True))
print(f"[1] distritos: {distritos_peru.shape} -> Junín {distritos_junin.shape}")


# ============================================================
# PROCESO 2: Descargar los focos de calor (NASA FIRMS)
# Insumo:    3 archivos CSV de FIRMS (uno por satélite, Sudamérica, 7 días)
# Resultado: focos_sudamerica (todos los focos, en una sola tabla)
# ============================================================
# Un foco de calor es un píxel de 375 m donde el satélite detectó una temperatura
# anómala compatible con fuego. NO es un incendio confirmado: también incluye
# quemas agrícolas y de pastos. Los focos se guardan en una carpeta con la fecha,
# así cada actualización no borra las anteriores.

carpeta_firms = DATOS / f"firms_{CORTE:%Y-%m-%d}"
carpeta_firms.mkdir(exist_ok=True)

partes = []
for satelite, url in URL_FIRMS.items():
    archivo = descargar(url, carpeta_firms / f"{satelite}_7d.csv")
    # pd.read_csv(): lee el CSV como tabla; acq_time se lee como texto para no perder los ceros ("0525")
    tabla = pd.read_csv(archivo, dtype={"acq_time": str})
    tabla["satelite"] = satelite
    partes.append(tabla)

# pd.concat(): apila las 3 tablas una debajo de otra (como copiar y pegar filas en Excel)
#   Antes: 3 tablas | Después: 1 tabla con todas las filas
focos_sudamerica = pd.concat(partes, ignore_index=True)
print(f"[2] focos en Sudamérica: {focos_sudamerica.shape}")


# ============================================================
# PROCESO 3: Limpiar los focos y quedarse con los de Junín
# Insumo:    focos_sudamerica + distritos_junin
# Resultado: focos_junin (focos dentro de Junín, con hora de Lima y distrito)
# ============================================================

# --- Paso 3.1: Hora de Lima ---
# FIRMS entrega la hora en UTC (hora universal) y como texto HHMM.
# pd.to_datetime(): junta fecha + hora de texto y la convierte en fecha-hora real
#   Antes: "2026-10-16" y "0525" | Después: 2026-10-16 05:25 (UTC)
focos_sudamerica["fecha_hora_utc"] = pd.to_datetime(
    focos_sudamerica["acq_date"] + " " + focos_sudamerica["acq_time"].str.zfill(4),
    format="%Y-%m-%d %H%M")
# Restar 5 horas: UTC -> hora de Lima
focos_sudamerica["fecha_hora_lima"] = focos_sudamerica["fecha_hora_utc"] - pd.Timedelta(hours=5)

# --- Paso 3.2: Prefiltro rápido con la caja que encierra a Junín ---
# .total_bounds: devuelve [lon mín, lat mín, lon máx, lat máx] de Junín; no modifica nada
lon_min, lat_min, lon_max, lat_max = distritos_junin.total_bounds
en_caja = focos_sudamerica[focos_sudamerica["longitude"].between(lon_min, lon_max)
                           & focos_sudamerica["latitude"].between(lat_min, lat_max)]

# --- Paso 3.3: Convertir la tabla en puntos ---
# gpd.points_from_xy(): crea un punto por fila a partir de longitud y latitud
#   Antes: dos columnas de números | Después: una columna "geometry" con puntos
puntos = gpd.GeoDataFrame(en_caja,
                          geometry=gpd.points_from_xy(en_caja["longitude"], en_caja["latitude"]),
                          crs="EPSG:4326")  # EPSG:4326 = coordenadas en grados (lat/lon)

# --- Paso 3.4: Unión espacial con los distritos ---
# gpd.sjoin(): a cada foco le pega el distrito en el que cae (como un BUSCARV, pero por ubicación)
#   Antes: focos sin distrito | Después: solo los focos dentro de Junín, con ubigeo y nombre
focos_junin = (gpd.sjoin(puntos, distritos_junin, predicate="within", how="inner")
               .drop(columns="index_right")
               .reset_index(drop=True))

# Marca (no elimina) los focos de confianza baja, para verlos distintos en el mapa
focos_junin["confianza_baja"] = focos_junin["confidence"].eq("low")
# Antigüedad del foco respecto al corte, en horas
focos_junin["horas_desde_foco"] = (CORTE - focos_junin["fecha_hora_lima"]).dt.total_seconds() / 3600

print(f"[3] en la caja de Junín: {len(en_caja):,} | dentro de distritos de Junín: {len(focos_junin):,} "
      f"| confianza baja: {focos_junin['confianza_baja'].sum():,}")
print(f"    periodo (hora de Lima): {focos_junin['fecha_hora_lima'].min():%Y-%m-%d %H:%M} "
      f"a {focos_junin['fecha_hora_lima'].max():%Y-%m-%d %H:%M}")


# ============================================================
# PROCESO 4: Agrupar los focos en incendios
# Insumo:    focos_junin
# Resultado: incendios (una fila por incendio) y la columna "incendio" en cada foco
# ============================================================
# Regla del profesor: dos focos son del mismo incendio si están a <= 1 km y a <= 48 h,
# directamente o a través de otros focos (como una cadena).

# --- Paso 4.1: Pasar a metros ---
# .to_crs(): cambia el sistema de coordenadas; en grados no se pueden medir distancias
#   Antes: lon/lat en grados | Después: x/y en metros (UTM 18 Sur)
focos_m = focos_junin.to_crs(epsg=CRS_METROS)
xy = np.column_stack([focos_m.geometry.x, focos_m.geometry.y])

# --- Paso 4.2: Pares de focos cercanos en el espacio ---
# cKDTree(...).query_pairs(): lista todos los pares de focos a <= 1,000 m
pares = cKDTree(xy).query_pairs(r=DIST_FOCOS_M, output_type="ndarray")

# --- Paso 4.3: De esos pares, quedarse con los cercanos en el tiempo (<= 48 h) ---
horas = (focos_junin["fecha_hora_lima"] - focos_junin["fecha_hora_lima"].min()).dt.total_seconds().to_numpy() / 3600
pares = pares[np.abs(horas[pares[:, 0]] - horas[pares[:, 1]]) <= HORAS_FOCOS]

# --- Paso 4.4: Grupos conectados = incendios ---
# connected_components(): recorre la "red" de pares y le da el mismo número a los focos encadenados
#   Antes: lista de pares | Después: un número de grupo por foco
n_focos = len(focos_junin)
red = coo_matrix((np.ones(len(pares)), (pares[:, 0], pares[:, 1])), shape=(n_focos, n_focos))
n_incendios, grupo = connected_components(red, directed=False)

# Nombre legible: INC-001 es el incendio con más focos
orden = pd.Series(grupo).value_counts().index
codigo = {g: f"INC-{i + 1:03d}" for i, g in enumerate(orden)}
focos_junin["incendio"] = [codigo[g] for g in grupo]
focos_m["incendio"] = focos_junin["incendio"].values

# --- Paso 4.5: Zona de focos (huella) ---
# .buffer(..., cap_style="square"): dibuja el cuadrado de 375 m del píxel alrededor de cada foco
# .dissolve(): une los cuadrados del mismo incendio en un solo polígono (como "agrupar" en una tabla dinámica)
#   Antes: un punto por foco | Después: un polígono por incendio
cuadrados = focos_m[["incendio", "geometry"]].copy()
cuadrados["geometry"] = focos_m.buffer(LADO_PIXEL_M / 2, cap_style="square")
huella_m = cuadrados.dissolve(by="incendio")
huella_m["huella_ha"] = huella_m.area / 10_000   # m² -> hectáreas

# --- Paso 4.6: Tabla de incendios ---
# .groupby().agg(): resume los focos de cada incendio (como una tabla dinámica)
#   Antes: una fila por foco | Después: una fila por incendio
def mas_frecuente(serie):
    """Devuelve el valor que más se repite (el distrito donde cae la mayoría de focos)."""
    return serie.mode().iat[0]

incendios = (focos_junin.groupby("incendio")
             .agg(focos=("incendio", "size"),
                  inicio=("fecha_hora_lima", "min"),
                  ultimo_foco=("fecha_hora_lima", "max"),
                  frp_total_mw=("frp", "sum"),   # FRP: potencia radiativa del fuego, en megavatios
                  lat=("latitude", "mean"),
                  lon=("longitude", "mean"),
                  distrito=("nombdist", mas_frecuente),
                  provincia=("nombprov", mas_frecuente),
                  ubigeo=("ubigeo", mas_frecuente),
                  distritos_con_focos=("nombdist", lambda s: ", ".join(sorted(s.unique())))))
incendios["huella_ha"] = huella_m["huella_ha"]
# np.select(): asigna el estado según las horas sin focos (como un SI anidado de Excel)
horas_sin_focos = (CORTE - incendios["ultimo_foco"]).dt.total_seconds() / 3600
incendios["estado"] = np.select([horas_sin_focos < 24, horas_sin_focos < 72],
                                ["Activo", "En observación"], default="Sin actividad reciente")
print(f"[4] {n_focos:,} focos -> {n_incendios} incendios | estados: {incendios['estado'].value_counts().to_dict()}")


# ============================================================
# PROCESO 5: Construir la zona de 5 km
# Insumo:    focos_m (focos en metros)
# Resultado: zonas_5km (un polígono por incendio) y zona_total (todas unidas)
# ============================================================

# .buffer(5000): dibuja un círculo de 5 km alrededor de cada foco
# .dissolve(by="incendio"): une los círculos de un mismo incendio
#   Antes: un punto por foco | Después: un polígono de 5 km por incendio
circulos = focos_m[["incendio", "geometry"]].copy()
circulos["geometry"] = focos_m.buffer(RADIO_M)
zonas_5km_m = circulos.dissolve(by="incendio")

# .union_all(): junta todas las zonas en una sola. Así, si dos incendios se superponen,
# las personas de la parte común se cuentan UNA sola vez en el total.
zona_total_m = gpd.GeoDataFrame(geometry=[zonas_5km_m.union_all()], crs=CRS_METROS)

# El raster de población está en grados (EPSG:4326): las zonas deben estar en el mismo sistema
zonas_5km = zonas_5km_m.to_crs(epsg=4326)
zona_total = zona_total_m.to_crs(epsg=4326)
huella = huella_m.to_crs(epsg=4326)
print(f"[5] zona total de 5 km: {zona_total_m.area.iat[0] / 1e6:,.0f} km²")


# ============================================================
# PROCESO 6: Contar la población con WorldPop (dato raster)
# Insumo:    raster WorldPop 2026 (100 m) + zonas_5km + zona_total + distritos_junin
# Resultado: población a 5 km por incendio, en total y por distrito
# ============================================================

descargar(URL_WORLDPOP, ARCHIVO_WORLDPOP)


def poblacion(poligonos):
    """Suma las personas de WorldPop dentro de cada polígono (estadística zonal, como en el notebook)."""
    # zonal_stats(): recorta el raster con cada polígono y suma sus celdas
    #   Antes: polígonos sin población | Después: una suma de personas por polígono
    #   nodata explícito: así las celdas "sin dato" (-99999) no entran a la suma
    resultado = zonal_stats(poligonos, str(ARCHIVO_WORLDPOP), stats=["sum"], nodata=NODATA_WORLDPOP)
    return [round(r["sum"] or 0) for r in resultado]


# --- Paso 6.1: Por incendio (las zonas pueden superponerse: NO sumar esta columna) ---
zonas_5km["poblacion_5km"] = poblacion(zonas_5km)
incendios["poblacion_5km"] = zonas_5km["poblacion_5km"]

# --- Paso 6.2: Total sin doble conteo ---
poblacion_total_5km = poblacion(zona_total)[0]

# --- Paso 6.2b: Total estricto: solo incendios con 3 o más focos ---
# Un foco aislado en una ciudad puede ser una fuente de calor urbana, no un incendio forestal.
# .loc[...]: elige solo las zonas de esos incendios; .union_all() las junta sin doble conteo
#   Antes: zonas de todos los incendios | Después: una zona con los incendios de >= 3 focos
incendios_3focos = incendios.index[incendios["focos"] >= MIN_FOCOS]
zona_3focos = gpd.GeoDataFrame(geometry=[zonas_5km_m.loc[incendios_3focos].union_all()],
                               crs=CRS_METROS).to_crs(epsg=4326)
poblacion_3focos_5km = poblacion(zona_3focos)[0]
print(f"[6] incendios con >= {MIN_FOCOS} focos: {len(incendios_3focos)} "
      f"| población a <= 5 km de ellos: {poblacion_3focos_5km:,}")

# --- Paso 6.2c: Total para decidir hoy: incendios con actividad reciente (<= 72 h) ---
# Son los incendios "Activos" (último foco < 24 h) o "En observación" (24–72 h)
incendios_recientes = incendios.index[incendios["estado"] != "Sin actividad reciente"]
zona_recientes = gpd.GeoDataFrame(geometry=[zonas_5km_m.loc[incendios_recientes].union_all()],
                                  crs=CRS_METROS).to_crs(epsg=4326)
poblacion_recientes_5km = poblacion(zona_recientes)[0]
print(f"[6] incendios con actividad en <= {HORAS_RECIENTE} h: {len(incendios_recientes)} "
      f"| población a <= 5 km de ellos: {poblacion_recientes_5km:,}")

# --- Paso 6.3: Por distrito de Junín ---
# gpd.overlay(..., how="intersection"): recorta la zona de 5 km con cada distrito
#   Antes: una zona total | Después: un pedazo de zona por distrito
pedazos = gpd.overlay(distritos_junin, zona_total, how="intersection")
pedazos["poblacion_expuesta"] = poblacion(pedazos)
distritos_junin["poblacion_total"] = poblacion(distritos_junin)
distritos_junin = distritos_junin.merge(pedazos[["ubigeo", "poblacion_expuesta"]], on="ubigeo", how="left")
distritos_junin["poblacion_expuesta"] = distritos_junin["poblacion_expuesta"].fillna(0).astype(int)
distritos_junin["pct_expuesta"] = (100 * distritos_junin["poblacion_expuesta"]
                                   / distritos_junin["poblacion_total"].replace(0, np.nan)).round(1)
print(f"[6] población a <= 5 km (sin doble conteo): {poblacion_total_5km:,} "
      f"| solo dentro de Junín: {distritos_junin['poblacion_expuesta'].sum():,} "
      f"| población de Junín: {distritos_junin['poblacion_total'].sum():,}")


# ============================================================
# PROCESO 7: Establecimientos de salud a 5 km o menos (RENIPRESS)
# Insumo:    CSV RENIPRESS (SUSALUD) + focos_m
# Resultado: establecimientos_expuestos (con distancia al foco más cercano)
# ============================================================

descargar(URL_RENIPRESS, ARCHIVO_RENIPRESS)

# El archivo usa ";" como separador y UTF-8 con BOM (por eso "utf-8-sig")
renipress = pd.read_csv(ARCHIVO_RENIPRESS, sep=";", encoding="utf-8-sig", dtype=str)
# pd.to_numeric(): convierte las coordenadas de texto a número; si no se puede, deja vacío (NaN)
renipress["lat"] = pd.to_numeric(renipress["NORTE"], errors="coerce")
renipress["lon"] = pd.to_numeric(renipress["ESTE"], errors="coerce")

# --- Paso 7.1: Filtros en cascada (se informa cuántos quedan en cada paso) ---
filtros = [
    ("registros en RENIPRESS", pd.Series(True, index=renipress.index)),
    ("estado ACTIVO en el registro", renipress["ESTADO"].eq("ACTIVO")),
    ("establecimiento de salud (no servicio de apoyo)", renipress["TIPO_ESTABLECIMIENTO"].str.startswith("ESTABLECIMIENTO")),
    ("con categoría asignada (no '0')", renipress["CATEGORIA"].ne("0")),
    ("con coordenadas dentro del Perú", renipress["lat"].between(-18.4, -0.01) & renipress["lon"].between(-81.4, -68.6)),
]
conserva = pd.Series(True, index=renipress.index)
print("[7] filtros de RENIPRESS (todo el Perú):")
for descripcion, condicion in filtros:
    conserva &= condicion
    print(f"    {descripcion:50s} {conserva.sum():>7,}")
establecimientos = renipress[conserva].copy()
establecimientos = gpd.GeoDataFrame(establecimientos,
                                    geometry=gpd.points_from_xy(establecimientos["lon"], establecimientos["lat"]),
                                    crs="EPSG:4326")
establecimientos["red"] = np.where(establecimientos["INSTITUCION"].eq("PRIVADO"), "Privada", "Pública")

# --- Paso 7.2: Distancia al foco más cercano (hasta 5 km) ---
# gpd.sjoin_nearest(): a cada establecimiento le pega el foco más cercano a <= 5 km y la distancia
#   Antes: establecimientos sin distancia | Después: solo los que tienen un foco a <= 5 km
expuestos_m = gpd.sjoin_nearest(establecimientos.to_crs(epsg=CRS_METROS),
                                focos_m[["incendio", "geometry"]],
                                how="inner", max_distance=RADIO_M, distance_col="dist_foco_m")
# Si hay empate entre dos focos, el establecimiento sale dos veces: se deja una sola fila
expuestos_m = expuestos_m.drop_duplicates(subset="COD_IPRESS")
establecimientos_expuestos = expuestos_m.to_crs(epsg=4326)
establecimientos_expuestos["dist_foco_km"] = (establecimientos_expuestos["dist_foco_m"] / 1000).round(2)

# Conteo de establecimientos por incendio
incendios["establecimientos_5km"] = (establecimientos_expuestos.groupby("incendio").size()
                                     .reindex(incendios.index).fillna(0).astype(int))
# Establecimientos de Junín (para mostrarlos todos en el mapa)
establecimientos_junin = gpd.sjoin(establecimientos, distritos_junin[["ubigeo", "geometry"]],
                                   predicate="within", how="inner").drop(columns=["index_right", "ubigeo"])
# Establecimientos dentro de la zona estricta (incendios con >= 3 focos)
salud_3focos = gpd.sjoin(establecimientos_expuestos.drop(columns="index_right"), zona_3focos,
                         predicate="within", how="inner")
# Establecimientos cerca de incendios con actividad reciente (<= 72 h)
salud_recientes = gpd.sjoin(establecimientos_expuestos.drop(columns="index_right"), zona_recientes,
                            predicate="within", how="inner")
print(f"    establecimientos a <= 5 km de incendios con actividad en <= {HORAS_RECIENTE} h: "
      f"{len(salud_recientes):,}")
print(f"    establecimientos en Junín: {len(establecimientos_junin):,} "
      f"| a <= 5 km de un foco: {len(establecimientos_expuestos):,} "
      f"(públicos: {establecimientos_expuestos['red'].eq('Pública').sum()}) "
      f"| a <= 5 km de incendios con >= {MIN_FOCOS} focos: {len(salud_3focos):,}")


# ============================================================
# PROCESO 8: Guardar las tablas de resultados
# Insumo:    incendios, distritos_junin, establecimientos_expuestos
# Resultado: 3 archivos CSV en salidas/
# ============================================================

incendios = incendios.sort_values("poblacion_5km", ascending=False)
tabla_incendios = incendios.reset_index()[[
    "incendio", "estado", "provincia", "distrito", "ubigeo", "distritos_con_focos", "focos",
    "inicio", "ultimo_foco", "huella_ha", "frp_total_mw", "poblacion_5km", "establecimientos_5km", "lat", "lon"]]
tabla_incendios = tabla_incendios.round({"huella_ha": 0, "frp_total_mw": 1, "lat": 4, "lon": 4})
tabla_incendios.to_csv(SALIDAS / "incendios_junin.csv", index=False, encoding="utf-8-sig")

tabla_distritos = (distritos_junin.drop(columns="geometry")
                   .query("poblacion_expuesta > 0")
                   .sort_values("poblacion_expuesta", ascending=False))
tabla_distritos.to_csv(SALIDAS / "distritos_poblacion_expuesta.csv", index=False, encoding="utf-8-sig")

columnas_salud = ["COD_IPRESS", "NOMBRE", "INSTITUCION", "red", "CATEGORIA", "TIPO_ESTABLECIMIENTO",
                  "DEPARTAMENTO", "PROVINCIA", "DISTRITO", "UBIGEO", "RED", "MICRORED", "ESTADO",
                  "incendio", "dist_foco_km", "lat", "lon"]
tabla_salud = establecimientos_expuestos[columnas_salud].sort_values("dist_foco_km")
tabla_salud.to_csv(SALIDAS / "establecimientos_expuestos.csv", index=False, encoding="utf-8-sig")
print(f"[8] tablas guardadas en {SALIDAS}")
print(tabla_incendios.head(8)[["incendio", "estado", "distrito", "focos", "huella_ha",
                               "poblacion_5km", "establecimientos_5km"]].to_string(index=False))


# ============================================================
# PROCESO 9: Ficha territorial para el Gobierno Regional
# Insumo:    focos_junin, incendios, distritos_junin, establecimientos_expuestos
# Resultado: ficha_distritos_junin.csv (prioridad por distrito) y salud_por_red_junin.csv
# ============================================================
# Enfoque territorial: el Gobierno Regional coordina con las municipalidades (distritos y
# provincias) y con su red de salud (DIRESA -> redes -> microrredes). Por eso los
# resultados se resumen por esas unidades, y no solo por incendio.

# --- Paso 9.1: Focos por distrito (7 días y últimas 72 h) ---
# .groupby("ubigeo").agg(): cuenta los focos de cada distrito (como una tabla dinámica)
#   Antes: una fila por foco | Después: una fila por distrito con dos conteos
focos_junin["reciente"] = focos_junin["horas_desde_foco"] <= HORAS_RECIENTE
focos_por_distrito = focos_junin.groupby("ubigeo").agg(focos_7d=("ubigeo", "size"),
                                                       focos_72h=("reciente", "sum"))

# --- Paso 9.2: Incendios prioritarios por distrito ---
# Prioritario = incendio con >= 3 focos y con actividad en las últimas 72 h
prioritarios = incendios.index[(incendios["focos"] >= MIN_FOCOS)
                               & (incendios["estado"] != "Sin actividad reciente")]
incendios_por_distrito = (focos_junin[focos_junin["incendio"].isin(prioritarios)]
                          .groupby("ubigeo")["incendio"].nunique().rename("incendios_prioritarios"))

# --- Paso 9.3: Establecimientos de salud expuestos ubicados en cada distrito ---
salud_por_distrito = establecimientos_expuestos.groupby("UBIGEO").size().rename("establecimientos_expuestos")

# --- Paso 9.4: Unir todo en una ficha por distrito ---
# .join(): pega columnas usando el UBIGEO como llave (como BUSCARV); sin dato -> 0
#   Antes: 4 tablas separadas | Después: una fila por distrito con todos los indicadores
ficha = (distritos_junin.drop(columns="geometry").set_index("ubigeo")
         .join([focos_por_distrito, incendios_por_distrito, salud_por_distrito]))
conteos = ["focos_7d", "focos_72h", "incendios_prioritarios", "establecimientos_expuestos"]
ficha[conteos] = ficha[conteos].fillna(0).astype(int)

# --- Paso 9.5: Semáforo de prioridad (regla simple y explicable) ---
# np.select(): como un SI anidado de Excel; gana la primera condición que se cumple
ficha["prioridad"] = np.select(
    [ficha["incendios_prioritarios"] > 0,   # Alta: incendio de >= 3 focos con actividad en 72 h
     ficha["focos_72h"] > 0,                # Media: focos recientes, pero aislados
     ficha["focos_7d"] > 0],                # Baja: solo focos de hace más de 72 h
    ["1 · Alta", "2 · Media", "3 · Baja"], default="Sin focos")
ficha = ficha.sort_values(["prioridad", "poblacion_expuesta"], ascending=[True, False])
ficha.reset_index().to_csv(SALIDAS / "ficha_distritos_junin.csv", index=False, encoding="utf-8-sig")

# .map(): copia la prioridad y los conteos a la tabla con geometría (para colorear el mapa)
for col in ["prioridad"] + conteos:
    distritos_junin[col] = distritos_junin["ubigeo"].map(ficha[col])
print(f"[9] prioridad por distrito: {ficha['prioridad'].value_counts().sort_index().to_dict()}")

# --- Paso 9.6: Establecimientos expuestos de Junín por red de salud ---
# Junín = UBIGEO que empieza con "12". Los privados no pertenecen a una red de la DIRESA.
salud_junin = establecimientos_expuestos[establecimientos_expuestos["UBIGEO"].str.startswith("12")].copy()
salud_junin["RED"] = salud_junin["RED"].fillna("Sin red (privado u otro)")
salud_junin["cerca_de_incendio_reciente"] = salud_junin["COD_IPRESS"].isin(salud_recientes["COD_IPRESS"])
salud_por_red = (salud_junin.groupby("RED")
                 .agg(establecimientos_a_5km=("COD_IPRESS", "size"),
                      cerca_de_incendios_recientes=("cerca_de_incendio_reciente", "sum"),
                      distancia_minima_km=("dist_foco_km", "min"))
                 .sort_values("cerca_de_incendios_recientes", ascending=False))
salud_por_red.to_csv(SALIDAS / "salud_por_red_junin.csv", encoding="utf-8-sig")
print(salud_por_red.head(6).to_string())


# ============================================================
# PROCESO 10: Mapa interactivo (Folium)
# Insumo:    todas las capas anteriores + ficha territorial (Proceso 9)
# Resultado: salidas/mapa_incendios_junin.html (capas que se activan y desactivan)
# ============================================================

def a_texto(gdf):
    """Convierte las fechas a texto: el formato del mapa web (GeoJSON) no acepta fechas."""
    gdf = gdf.copy()
    for col in gdf.columns:
        if pd.api.types.is_datetime64_any_dtype(gdf[col]):
            gdf[col] = gdf[col].dt.strftime("%Y-%m-%d %H:%M")
    return gdf


# folium.Map(): crea el mapa base centrado en Junín
lat_centro = (lat_min + lat_max) / 2
lon_centro = (lon_min + lon_max) / 2
mapa = folium.Map(location=[lat_centro, lon_centro], zoom_start=8, tiles="CartoDB positron",
                  control_scale=True)
folium.TileLayer("OpenStreetMap", name="Calles (OpenStreetMap)").add_to(mapa)
folium.TileLayer(tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
                 attr="Esri World Imagery", name="Imagen satelital (Esri)").add_to(mapa)

# Paleta: cada tema tiene su propio color para que no se confundan
#   morado = territorio (distritos, provincias, prioridad) · rojo/naranja/amarillo = fuego
#   azul = salud · gris punteado = área de exposición de 5 km
COLOR_PRIORIDAD = {"1 · Alta": "#54278f", "2 · Media": "#9e9ac8", "3 · Baja": "#dadaeb", "Sin focos": "#ffffff"}

# --- Capa 1: Distritos, coloreados según la prioridad territorial ---
# .simplify(): quita vértices para que el mapa web pese menos (tolerancia ~100 m)
distritos_web = distritos_junin.copy()
distritos_web["geometry"] = distritos_web.simplify(0.001)
folium.GeoJson(
    distritos_web, name="Distritos · prioridad territorial",
    style_function=lambda f: {"color": "#3f007d", "weight": 1.2,
                              "fillColor": COLOR_PRIORIDAD[f["properties"]["prioridad"]],
                              "fillOpacity": 0 if f["properties"]["prioridad"] == "Sin focos" else 0.4},
    highlight_function=lambda f: {"weight": 3, "color": "#3f007d"},
    tooltip=folium.GeoJsonTooltip(
        fields=["nombdist", "nombprov", "ubigeo", "prioridad", "focos_72h", "poblacion_expuesta",
                "establecimientos_expuestos", "poblacion_total"],
        aliases=["Distrito", "Provincia", "UBIGEO", "Prioridad", "Focos en las últimas 72 h",
                 "Personas a ≤ 5 km de focos", "Establecimientos de salud a ≤ 5 km",
                 "Población 2026 (WorldPop)"], localize=True),
).add_to(mapa)

# --- Capa 2: Límites provinciales (línea gruesa, sin relleno) ---
# .dissolve(by="nombprov"): une los distritos de cada provincia en un solo polígono
#   Antes: 124 distritos | Después: 9 provincias
provincias = distritos_junin[["nombprov", "geometry"]].dissolve(by="nombprov").reset_index()
provincias["geometry"] = provincias.simplify(0.001)
folium.GeoJson(provincias, name="Provincias (límite)",
               style_function=lambda f: {"color": "#3f007d", "weight": 3.5, "fill": False}).add_to(mapa)

# --- Capa 3: Área de exposición de 5 km (NO es el avance del fuego) ---
zonas_web = zonas_5km.reset_index()[["incendio", "poblacion_5km", "geometry"]]
zonas_web["geometry"] = zonas_web.simplify(0.0005)
folium.GeoJson(
    zonas_web, name="Área de exposición: 5 km alrededor de los focos",
    style_function=lambda f: {"color": "#4d4d4d", "weight": 1, "dashArray": "5,4",
                              "fillColor": "#bdbdbd", "fillOpacity": 0.15},
    tooltip=folium.GeoJsonTooltip(fields=["incendio", "poblacion_5km"],
                                  aliases=["Área de exposición del incendio", "Personas a ≤ 5 km"],
                                  localize=True),
).add_to(mapa)

# --- Capa 4: Zona de focos (huella de píxeles; NO es un perímetro confirmado) ---
huella_web = huella.reset_index()[["incendio", "huella_ha", "geometry"]].round({"huella_ha": 0})
folium.GeoJson(
    huella_web, name="Zona de focos (estimada, no es perímetro)",
    style_function=lambda f: {"color": "#8b0000", "weight": 0.6, "fillColor": "#d7301f", "fillOpacity": 0.55},
    tooltip=folium.GeoJsonTooltip(fields=["incendio", "huella_ha"],
                                  aliases=["Incendio", "Huella VIIRS (ha)"], localize=True),
).add_to(mapa)

# --- Capa 5: Focos de calor (detecciones satelitales) ---
capa_focos = folium.FeatureGroup(name="Focos de calor VIIRS (detección satelital)")
for _, foco in a_texto(focos_junin).iterrows():
    # color según antigüedad: rojo < 24 h, naranja 24–72 h, amarillo > 72 h
    color = "#d7191c" if foco["horas_desde_foco"] < 24 else "#fd8d3c" if foco["horas_desde_foco"] < 72 else "#fecc5c"
    folium.CircleMarker(
        location=[foco["latitude"], foco["longitude"]], radius=3, color=color, weight=0.5,
        fill=True, fill_color=color,
        fill_opacity=0.35 if foco["confianza_baja"] else 0.9,   # los de confianza baja, más transparentes
        tooltip=(f"Foco {foco['fecha_hora_lima']} (Lima) · {foco['satelite']} · FRP {foco['frp']:.1f} MW"
                 f" · confianza {foco['confidence']} · {foco['incendio']}"),
    ).add_to(capa_focos)
capa_focos.add_to(mapa)

# --- Capa 6: Incendios con >= 3 focos (ícono de fuego; clic = resumen) ---
# folium.DivIcon(): un ícono hecho con HTML; aquí, un emoji de fuego más grande si es más reciente
capa_incendios = folium.FeatureGroup(name=f"Incendios con ≥ {MIN_FOCOS} focos (clic: resumen)")
for _, inc in a_texto(tabla_incendios[tabla_incendios["focos"] >= MIN_FOCOS]).iterrows():
    tam = {"Activo": 26, "En observación": 20}.get(inc["estado"], 14)
    resumen = (f"<b>{inc['incendio']}</b> · {inc['estado']}<br>"
               f"{inc['distrito']} ({inc['provincia']}) · UBIGEO {inc['ubigeo']}<br>"
               f"Focos: {inc['focos']} · Huella VIIRS: {inc['huella_ha']:,.0f} ha<br>"
               f"Del {inc['inicio']} al {inc['ultimo_foco']}<br>"
               f"<b>Personas a ≤ 5 km: {inc['poblacion_5km']:,}</b><br>"
               f"Establecimientos de salud a ≤ 5 km: {inc['establecimientos_5km']}")
    folium.Marker(location=[inc["lat"], inc["lon"]],
                  icon=folium.DivIcon(html=f'<div style="font-size:{tam}px;line-height:1">🔥</div>',
                                      icon_size=(tam, tam), icon_anchor=(tam // 2, tam // 2)),
                  tooltip=f"{inc['incendio']} · {inc['estado']} · {inc['distrito']}",
                  popup=folium.Popup(resumen, max_width=320)).add_to(capa_incendios)
capa_incendios.add_to(mapa)


def icono_salud(publico):
    """Cuadrado con una cruz: azul relleno = pública; blanco con borde azul = privada."""
    fondo, cruz = ("#1565c0", "white") if publico else ("white", "#1565c0")
    return folium.DivIcon(
        html=(f'<div style="background:{fondo};color:{cruz};border:2px solid #1565c0;border-radius:3px;'
              f'width:14px;height:14px;font:bold 14px/14px Arial;text-align:center;'
              f'box-shadow:0 0 3px white">+</div>'),
        icon_size=(18, 18), icon_anchor=(9, 9))


# --- Capa 7: Establecimientos de salud a ≤ 5 km de un foco ---
capa_salud_exp = folium.FeatureGroup(name="Salud: establecimientos a ≤ 5 km de un foco")
for _, es in tabla_salud.iterrows():
    folium.Marker(
        location=[es["lat"], es["lon"]], icon=icono_salud(es["red"] == "Pública"),
        tooltip=(f"{es['NOMBRE']} · categoría {es['CATEGORIA']} · {es['INSTITUCION']}"
                 f" · red {es['RED']} · a {es['dist_foco_km']} km del foco más cercano"),
    ).add_to(capa_salud_exp)
capa_salud_exp.add_to(mapa)

# --- Capa 8: Todos los establecimientos de salud de Junín (apagada al inicio) ---
capa_salud = folium.FeatureGroup(name="Salud: todos los establecimientos de Junín", show=False)
for _, es in establecimientos_junin.iterrows():
    folium.CircleMarker(location=[es["lat"], es["lon"]], radius=2.5, color="#1565c0", weight=0.5,
                        fill=True, fill_opacity=0.8,
                        tooltip=f"{es['NOMBRE']} · {es['CATEGORIA']} · {es['INSTITUCION']} · {es['ESTADO']}",
                        ).add_to(capa_salud)
capa_salud.add_to(mapa)

# folium.LayerControl(): el panel para activar y desactivar capas
folium.LayerControl(collapsed=False).add_to(mapa)

# --- Panel "Para decidir hoy" (lo primero que ve el gobernador) ---
alta = ficha[ficha["prioridad"] == "1 · Alta"].head(10)
filas = "".join(
    f"<tr><td>{d.nombdist.title()}</td><td>{d.nombprov.title()}</td>"
    f"<td align='right'>{d.focos_72h}</td><td align='right'>{d.poblacion_expuesta:,}</td>"
    f"<td align='right'>{d.establecimientos_expuestos}</td></tr>"
    for d in alta.itertuples())
# Solo redes de la DIRESA (se excluyen los privados, que no pertenecen a una red)
redes_diresa = salud_por_red[~salud_por_red.index.str.contains("Sin red|NO PERTENECE", case=False)]
redes = ", ".join(f"{red.title()} ({n})" for red, n in
                  redes_diresa["cerca_de_incendios_recientes"].head(4).items() if n > 0)
panel = f"""
<details open style="position:fixed; top:80px; left:12px; z-index:9999; background:white; padding:8px 12px;
         border:2px solid #3f007d; border-radius:6px; font:12px/1.4 Arial; max-width:390px; max-height:55vh; overflow:auto;">
<summary style="cursor:pointer; font-weight:bold; color:#3f007d; font-size:13px">
Para decidir hoy · Gobierno Regional de Junín</summary>
<div style="margin-top:6px">
<b>{len(incendios_recientes)}</b> incendios con actividad en las últimas 72 h ·
<b>{poblacion_recientes_5km:,}</b> personas a ≤ 5 km ·
<b>{len(salud_recientes)}</b> establecimientos de salud a ≤ 5 km<br><br>
<b>Distritos de prioridad alta</b> (incendio de ≥ {MIN_FOCOS} focos con actividad en 72 h):
<table style="border-collapse:collapse; width:100%; margin-top:4px">
<tr style="border-bottom:1px solid #999"><th align="left">Distrito</th><th align="left">Provincia</th>
<th>Focos 72 h</th><th>Personas ≤ 5 km</th><th>Salud ≤ 5 km</th></tr>
{filas}
</table><br>
<b>Redes de salud a alertar</b> (establecimientos cerca de incendios recientes): {redes or "ninguna"}<br>
<i>Uso: coordinar con las municipalidades de estos distritos y alertar a las redes de salud.
Detalle: salidas/ficha_distritos_junin.csv y salud_por_red_junin.csv</i>
</div></details>"""
mapa.get_root().html.add_child(folium.Element(panel))

# --- Leyenda: cómo leer el mapa ---
leyenda = f"""
<details open style="position:fixed; bottom:24px; left:12px; z-index:9999; background:white; padding:8px 12px;
         border:1px solid #999; border-radius:6px; font:12px/1.55 Arial; max-width:340px;">
<summary style="cursor:pointer; font-weight:bold">Cómo leer el mapa</summary>
<b>Territorio</b> (morado): <span style="border-bottom:3px solid #3f007d">&nbsp;provincia&nbsp;</span>
<span style="border-bottom:1px solid #3f007d">&nbsp;distrito&nbsp;</span> · prioridad:
<span style="background:#54278f;opacity:.6">&nbsp;&nbsp;&nbsp;</span> alta
<span style="background:#9e9ac8;opacity:.6">&nbsp;&nbsp;&nbsp;</span> media
<span style="background:#dadaeb">&nbsp;&nbsp;&nbsp;</span> baja<br>
<b>Fuego</b>: 🔥 incendio de ≥ {MIN_FOCOS} focos (más grande = más reciente) ·
<span style="color:#d7191c">●</span> foco &lt; 24 h <span style="color:#fd8d3c">●</span> 24–72 h
<span style="color:#fecc5c">●</span> &gt; 72 h ·
<span style="background:#d7301f;opacity:.7">&nbsp;&nbsp;&nbsp;</span> zona de focos (píxeles de 375 m; <b>no es perímetro</b>)<br>
<b>Exposición</b>: <span style="border:1.5px dashed #4d4d4d">&nbsp;&nbsp;&nbsp;&nbsp;</span>
5 km alrededor de los focos. <b>No es el avance del fuego</b>: es el radio donde se cuenta a la población.<br>
<b>Salud</b>: <span style="background:#1565c0;color:white;font-weight:bold;padding:0 3px">+</span> pública
<span style="border:1.5px solid #1565c0;color:#1565c0;font-weight:bold;padding:0 3px">+</span> privada (a ≤ 5 km de un foco)<br>
Personas a ≤ 5 km: <b>{poblacion_total_5km:,}</b> de cualquier foco · <b>{poblacion_3focos_5km:,}</b> de incendios
con ≥ {MIN_FOCOS} focos (WorldPop 2026, sin doble conteo)<br>
<i>Un foco de calor no es un incendio confirmado; personas cercanas no son damnificadas.</i><br>
Corte {CORTE:%d-%m-%Y %H:%M} (Lima) · Fuentes: NASA FIRMS · WorldPop · RENIPRESS (30-09-2026) · INEI
</details>"""
mapa.get_root().html.add_child(folium.Element(leyenda))

mapa.save(SALIDAS / "mapa_incendios_junin.html")
print(f"[10] mapa interactivo: {SALIDAS / 'mapa_incendios_junin.html'}")


# ============================================================
# PROCESO 11: Mapa estático para el README
# Insumo:    distritos_junin, zona_total, focos_junin, establecimientos_expuestos
# Resultado: salidas/mapa_junin.png
# ============================================================

fig, ax = plt.subplots(figsize=(8, 8.5))
# Coroplético: los distritos se colorean según la población expuesta
distritos_junin.plot(ax=ax, column="poblacion_expuesta", cmap="Purples", linewidth=0.3,
                     edgecolor="#777777", legend=True,
                     legend_kwds={"label": "Personas a ≤ 5 km de un foco (por distrito)", "shrink": 0.5})
zona_total.boundary.plot(ax=ax, color="#4d4d4d", linewidth=0.6, linestyle="--")
focos_junin.plot(ax=ax, color="#b10026", markersize=2)
establecimientos_expuestos.plot(ax=ax, color="#1565c0", marker="P", markersize=22, edgecolor="white",
                                linewidth=0.3)
ax.set_title(f"Junín: focos de calor VIIRS ({focos_junin['fecha_hora_lima'].min():%d %b} – "
             f"{focos_junin['fecha_hora_lima'].max():%d %b %Y}), zona de 5 km\n"
             f"y establecimientos de salud expuestos (cruz azul)", fontsize=10)
ax.axis("off")
plt.tight_layout()
plt.savefig(SALIDAS / "mapa_junin.png", dpi=150)
print(f"[11] mapa estático: {SALIDAS / 'mapa_junin.png'}")
