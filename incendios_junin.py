# ============================================================
# INCENDIOS EN JUNÍN · ¿A QUIÉN AFECTA Y PUEDE RESPONDER LA RED DE SALUD?
# Preguntas:
#   1. ¿Cuáles de los focos de calor de los últimos 7 días son probablemente incendios forestales?
#   2. ¿Cuántas personas y qué establecimientos de salud están a 5 km o menos de ellos?
#   3. ¿Los establecimientos expuestos tienen señal para telemedicina y un hospital a distancia
#      razonable para trasladar heridos?
# Camino D del profesor + valor agregado: criterios de evidencia y red de salud.
#
# Antes de este script (una vez al mes):  python red_salud_junin.py
# Cómo ejecutar (cada actualización):     python incendios_junin.py
# Resultados:                             carpeta salidas/ (mapas, tablas y página web)
# ============================================================

# --- Paquetes ---
# Instalar (una sola vez por computadora):
#   uv pip install -r requirements.txt
# Cargar (en cada ejecución): las líneas import de abajo.
import json                          # json: guarda los datos de la página web en formato JSON
import sys
import warnings
from datetime import datetime, timedelta, timezone
from pathlib import Path

import folium                        # folium: arma mapas interactivos (HTML) desde Python
import geopandas as gpd              # geopandas: tablas con geometría (puntos, polígonos); pandas + mapas
import matplotlib.pyplot as plt      # matplotlib: mapa estático (imagen PNG) para el README
import numpy as np                   # numpy: cálculos con arreglos de números (horas, distancias)
import pandas as pd                  # pandas: tablas de filas y columnas, como Excel pero con código
import planetary_computer            # planetary_computer: da permiso temporal para leer datos de Microsoft
import pystac_client                 # pystac_client: busca imágenes y mapas en catálogos STAC
import rasterio                      # rasterio: lee rasters (WorldCover, luces nocturnas)
import requests                      # requests: descarga archivos de internet
from rasterio.windows import Window  # rasterio: lee solo un recorte del raster, no el archivo entero
from rasterstats import zonal_stats  # rasterstats: resume un raster dentro de polígonos (estadística zonal)
from scipy.sparse import coo_matrix  # scipy: guarda qué focos están conectados con cuáles
from scipy.sparse.csgraph import connected_components  # scipy: encuentra los grupos de focos conectados
from scipy.spatial import cKDTree    # scipy: busca rápido los focos vecinos a 1 km o menos

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8")  # para que las tildes se vean bien en la consola de Windows


# ============================================================
# PROCESO 0: Configuración
# Insumo:    nada (son decisiones del análisis)
# Resultado: rutas, direcciones web, parámetros y funciones de apoyo en un solo lugar
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

# Parámetros del método
RADIO_M = 5_000          # zona de exposición: 5 km alrededor de cada foco
DIST_FOCOS_M = 1_000     # dos focos son del mismo incendio si están a <= 1 km ...
HORAS_FOCOS = 48         # ... y a <= 48 horas entre sí (regla del profesor)
LADO_PIXEL_M = 375       # tamaño del píxel VIIRS: cada foco representa un cuadrado de 375 m
HORAS_RECIENTE = 72      # "actividad reciente": último foco hace 72 h o menos
CRS_METROS = 32718       # UTM zona 18 Sur: sistema de coordenadas en METROS (Junín cae en esta zona)
NODATA_WORLDPOP = -99999

# Criterios de evidencia (Proceso 5)
MIN_FOCOS = 3            # un incendio forestal probable tiene al menos 3 focos ...
MIN_PASADAS = 2          # ... vistos en 2 o más pasadas, o por 2 o más satélites ...
PCT_VEGETACION = 50      # ... y la mayoría en vegetación natural (bosque, arbustal, pastizal)
UMBRAL_LUZ = 5           # radiancia nocturna (nW/cm²/sr) desde la que un píxel es ciudad o industria
VENTANA_WORLDCOVER = 37  # 37 píxeles de 10 m ≈ 375 m: el tamaño del píxel VIIRS

# Red de salud (Proceso 9): reglas de alerta
MIN_TRASLADO_LARGO = 120  # minutos al hospital de referencia desde los que el traslado es largo

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
# Cobertura del suelo: ESA WorldCover 2021 (10 m), desde Microsoft Planetary Computer (recurso del profesor)
STAC_PLANETARY = "https://planetarycomputer.microsoft.com/api/stac/v1"
ARCHIVO_COBERTURA_FOCOS = DATOS / "worldcover_focos.csv"   # memoria: cada foco se lee una sola vez
# Luces nocturnas VIIRS 2021 (1 km): el raster del notebook del profesor
URL_LUCES = ("https://huggingface.co/datasets/aquiro1994/ds-python-up/resolve/main/"
             "_data/VIIRS_NTL_Peru_YearlyComposite_2021.tif")
ARCHIVO_LUCES = DATOS / "VIIRS_NTL_Peru_2021.tif"
# Red de salud: capa base que produce red_salud_junin.py
ARCHIVO_RED_NODOS = SALIDAS / "red_salud_nodos.csv"
ARCHIVO_RED_ENLACES = SALIDAS / "red_salud_enlaces.csv"


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


def a_texto(tabla):
    """Convierte las fechas a texto: los formatos web (GeoJSON, JSON) no aceptan fechas."""
    tabla = tabla.copy()
    for col in tabla.columns:
        if pd.api.types.is_datetime64_any_dtype(tabla[col]):
            tabla[col] = tabla[col].dt.strftime("%Y-%m-%d %H:%M")
    return tabla


if not ARCHIVO_RED_NODOS.exists():
    sys.exit("Falta la capa base de la red de salud. Ejecute primero:  python red_salud_junin.py")

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
# anómala compatible con fuego. NO es un incendio confirmado: el Proceso 5 lo evalúa.
# Los focos se guardan en una carpeta con la fecha: cada actualización no borra las anteriores.

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

# Marca (no elimina) los focos de confianza baja
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
red_focos = coo_matrix((np.ones(len(pares)), (pares[:, 0], pares[:, 1])), shape=(n_focos, n_focos))
n_incendios, grupo = connected_components(red_focos, directed=False)

# Nombre legible: INC-001 es el incendio con más focos
orden = pd.Series(grupo).value_counts().index
codigo = {g: f"INC-{i + 1:03d}" for i, g in enumerate(orden)}
focos_junin["incendio"] = [codigo[g] for g in grupo]
focos_m["incendio"] = focos_junin["incendio"].values

# --- Paso 4.5: Zona de focos (huella) ---
# .buffer(..., cap_style="square"): dibuja el cuadrado de 375 m del píxel alrededor de cada foco
# .dissolve(): une los cuadrados del mismo incendio en un solo polígono
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
incendios["estado"] = np.select([horas_sin_focos < 24, horas_sin_focos < HORAS_RECIENTE],
                                ["Activo", "En observación"], default="Sin actividad reciente")
print(f"[4] {n_focos:,} focos -> {n_incendios} incendios | estados: {incendios['estado'].value_counts().to_dict()}")


# ============================================================
# PROCESO 5: Evaluar la evidencia de cada incendio
# Insumo:    focos_junin, incendios + ESA WorldCover 2021 (10 m) + luces nocturnas VIIRS 2021
# Resultado: columna "evidencia" en cada incendio
# ============================================================
# Un foco de calor no es un incendio confirmado. Se combinan cuatro pistas de datos abiertos:
#   1. Persistencia: el fuego se ve en 2 o más pasadas del satélite (no es un destello)
#   2. Coincidencia: lo ven 2 o más de los 3 satélites
#   3. Cobertura del suelo (WorldCover): ¿arde bosque o pastizal, un cultivo o una ciudad?
#   4. Luz nocturna (el raster del notebook del profesor): un foco en un píxel muy iluminado
#      sugiere ciudad o industria, no vegetación
# Límite: estas pistas reducen las falsas alarmas, pero no recuperan incendios ocultos por nubes.

CLASES_WORLDCOVER = {10: "Bosque", 20: "Arbustal", 30: "Pastizal", 40: "Cultivo", 50: "Urbano",
                     60: "Suelo desnudo", 70: "Nieve", 80: "Agua", 90: "Humedal", 95: "Manglar",
                     100: "Musgo", 0: "Sin dato"}
VEGETACION_NATURAL = [10, 20, 30, 90, 95, 100]


def leer_worldcover(puntos_por_leer):
    """Clase de cobertura más frecuente en un cuadrado de ~375 m alrededor de cada punto."""
    # pystac_client.Client.open(): abre el catálogo de Planetary Computer
    # planetary_computer.sign_inplace: agrega el permiso temporal para leer los archivos
    catalogo = pystac_client.Client.open(STAC_PLANETARY, modifier=planetary_computer.sign_inplace)
    caja = [puntos_por_leer["longitude"].min(), puntos_por_leer["latitude"].min(),
            puntos_por_leer["longitude"].max(), puntos_por_leer["latitude"].max()]
    # .search(): busca las teselas de WorldCover que tocan la caja; se usan las de 2021 (versión 2)
    teselas = [t for t in catalogo.search(collections=["esa-worldcover"], bbox=caja).items() if "2021" in t.id]
    medio = VENTANA_WORLDCOVER // 2
    leidos = []
    for tesela in teselas:
        with rasterio.open(tesela.assets["map"].href) as fuente:
            b = fuente.bounds
            dentro = puntos_por_leer[puntos_por_leer["longitude"].between(b.left, b.right)
                                     & puntos_por_leer["latitude"].between(b.bottom, b.top)]
            for _, p in dentro.iterrows():
                fila, columna = fuente.index(p["longitude"], p["latitude"])
                col0, fila0 = max(columna - medio, 0), max(fila - medio, 0)
                # Window(): solo se descarga un recorte de 37 × 37 píxeles, no la tesela entera
                ventana = Window(col0, fila0, min(VENTANA_WORLDCOVER, fuente.width - col0),
                                 min(VENTANA_WORLDCOVER, fuente.height - fila0))
                bloque = fuente.read(1, window=ventana)
                validos = bloque[bloque > 0]   # 0 = sin dato
                # np.bincount().argmax(): la clase que más se repite (como CONTAR.SI y elegir la mayor)
                leidos.append((p["clave"], int(np.bincount(validos).argmax()) if validos.size else 0))
    return pd.DataFrame(leidos, columns=["clave", "clase_worldcover"])


# --- Paso 5.1: Cobertura del suelo de cada foco (con memoria para no releer) ---
focos_junin["clave"] = (focos_junin["latitude"].round(5).astype(str) + "_"
                        + focos_junin["longitude"].round(5).astype(str))
memoria = (pd.read_csv(ARCHIVO_COBERTURA_FOCOS, dtype={"clave": str}) if ARCHIVO_COBERTURA_FOCOS.exists()
           else pd.DataFrame(columns=["clave", "clase_worldcover"]))
faltan = (focos_junin.loc[~focos_junin["clave"].isin(memoria["clave"]), ["clave", "latitude", "longitude"]]
          .drop_duplicates("clave"))
print(f"[5] WorldCover: {len(faltan):,} focos por leer ({len(memoria):,} ya estaban en memoria)")
if len(faltan):
    memoria = pd.concat([memoria, leer_worldcover(faltan)], ignore_index=True).drop_duplicates("clave")
    memoria.to_csv(ARCHIVO_COBERTURA_FOCOS, index=False)
focos_junin["clase_worldcover"] = (focos_junin["clave"].map(memoria.set_index("clave")["clase_worldcover"])
                                   .fillna(0).astype(int))

# --- Paso 5.2: Luz nocturna en cada foco ---
descargar(URL_LUCES, ARCHIVO_LUCES)
with rasterio.open(ARCHIVO_LUCES) as fuente:
    # .sample(): lee el valor del raster en cada punto (longitud, latitud)
    luz = np.array([v[0] for v in fuente.sample(zip(focos_junin["longitude"], focos_junin["latitude"]))], dtype=float)
# Fuera del Perú el raster tiene NaN (la lección del notebook): se reemplaza por 0
focos_junin["luz_nocturna"] = np.nan_to_num(luz, nan=0.0)

# --- Paso 5.3: Entorno de cada foco ---
# np.select(): como un SI anidado; gana la primera condición que se cumple
focos_junin["entorno"] = np.select(
    [focos_junin["clase_worldcover"].eq(50) | focos_junin["luz_nocturna"].ge(UMBRAL_LUZ),
     focos_junin["clase_worldcover"].isin(VEGETACION_NATURAL),
     focos_junin["clase_worldcover"].eq(40)],
    ["Urbano o industrial", "Vegetación natural", "Cultivo"], default="Otro")

# --- Paso 5.4: Pistas por incendio ---
# Una "pasada" es la fecha local + día o noche: dos pasadas distintas = fuego visto en momentos distintos
focos_junin["pasada"] = focos_junin["fecha_hora_lima"].dt.strftime("%Y-%m-%d") + focos_junin["daynight"]
pistas = focos_junin.groupby("incendio").agg(
    pasadas=("pasada", "nunique"),
    satelites=("satelite", "nunique"),
    pct_vegetacion=("entorno", lambda s: round(100 * s.eq("Vegetación natural").mean())),
    pct_cultivo=("entorno", lambda s: round(100 * s.eq("Cultivo").mean())),
    pct_urbano=("entorno", lambda s: round(100 * s.eq("Urbano o industrial").mean())),
    pct_confianza_baja=("confianza_baja", lambda s: round(100 * s.mean())),
    cobertura_principal=("clase_worldcover", lambda s: CLASES_WORLDCOVER.get(int(s.mode().iat[0]), "Otro")))
incendios = incendios.join(pistas)

# --- Paso 5.5: Nivel de evidencia (regla explícita, en orden) ---
incendios["evidencia"] = np.select(
    [incendios["pct_urbano"] >= 50,
     (incendios["focos"] >= MIN_FOCOS)
     & ((incendios["pasadas"] >= MIN_PASADAS) | (incendios["satelites"] >= MIN_PASADAS))
     & (incendios["pct_vegetacion"] >= PCT_VEGETACION)
     & (incendios["pct_confianza_baja"] < 50),
     incendios["pct_cultivo"] >= 50],
    ["Probable fuente no forestal", "Incendio forestal probable", "Probable quema agropecuaria"],
    default="Detección por verificar")
focos_junin["evidencia"] = focos_junin["incendio"].map(incendios["evidencia"])
focos_m["evidencia"] = focos_junin["evidencia"].values
print(f"    evidencia: {incendios['evidencia'].value_counts().to_dict()}")
print(f"    focos por entorno: {focos_junin['entorno'].value_counts().to_dict()}")


# ============================================================
# PROCESO 6: Construir la zona de 5 km
# Insumo:    focos_m (focos en metros) + evidencia
# Resultado: zonas_5km (un polígono por incendio) y zonas unidas sin doble conteo
# ============================================================

# .buffer(5000): dibuja un círculo de 5 km alrededor de cada foco
# .dissolve(by="incendio"): une los círculos de un mismo incendio
#   Antes: un punto por foco | Después: un polígono de 5 km por incendio
circulos = focos_m[["incendio", "geometry"]].copy()
circulos["geometry"] = focos_m.buffer(RADIO_M)
zonas_5km_m = circulos.dissolve(by="incendio")


def zona_unida(codigos):
    """Une las zonas de 5 km de los incendios indicados en una sola (sin doble conteo), en grados."""
    # .union_all(): si dos zonas se superponen, la parte común queda UNA sola vez
    geometrias = [zonas_5km_m.loc[codigos].union_all()] if len(codigos) else []
    return gpd.GeoDataFrame(geometry=geometrias, crs=CRS_METROS).to_crs(epsg=4326)


incendios_probables = incendios.index[incendios["evidencia"] == "Incendio forestal probable"]
# "Para decidir hoy": incendios forestales probables con actividad en las últimas 72 h
incendios_recientes = incendios.index[(incendios["evidencia"] == "Incendio forestal probable")
                                      & (incendios["estado"] != "Sin actividad reciente")]
zona_total = zona_unida(incendios.index)          # todos los focos (para comparar)
zona_probables = zona_unida(incendios_probables)
zona_recientes = zona_unida(incendios_recientes)

# El raster de población está en grados (EPSG:4326): las zonas deben estar en el mismo sistema
zonas_5km = zonas_5km_m.to_crs(epsg=4326)
huella = huella_m.to_crs(epsg=4326)
print(f"[6] incendios forestales probables: {len(incendios_probables)} "
      f"| con actividad en <= {HORAS_RECIENTE} h: {len(incendios_recientes)}")


# ============================================================
# PROCESO 7: Contar la población con WorldPop (dato raster)
# Insumo:    raster WorldPop 2026 (100 m) + zonas de 5 km + distritos_junin
# Resultado: población a 5 km por incendio, en totales sin doble conteo y por distrito
# ============================================================

descargar(URL_WORLDPOP, ARCHIVO_WORLDPOP)


def poblacion(poligonos):
    """Suma las personas de WorldPop dentro de cada polígono (estadística zonal, como en el notebook)."""
    if len(poligonos) == 0:
        return []
    # zonal_stats(): recorta el raster con cada polígono y suma sus celdas
    #   Antes: polígonos sin población | Después: una suma de personas por polígono
    #   nodata explícito: así las celdas "sin dato" (-99999) no entran a la suma
    resultado = zonal_stats(poligonos, str(ARCHIVO_WORLDPOP), stats=["sum"], nodata=NODATA_WORLDPOP)
    return [round(r["sum"] or 0) for r in resultado]


def poblacion_unica(zona):
    """Población de una zona unida (0 si la zona está vacía)."""
    valores = poblacion(zona)
    return valores[0] if valores else 0


# --- Paso 7.1: Por incendio (las zonas pueden superponerse: NO sumar esta columna) ---
zonas_5km["poblacion_5km"] = poblacion(zonas_5km)
incendios["poblacion_5km"] = zonas_5km["poblacion_5km"]

# --- Paso 7.2: Tres totales sin doble conteo ---
poblacion_total_5km = poblacion_unica(zona_total)              # cualquier foco
poblacion_probables_5km = poblacion_unica(zona_probables)      # incendios forestales probables
poblacion_recientes_5km = poblacion_unica(zona_recientes)      # ... con actividad en 72 h
print(f"[7] personas a <= 5 km: cualquier foco {poblacion_total_5km:,} | incendios probables "
      f"{poblacion_probables_5km:,} | probables con actividad en 72 h {poblacion_recientes_5km:,}")

# --- Paso 7.3: Por distrito de Junín (zona de incendios forestales probables) ---
# gpd.overlay(..., how="intersection"): recorta la zona con cada distrito
#   Antes: una zona unida | Después: un pedazo de zona por distrito
distritos_junin["poblacion_total"] = poblacion(distritos_junin)
if len(zona_probables):
    pedazos = gpd.overlay(distritos_junin[["ubigeo", "geometry"]], zona_probables, how="intersection")
    pedazos["poblacion_expuesta"] = poblacion(pedazos)
    distritos_junin = distritos_junin.merge(pedazos[["ubigeo", "poblacion_expuesta"]], on="ubigeo", how="left")
else:
    distritos_junin["poblacion_expuesta"] = 0
distritos_junin["poblacion_expuesta"] = distritos_junin["poblacion_expuesta"].fillna(0).astype(int)
distritos_junin["pct_expuesta"] = (100 * distritos_junin["poblacion_expuesta"]
                                   / distritos_junin["poblacion_total"].replace(0, np.nan)).round(1)
print(f"    población de Junín: {distritos_junin['poblacion_total'].sum():,} "
      f"| expuesta a incendios probables, dentro de Junín: {distritos_junin['poblacion_expuesta'].sum():,}")


# ============================================================
# PROCESO 8: Establecimientos de salud a 5 km o menos (RENIPRESS)
# Insumo:    CSV RENIPRESS (SUSALUD) + focos_m
# Resultado: establecimientos_expuestos (a <= 5 km de un incendio forestal probable)
# ============================================================

descargar(URL_RENIPRESS, ARCHIVO_RENIPRESS)

# El archivo usa ";" como separador y UTF-8 con BOM (por eso "utf-8-sig")
renipress = pd.read_csv(ARCHIVO_RENIPRESS, sep=";", encoding="utf-8-sig", dtype=str)
# pd.to_numeric(): convierte las coordenadas de texto a número; si no se puede, deja vacío (NaN)
renipress["lat"] = pd.to_numeric(renipress["NORTE"], errors="coerce")
renipress["lon"] = pd.to_numeric(renipress["ESTE"], errors="coerce")

# --- Paso 8.1: Filtros en cascada (se informa cuántos quedan en cada paso) ---
filtros = [
    ("registros en RENIPRESS", pd.Series(True, index=renipress.index)),
    ("estado ACTIVO en el registro", renipress["ESTADO"].eq("ACTIVO")),
    ("establecimiento de salud (no servicio de apoyo)", renipress["TIPO_ESTABLECIMIENTO"].str.startswith("ESTABLECIMIENTO")),
    ("con categoría asignada (no '0')", renipress["CATEGORIA"].ne("0")),
    ("con coordenadas dentro del Perú", renipress["lat"].between(-18.4, -0.01) & renipress["lon"].between(-81.4, -68.6)),
]
conserva = pd.Series(True, index=renipress.index)
print("[8] filtros de RENIPRESS (todo el Perú):")
for descripcion, condicion in filtros:
    conserva &= condicion
    print(f"    {descripcion:50s} {conserva.sum():>7,}")
establecimientos = renipress[conserva].copy()
establecimientos = gpd.GeoDataFrame(establecimientos,
                                    geometry=gpd.points_from_xy(establecimientos["lon"], establecimientos["lat"]),
                                    crs="EPSG:4326")
establecimientos["red"] = np.where(establecimientos["INSTITUCION"].eq("PRIVADO"), "Privada", "Pública")
establecimientos_m = establecimientos.to_crs(epsg=CRS_METROS)


def mas_cercano_a(focos_elegidos):
    """Establecimientos con un foco a <= 5 km y la distancia a ese foco (uno por establecimiento)."""
    if len(focos_elegidos) == 0:
        return establecimientos_m.iloc[0:0].assign(incendio=None, dist_foco_m=np.nan)
    # gpd.sjoin_nearest(): pega a cada establecimiento el foco más cercano a <= 5 km y la distancia
    #   Antes: establecimientos sin distancia | Después: solo los que tienen un foco a <= 5 km
    cercanos = gpd.sjoin_nearest(establecimientos_m, focos_elegidos[["incendio", "geometry"]],
                                 how="inner", max_distance=RADIO_M, distance_col="dist_foco_m")
    # Si hay empate entre dos focos, el establecimiento sale dos veces: se deja una sola fila
    return cercanos.drop_duplicates(subset="COD_IPRESS").drop(columns="index_right")


# --- Paso 8.2: Cerca de cualquier foco (para comparar) y cerca de incendios probables ---
cerca_cualquier = mas_cercano_a(focos_m)
incendios["establecimientos_5km"] = (cerca_cualquier.groupby("incendio").size()
                                     .reindex(incendios.index).fillna(0).astype(int))
expuestos_m = mas_cercano_a(focos_m[focos_m["incendio"].isin(incendios_probables)])
establecimientos_expuestos = expuestos_m.to_crs(epsg=4326)
establecimientos_expuestos["dist_foco_km"] = (establecimientos_expuestos["dist_foco_m"] / 1000).round(2)
# Cerca de incendios probables con actividad en 72 h
salud_recientes = (gpd.sjoin(establecimientos_expuestos, zona_recientes, predicate="within", how="inner")
                   if len(zona_recientes) else establecimientos_expuestos.iloc[0:0])
# Establecimientos de Junín (para mostrarlos todos en el mapa)
establecimientos_junin = gpd.sjoin(establecimientos, distritos_junin[["ubigeo", "geometry"]],
                                   predicate="within", how="inner").drop(columns=["index_right", "ubigeo"])
print(f"    establecimientos en Junín: {len(establecimientos_junin):,} | a <= 5 km de cualquier foco: "
      f"{len(cerca_cualquier):,} | de un incendio probable: {len(establecimientos_expuestos):,} "
      f"| de uno con actividad en 72 h: {len(salud_recientes):,}")


# ============================================================
# PROCESO 9: La red de salud ante el fuego
# Insumo:    red_salud_nodos.csv (capa base: grafo, tiempo de viaje y 4G) + focos de incendios probables
# Resultado: alerta por establecimiento de la red y resumen por hospital de referencia
# ============================================================
# Pregunta: si un establecimiento expuesto recibe heridos, ¿puede hacer una teleconsulta (4G)
# y trasladarlos a su hospital de referencia en un tiempo razonable?

# pd.read_csv(): lee la capa base que produjo red_salud_junin.py (no se recalcula aquí)
red = pd.read_csv(ARCHIVO_RED_NODOS, encoding="utf-8-sig",
                  dtype={"COD_IPRESS": str, "padre": str, "hospital_ref": str, "UBIGEO": str})
red = gpd.GeoDataFrame(red, geometry=gpd.points_from_xy(red["lon"], red["lat"]), crs="EPSG:4326")
red_m = red.to_crs(epsg=CRS_METROS)

# --- Paso 9.1: Distancia de cada nodo al foco de incendio probable más cercano ---
focos_probables_m = focos_m[focos_m["incendio"].isin(incendios_probables)]
if len(focos_probables_m):
    cercania = (gpd.sjoin_nearest(red_m[["COD_IPRESS", "geometry"]], focos_probables_m[["incendio", "geometry"]],
                                  how="left", distance_col="dist_m")
                .drop_duplicates("COD_IPRESS").set_index("COD_IPRESS"))
    red["dist_incendio_km"] = red["COD_IPRESS"].map(cercania["dist_m"] / 1000).round(2)
    red["incendio_cercano"] = red["COD_IPRESS"].map(cercania["incendio"])
else:
    red["dist_incendio_km"], red["incendio_cercano"] = np.nan, None
red["expuesto"] = red["dist_incendio_km"] <= RADIO_M / 1000
red["expuesto_reciente"] = (red["COD_IPRESS"].isin(
    gpd.sjoin(red[["COD_IPRESS", "geometry"]], zona_recientes, predicate="within")["COD_IPRESS"])
    if len(zona_recientes) else False)

# --- Paso 9.2: Alerta por nodo (regla explícita) ---
#   Crítico:  primer nivel expuesto, SIN 4G garantizado Y a más de 2 h de su hospital
#   Atención: primer nivel expuesto, sin 4G O a más de 2 h
#   Expuesto: expuesto, con 4G y hospital a menos de 2 h (o es un hospital)
primer_nivel = red["CATEGORIA"].str.startswith("I-")
sin_4g = red["conectividad"].ne("4G garantizado")
traslado_largo = red["min_a_hospital"] > MIN_TRASLADO_LARGO
red["alerta"] = np.select(
    [red["expuesto"] & primer_nivel & sin_4g & traslado_largo,
     red["expuesto"] & primer_nivel & (sin_4g | traslado_largo),
     red["expuesto"]],
    ["Crítico", "Atención", "Expuesto"], default="Sin exposición")

# --- Paso 9.3: Ramas del árbol que conectan a los nodos expuestos con su hospital ---
padre_de = red.set_index("COD_IPRESS")["padre"].to_dict()
en_ruta = set()
for codigo_nodo in red.loc[red["expuesto"], "COD_IPRESS"]:
    actual = codigo_nodo
    while isinstance(actual, str) and actual not in en_ruta:   # sube hasta la raíz del árbol
        en_ruta.add(actual)
        actual = padre_de.get(actual)
red["en_ruta_expuesta"] = red["COD_IPRESS"].isin(en_ruta)

# --- Paso 9.4: Guardar alertas y resumen por hospital de referencia ---
orden_alerta = {"Crítico": 0, "Atención": 1, "Expuesto": 2, "Sin exposición": 3}
alertas = (red[red["expuesto"]].assign(orden=red["alerta"].map(orden_alerta))
           .sort_values(["orden", "min_a_hospital"], ascending=[True, False]))
columnas_alerta = ["alerta", "COD_IPRESS", "NOMBRE", "sistema", "CATEGORIA", "rol", "RED", "MICRORED", "DISTRITO",
                   "nombre_padre", "nombre_hospital_ref", "min_a_hospital", "conectividad", "cobertura_4g_pct",
                   "incendio_cercano", "dist_incendio_km", "expuesto_reciente"]
alertas[columnas_alerta].to_csv(SALIDAS / "red_salud_alertas.csv", index=False, encoding="utf-8-sig")
por_hospital = (red[red["expuesto"] & red["nombre_hospital_ref"].notna()]
                .groupby(["sistema", "nombre_hospital_ref"])
                .agg(expuestos=("COD_IPRESS", "size"),
                     criticos=("alerta", lambda s: s.eq("Crítico").sum()),
                     atencion=("alerta", lambda s: s.eq("Atención").sum()),
                     recientes=("expuesto_reciente", "sum"),
                     minutos_mediana=("min_a_hospital", "median"))
                .sort_values(["criticos", "expuestos"], ascending=False))
por_hospital.to_csv(SALIDAS / "red_salud_por_hospital.csv", encoding="utf-8-sig")
print(f"[9] red de salud: {len(red)} nodos | expuestos: {int(red['expuesto'].sum())} "
      f"| alertas: {red.loc[red['expuesto'], 'alerta'].value_counts().to_dict()}")
print(por_hospital.round(0).to_string())


# ============================================================
# PROCESO 10: Guardar las tablas de resultados
# Insumo:    incendios, distritos_junin, establecimientos_expuestos
# Resultado: archivos CSV en salidas/
# ============================================================

incendios = incendios.sort_values("poblacion_5km", ascending=False)
tabla_incendios = incendios.reset_index()[[
    "incendio", "evidencia", "estado", "provincia", "distrito", "ubigeo", "distritos_con_focos", "focos",
    "pasadas", "satelites", "cobertura_principal", "pct_vegetacion", "pct_cultivo", "pct_urbano",
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
print(f"[10] tablas guardadas en {SALIDAS}")
print(tabla_incendios.head(8)[["incendio", "evidencia", "estado", "distrito", "focos",
                               "poblacion_5km", "establecimientos_5km"]].to_string(index=False))


# ============================================================
# PROCESO 11: Ficha territorial para el Gobierno Regional
# Insumo:    focos_junin, incendios, distritos_junin, establecimientos_expuestos, red
# Resultado: ficha_distritos_junin.csv (prioridad por distrito)
# ============================================================
# Enfoque territorial: el Gobierno Regional coordina con las municipalidades (distritos y
# provincias) y con su red de salud. Por eso los resultados se resumen por distrito.

# --- Paso 11.1: Focos por distrito (7 días y últimas 72 h) ---
# .groupby("ubigeo").agg(): cuenta los focos de cada distrito (como una tabla dinámica)
focos_junin["reciente"] = focos_junin["horas_desde_foco"] <= HORAS_RECIENTE
focos_por_distrito = focos_junin.groupby("ubigeo").agg(focos_7d=("ubigeo", "size"),
                                                       focos_72h=("reciente", "sum"))

# --- Paso 11.2: Incendios forestales probables con actividad reciente, por distrito ---
incendios_por_distrito = (focos_junin[focos_junin["incendio"].isin(incendios_recientes)]
                          .groupby("ubigeo")["incendio"].nunique().rename("incendios_prioritarios"))

# --- Paso 11.3: Establecimientos expuestos y puntos críticos de la red, por distrito ---
salud_por_distrito = establecimientos_expuestos.groupby("UBIGEO").size().rename("establecimientos_expuestos")
criticos_por_distrito = red[red["alerta"].eq("Crítico")].groupby("UBIGEO").size().rename("puntos_criticos")

# --- Paso 11.4: Unir todo en una ficha por distrito ---
# .join(): pega columnas usando el UBIGEO como llave (como BUSCARV); sin dato -> 0
#   Antes: 5 tablas separadas | Después: una fila por distrito con todos los indicadores
ficha = (distritos_junin.drop(columns="geometry").set_index("ubigeo")
         .join([focos_por_distrito, incendios_por_distrito, salud_por_distrito, criticos_por_distrito]))
conteos = ["focos_7d", "focos_72h", "incendios_prioritarios", "establecimientos_expuestos", "puntos_criticos"]
ficha[conteos] = ficha[conteos].fillna(0).astype(int)

# --- Paso 11.5: Semáforo de prioridad (regla simple y explicable) ---
# np.select(): como un SI anidado de Excel; gana la primera condición que se cumple
ficha["prioridad"] = np.select(
    [ficha["incendios_prioritarios"] > 0,   # Alta: incendio forestal probable con actividad en 72 h
     ficha["focos_72h"] > 0,                # Media: focos recientes, sin evidencia suficiente
     ficha["focos_7d"] > 0],                # Baja: solo focos de hace más de 72 h
    ["1 · Alta", "2 · Media", "3 · Baja"], default="Sin focos")
ficha = ficha.sort_values(["prioridad", "poblacion_expuesta"], ascending=[True, False])
ficha.reset_index().to_csv(SALIDAS / "ficha_distritos_junin.csv", index=False, encoding="utf-8-sig")

# .map(): copia la prioridad y los conteos a la tabla con geometría (para colorear el mapa)
for col in ["prioridad"] + conteos:
    distritos_junin[col] = distritos_junin["ubigeo"].map(ficha[col])
print(f"[11] prioridad por distrito: {ficha['prioridad'].value_counts().sort_index().to_dict()}")


# ============================================================
# PROCESO 12: Mapa interactivo (Folium)
# Insumo:    todas las capas anteriores
# Resultado: salidas/mapa_incendios_junin.html (capas que se activan y desactivan)
# ============================================================
# Paleta: cada tema tiene su propio color para que no se confundan
#   morado = territorio · rojo/naranja/amarillo = fuego · verde azulado = red de salud
#   (DIRESA: línea continua y cuadrado; EsSalud: línea punteada y círculo) · gris punteado = 5 km
COLOR_PRIORIDAD = {"1 · Alta": "#54278f", "2 · Media": "#9e9ac8", "3 · Baja": "#dadaeb", "Sin focos": "#ffffff"}
COLOR_SALUD = "#008f8f"
LINEA_SISTEMA = {"DIRESA Junín": None, "EsSalud": "6,4"}      # EsSalud: línea punteada
COLOR_ALERTA = {"Crítico": "#b3261e", "Atención": "#d98c00", "Expuesto": "#4d4d4d"}

# folium.Map(): crea el mapa base centrado en Junín
mapa = folium.Map(location=[(lat_min + lat_max) / 2, (lon_min + lon_max) / 2], zoom_start=8,
                  tiles="CartoDB positron", control_scale=True)
folium.TileLayer("OpenStreetMap", name="Calles (OpenStreetMap)").add_to(mapa)
folium.TileLayer(tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
                 attr="Esri World Imagery", name="Imagen satelital (Esri)").add_to(mapa)

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
                "establecimientos_expuestos", "puntos_criticos", "poblacion_total"],
        aliases=["Distrito", "Provincia", "UBIGEO", "Prioridad", "Focos en las últimas 72 h",
                 "Personas a ≤ 5 km de incendios probables", "Establecimientos de salud a ≤ 5 km",
                 "Puntos críticos de la red", "Población 2026 (WorldPop)"], localize=True),
).add_to(mapa)

# --- Capa 2: Límites provinciales (línea gruesa, sin relleno) ---
# .dissolve(by="nombprov"): une los distritos de cada provincia en un solo polígono
provincias = distritos_junin[["nombprov", "geometry"]].dissolve(by="nombprov").reset_index()
provincias["geometry"] = provincias.simplify(0.001)
folium.GeoJson(provincias, name="Provincias (límite)",
               style_function=lambda f: {"color": "#3f007d", "weight": 3.5, "fill": False}).add_to(mapa)

# --- Capa 3: Área de exposición de 5 km de los incendios probables (NO es el avance del fuego) ---
zonas_web = zonas_5km.loc[incendios_probables].reset_index()[["incendio", "poblacion_5km", "geometry"]]
zonas_web["geometry"] = zonas_web.simplify(0.0005)
folium.GeoJson(
    zonas_web, name="Área de exposición: 5 km de incendios probables",
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
    # Color según antigüedad (rojo < 24 h, naranja 24–72 h, amarillo > 72 h) solo si el foco es de un
    # incendio forestal probable; gris si no lo es (mismo código de colores que la leyenda)
    if foco["evidencia"] != "Incendio forestal probable":
        color = "#a3a0ad"
    else:
        color = "#d7191c" if foco["horas_desde_foco"] < 24 else "#fd8d3c" if foco["horas_desde_foco"] < 72 else "#fecc5c"
    folium.CircleMarker(
        location=[foco["latitude"], foco["longitude"]], radius=3, color=color, weight=0.5,
        fill=True, fill_color=color, fill_opacity=0.9,
        tooltip=(f"Foco {foco['fecha_hora_lima']} (Lima) · {foco['satelite']} · FRP {foco['frp']:.1f} MW"
                 f" · {CLASES_WORLDCOVER.get(foco['clase_worldcover'], 'Otro')} · {foco['incendio']}: {foco['evidencia']}"),
    ).add_to(capa_focos)
capa_focos.add_to(mapa)

# --- Capa 6: Incendios según su evidencia (clic = resumen) ---
# folium.DivIcon(): un ícono hecho con HTML
ICONO_EVIDENCIA = {"Incendio forestal probable": "🔥", "Probable quema agropecuaria": "🌾",
                   "Probable fuente no forestal": "🏭"}
capa_incendios = folium.FeatureGroup(name="Incendios: forestal probable 🔥, quema agropecuaria 🌾, no forestal 🏭")
for _, inc in a_texto(tabla_incendios[tabla_incendios["evidencia"].isin(ICONO_EVIDENCIA)]).iterrows():
    tam = {"Activo": 26, "En observación": 20}.get(inc["estado"], 15)
    resumen = (f"<b>{inc['incendio']}</b> · {inc['evidencia']}<br>{inc['estado']} · "
               f"{inc['distrito']} ({inc['provincia']})<br>"
               f"{inc['focos']} focos en {inc['pasadas']} pasadas, {inc['satelites']} satélites · "
               f"cobertura principal: {inc['cobertura_principal']}<br>"
               f"Huella VIIRS: {inc['huella_ha']:,.0f} ha · del {inc['inicio']} al {inc['ultimo_foco']}<br>"
               f"<b>Personas a ≤ 5 km: {inc['poblacion_5km']:,}</b> · "
               f"establecimientos de salud a ≤ 5 km: {inc['establecimientos_5km']}")
    folium.Marker(location=[inc["lat"], inc["lon"]],
                  icon=folium.DivIcon(html=f'<div style="font-size:{tam}px;line-height:1">'
                                           f'{ICONO_EVIDENCIA[inc["evidencia"]]}</div>',
                                      icon_size=(tam, tam), icon_anchor=(tam // 2, tam // 2)),
                  tooltip=f"{inc['incendio']} · {inc['evidencia']} · {inc['distrito']}",
                  popup=folium.Popup(resumen, max_width=340)).add_to(capa_incendios)
capa_incendios.add_to(mapa)

# --- Capa 7: Ramas de la red que conectan a los establecimientos expuestos con su hospital ---
coord_nodo = red.set_index("COD_IPRESS")[["lat", "lon"]]


def dibujar_enlaces(nodos_hijos, grupo):
    """Dibuja una línea de cada establecimiento a su referencia (color según el sistema)."""
    for _, n in nodos_hijos[nodos_hijos["padre"].notna()].iterrows():
        destino = coord_nodo.loc[n["padre"]]
        folium.PolyLine([[n["lat"], n["lon"]], [destino["lat"], destino["lon"]]],
                        color=COLOR_SALUD, weight=1.6, opacity=0.8, dash_array=LINEA_SISTEMA[n["sistema"]],
                        tooltip=f"{n['NOMBRE']} → {n['nombre_padre']} · {n['min_a_padre']:.0f} min").add_to(grupo)


capa_ramas = folium.FeatureGroup(name="Red de referencia: ramas de los establecimientos expuestos")
dibujar_enlaces(red[red["en_ruta_expuesta"]], capa_ramas)
capa_ramas.add_to(mapa)
capa_red = folium.FeatureGroup(name="Red de referencia completa (DIRESA y EsSalud)", show=False)
dibujar_enlaces(red, capa_red)
capa_red.add_to(mapa)


def icono_nodo(sistema, alerta):
    """Cuadrado (DIRESA) o círculo (EsSalud) con cruz; el borde grueso marca la alerta."""
    color = COLOR_SALUD
    borde = COLOR_ALERTA.get(alerta, color)
    radio = "50%" if sistema == "EsSalud" else "3px"
    grosor = 3 if alerta in ("Crítico", "Atención") else 1.5
    return folium.DivIcon(
        html=(f'<div style="background:{color};color:white;border:{grosor}px solid {borde};border-radius:{radio};'
              f'width:14px;height:14px;font:bold 13px/14px Arial;text-align:center;'
              f'box-shadow:0 0 3px white">+</div>'),
        icon_size=(20, 20), icon_anchor=(10, 10))


# --- Capa 8: Establecimientos de la red a ≤ 5 km de un incendio probable (con su alerta) ---
capa_alertas = folium.FeatureGroup(name="Salud: red pública expuesta (borde rojo = crítico, ámbar = atención)")
for _, n in red[red["expuesto"]].iterrows():
    folium.Marker(
        location=[n["lat"], n["lon"]], icon=icono_nodo(n["sistema"], n["alerta"]),
        tooltip=(f"<b>{n['alerta']}</b> · {n['NOMBRE']} ({n['CATEGORIA']}, {n['sistema']})<br>"
                 f"Refiere a: {n['nombre_padre']} · hospital: {n['nombre_hospital_ref']}<br>"
                 f"Traslado estimado al hospital: {n['min_a_hospital']:.0f} min · {n['conectividad']}"
                 f" ({n['cobertura_4g_pct']:.0f} % del centro poblado con 4G)<br>"
                 f"A {n['dist_incendio_km']} km de {n['incendio_cercano']}"),
    ).add_to(capa_alertas)
capa_alertas.add_to(mapa)

# --- Capa 9: Todos los establecimientos de salud de Junín (apagada al inicio) ---
capa_salud = folium.FeatureGroup(name="Salud: todos los establecimientos de Junín", show=False)
for _, es in establecimientos_junin.iterrows():
    folium.CircleMarker(location=[es["lat"], es["lon"]], radius=2.5, color="#008f8f", weight=0.5,
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
    f"<td align='right'>{d.establecimientos_expuestos}</td><td align='right'>{d.puntos_criticos}</td></tr>"
    for d in alta.itertuples())
hospitales = ", ".join(f"{h.title()} ({int(f.expuestos)}, {int(f.criticos)} críticos)"
                       for (_, h), f in por_hospital.head(4).iterrows())
n_criticos = int(red["alerta"].eq("Crítico").sum())
panel = f"""
<details open style="position:fixed; top:80px; left:12px; z-index:9999; background:white; padding:8px 12px;
         border:2px solid #3f007d; border-radius:6px; font:12px/1.4 Arial; max-width:410px; max-height:55vh; overflow:auto;">
<summary style="cursor:pointer; font-weight:bold; color:#3f007d; font-size:13px">
Para decidir hoy · Gobierno Regional de Junín</summary>
<div style="margin-top:6px">
<b>{len(incendios_recientes)}</b> incendios forestales probables con actividad en las últimas 72 h ·
<b>{poblacion_recientes_5km:,}</b> personas a ≤ 5 km ·
<b>{len(salud_recientes)}</b> establecimientos de salud a ≤ 5 km ·
<b>{n_criticos}</b> puntos críticos de la red<br><br>
<b>Distritos de prioridad alta</b> (incendio forestal probable con actividad en 72 h):
<table style="border-collapse:collapse; width:100%; margin-top:4px">
<tr style="border-bottom:1px solid #999"><th align="left">Distrito</th><th align="left">Provincia</th>
<th>Focos 72 h</th><th>Personas ≤ 5 km</th><th>Salud ≤ 5 km</th><th>Críticos</th></tr>
{filas}
</table><br>
<b>Hospitales que recibirían heridos</b> (establecimientos expuestos que refieren a cada uno): {hospitales or "ninguno"}<br>
<i>Crítico = puesto o centro de salud expuesto, sin 4G garantizado y a más de 2 h de su hospital.
Detalle: salidas/red_salud_alertas.csv</i>
</div></details>"""
mapa.get_root().html.add_child(folium.Element(panel))

# --- Leyenda: cómo leer el mapa ---
leyenda = f"""
<details open style="position:fixed; bottom:24px; left:12px; z-index:9999; background:white; padding:8px 12px;
         border:1px solid #999; border-radius:6px; font:12px/1.55 Arial; max-width:360px;">
<summary style="cursor:pointer; font-weight:bold">Cómo leer el mapa</summary>
<b>Territorio</b> (morado): <span style="border-bottom:3px solid #3f007d">&nbsp;provincia&nbsp;</span>
<span style="border-bottom:1px solid #3f007d">&nbsp;distrito&nbsp;</span>. Prioridad del distrito:<br>
<span style="background:#54278f;opacity:.6">&nbsp;&nbsp;&nbsp;</span> <b>alta</b>: incendio forestal probable con actividad en 72 h<br>
<span style="background:#9e9ac8;opacity:.6">&nbsp;&nbsp;&nbsp;</span> <b>media</b>: focos en 72 h sin evidencia de incendio forestal<br>
<span style="background:#dadaeb">&nbsp;&nbsp;&nbsp;</span> <b>baja</b>: solo focos de hace más de 72 h · sin color: sin focos en 7 días<br>
<b>Fuego</b>: 🔥 incendio forestal probable · 🌾 probable quema agropecuaria · 🏭 probable fuente no forestal ·
foco de calor de un incendio probable, detectado <span style="color:#d7191c">●</span> hace menos de 24 h
<span style="color:#fd8d3c">●</span> hace 24 a 72 h <span style="color:#fecc5c">●</span> hace más de 72 h ·
<span style="color:#a3a0ad">●</span> gris: foco sin evidencia de incendio forestal<br>
<b>Exposición</b>: <span style="border:1.5px dashed #4d4d4d">&nbsp;&nbsp;&nbsp;&nbsp;</span>
5 km alrededor de los incendios probables. <b>No es el avance del fuego</b>.<br>
<b>Red de salud</b>: <span style="color:#008f8f;font-weight:bold">━ ■</span> DIRESA
<span style="color:#008f8f;font-weight:bold">╌ ●</span> EsSalud (línea = a quién refiere) · borde
<span style="color:#b3261e;font-weight:bold">■</span> crítico <span style="color:#d98c00;font-weight:bold">■</span> atención<br>
Personas a ≤ 5 km: <b>{poblacion_total_5km:,}</b> de cualquier foco · <b>{poblacion_probables_5km:,}</b> de incendios
forestales probables (WorldPop 2026, sin doble conteo)<br>
<i>Grafo aproximado con RENIPRESS; el tiempo de traslado es modelado (Malaria Atlas Project).</i><br>
Corte {CORTE:%d-%m-%Y %H:%M} (Lima) · Fuentes: NASA FIRMS · ESA WorldCover · WorldPop · RENIPRESS · OSIPTEL · INEI
</details>"""
mapa.get_root().html.add_child(folium.Element(leyenda))

mapa.save(SALIDAS / "mapa_incendios_junin.html")
print(f"[12] mapa interactivo: {SALIDAS / 'mapa_incendios_junin.html'}")


# ============================================================
# PROCESO 13: Mapa estático para el README
# Insumo:    distritos_junin, zona_probables, focos_junin, red
# Resultado: salidas/mapa_junin.png
# ============================================================

fig, ax = plt.subplots(figsize=(8, 8.5))
# Coroplético: los distritos se colorean según la población expuesta a incendios probables
distritos_junin.plot(ax=ax, column="poblacion_expuesta", cmap="Purples", linewidth=0.3,
                     edgecolor="#777777", legend=True,
                     legend_kwds={"label": "Personas a ≤ 5 km de incendios forestales probables", "shrink": 0.5})
if len(zona_probables):
    zona_probables.boundary.plot(ax=ax, color="#4d4d4d", linewidth=0.6, linestyle="--")
probable = focos_junin["evidencia"].eq("Incendio forestal probable")
focos_junin[~probable].plot(ax=ax, color="#bdbdbd", markersize=2)
focos_junin[probable].plot(ax=ax, color="#b10026", markersize=3)
for estado_alerta, color in COLOR_ALERTA.items():
    seleccion = red[red["alerta"].eq(estado_alerta)]
    ax.scatter(seleccion["lon"], seleccion["lat"], s=26, marker="P", color=color, edgecolor="white",
               linewidth=0.4, label=f"Salud: {estado_alerta.lower()} ({len(seleccion)})", zorder=3)
ax.legend(frameon=False, loc="lower left", fontsize=8)
ax.set_title(f"Junín: incendios forestales probables ({focos_junin['fecha_hora_lima'].min():%d %b} – "
             f"{focos_junin['fecha_hora_lima'].max():%d %b %Y}; rojo) y otros focos (gris),\n"
             f"zona de 5 km y establecimientos de la red pública expuestos", fontsize=10)
ax.axis("off")
plt.tight_layout()
plt.savefig(SALIDAS / "mapa_junin.png", dpi=150)
print(f"[13] mapa estático: {SALIDAS / 'mapa_junin.png'}")


# ============================================================
# PROCESO 14: Página web de presentación
# Insumo:    capas y tablas de los procesos anteriores + web/plantilla_web.html
# Resultado: salidas/web_junin.html (página con cifras, mapa propio, criterios y red de salud)
# ============================================================
# La página dibuja su propio mapa, sin mapas de fondo de internet. Por eso las capas se
# exportan en METROS (UTM 18 Sur): así el navegador las dibuja tal cual, sin deformarlas.

def capa_web(gdf, columnas, tolerancia_m):
    """Exporta una capa a GeoJSON en metros, simplificada y redondeada al metro."""
    # .to_crs() -> metros | .simplify() -> menos vértices | .set_precision(1) -> coordenadas enteras
    #   Antes: polígonos detallados en grados | Después: polígonos livianos en metros
    capa = gdf[columnas + ["geometry"]].to_crs(epsg=CRS_METROS)
    capa["geometry"] = capa.simplify(tolerancia_m).set_precision(1)
    return json.loads(capa.to_json(drop_id=True))


def puntos_web(tabla, lon="lon", lat="lat"):
    """Devuelve las coordenadas x, y en metros (UTM 18 Sur) de una tabla con longitud y latitud."""
    puntos_xy = gpd.GeoSeries(gpd.points_from_xy(tabla[lon], tabla[lat]), crs=4326).to_crs(epsg=CRS_METROS)
    return puntos_xy.x.round().astype(int).tolist(), puntos_xy.y.round().astype(int).tolist()


def fila_web(fila):
    """Convierte una fila de una tabla en un diccionario con fechas como texto."""
    return a_texto(fila.to_frame().T).iloc[0].to_dict()


# --- Paso 14.1: Puntos (focos, incendios y nodos de la red) ---
NIVELES = ["Incendio forestal probable", "Probable quema agropecuaria", "Probable fuente no forestal",
           "Detección por verificar"]
fx, fy = puntos_web(focos_junin, "longitude", "latitude")
focos_web = [[x, y, round(h, 1), round(frp, 1), NIVELES.index(ev), inc, CLASES_WORLDCOVER.get(c, "Otro")]
             for x, y, h, frp, ev, inc, c in
             zip(fx, fy, focos_junin["horas_desde_foco"], focos_junin["frp"], focos_junin["evidencia"],
                 focos_junin["incendio"], focos_junin["clase_worldcover"])]

incendios_mapa = a_texto(tabla_incendios[tabla_incendios["evidencia"] != "Detección por verificar"]).reset_index(drop=True)
incendios_mapa["x"], incendios_mapa["y"] = puntos_web(incendios_mapa)
incendios_web = incendios_mapa[["incendio", "evidencia", "estado", "provincia", "distrito", "focos", "pasadas",
                                "satelites", "cobertura_principal", "inicio", "ultimo_foco", "huella_ha",
                                "poblacion_5km", "establecimientos_5km", "x", "y"]].to_dict("records")

red_web = red.copy()
red_web["x"], red_web["y"] = puntos_web(red_web)
red_web = (red_web[["COD_IPRESS", "NOMBRE", "sistema", "CATEGORIA", "rol", "RED", "MICRORED", "DISTRITO", "padre",
                    "hospital_ref", "nombre_hospital_ref", "min_a_padre", "min_a_hospital", "conectividad",
                    "cobertura_4g_pct", "alerta", "expuesto", "expuesto_reciente", "en_ruta_expuesta",
                    "dist_incendio_km", "incendio_cercano", "x", "y"]]
           .replace({np.nan: None}).to_dict("records"))

# --- Paso 14.2: Cifras clave, tablas y casos para los hallazgos ---
conteo_evidencia = incendios["evidencia"].value_counts()
cifras = {
    "corte": f"{CORTE:%d/%m/%Y %H:%M}",
    "desde": f"{focos_junin['fecha_hora_lima'].min():%d/%m}",
    "hasta": f"{focos_junin['fecha_hora_lima'].max():%d/%m/%Y}",
    "focos": n_focos, "incendios": n_incendios,
    "probables": len(incendios_probables), "recientes": len(incendios_recientes),
    "poblacion_total": poblacion_total_5km, "poblacion_probables": poblacion_probables_5km,
    "poblacion_recientes": poblacion_recientes_5km, "poblacion_junin": int(distritos_junin["poblacion_total"].sum()),
    "salud_cualquier": len(cerca_cualquier), "salud_probables": len(establecimientos_expuestos),
    "salud_recientes": len(salud_recientes), "salud_junin": len(establecimientos_junin),
    "distritos_alta": int((ficha["prioridad"] == "1 · Alta").sum()), "distritos": len(ficha),
    "red_nodos": len(red), "red_expuestos": int(red["expuesto"].sum()),
    "criticos": int(red["alerta"].eq("Crítico").sum()), "atencion": int(red["alerta"].eq("Atención").sum()),
    "red_sin_4g": int(red["conectividad"].ne("4G garantizado").sum()),
    "min_traslado_largo": MIN_TRASLADO_LARGO, "umbral_luz": UMBRAL_LUZ,
}
evidencia_web = [{"nivel": nivel, "incendios": int(conteo_evidencia.get(nivel, 0)),
                  "focos": int(focos_junin["evidencia"].eq(nivel).sum())} for nivel in NIVELES]
alta_web = (ficha[ficha["prioridad"] == "1 · Alta"].reset_index()
            [["ubigeo", "nombdist", "nombprov", "focos_72h", "incendios_prioritarios",
              "poblacion_expuesta", "establecimientos_expuestos", "puntos_criticos"]].to_dict("records"))
hospitales_web = por_hospital.reset_index().replace({np.nan: None}).to_dict("records")
probables_tabla = tabla_incendios[tabla_incendios["evidencia"] == "Incendio forestal probable"]
mayor = fila_web(tabla_incendios.loc[tabla_incendios["focos"].idxmax()])      # el incendio con más focos
descartado = (fila_web(tabla_incendios[tabla_incendios["evidencia"] == "Probable fuente no forestal"].iloc[0])
              if (tabla_incendios["evidencia"] == "Probable fuente no forestal").any() else None)
# El foco "por verificar" con más gente cerca (la tabla ya está ordenada por población)
aislado = fila_web(tabla_incendios[tabla_incendios["evidencia"] == "Detección por verificar"].iloc[0])
criticos_web = (red[red["alerta"].isin(["Crítico", "Atención"])]
                .assign(orden=red["alerta"].map(orden_alerta))
                .sort_values(["orden", "min_a_hospital"], ascending=[True, False])
                [["COD_IPRESS", "NOMBRE", "sistema", "CATEGORIA", "DISTRITO", "alerta", "nombre_hospital_ref", "min_a_hospital",
                  "conectividad", "cobertura_4g_pct", "dist_incendio_km", "incendio_cercano", "expuesto_reciente"]]
                .replace({np.nan: None}).to_dict("records"))

datos_web = {
    "cifras": cifras,
    "distritos": capa_web(distritos_junin, ["ubigeo", "nombdist", "nombprov", "prioridad", "focos_72h", "focos_7d",
                                            "poblacion_expuesta", "poblacion_total", "establecimientos_expuestos",
                                            "puntos_criticos"], 150),
    "provincias": capa_web(provincias, ["nombprov"], 150),
    "zona_probables": capa_web(zona_probables, [], 150) if len(zona_probables) else {"type": "FeatureCollection", "features": []},
    "zona_recientes": capa_web(zona_recientes, [], 100) if len(zona_recientes) else {"type": "FeatureCollection", "features": []},
    "focos": focos_web, "niveles": NIVELES, "incendios": incendios_web, "red": red_web,
    "evidencia": evidencia_web, "alta": alta_web, "hospitales": hospitales_web, "criticos": criticos_web,
    "casos": {"mayor": mayor, "descartado": descartado, "aislado": aislado},
}

# --- Paso 14.3: Insertar los datos en la plantilla de la página ---
# .replace(): pone el texto JSON donde la plantilla dice __DATOS_WEB__
#   Antes: plantilla sin datos | Después: una sola página HTML con todo adentro
plantilla = (CARPETA / "web" / "plantilla_web.html").read_text(encoding="utf-8")
texto_json = json.dumps(datos_web, ensure_ascii=False, separators=(",", ":"),
                        default=lambda v: v.item() if hasattr(v, "item") else str(v))
(SALIDAS / "web_junin.html").write_text(plantilla.replace("__DATOS_WEB__", texto_json), encoding="utf-8")
print(f"[14] página web: {SALIDAS / 'web_junin.html'} ({len(texto_json) / 1e6:.2f} MB de datos)")
