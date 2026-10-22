# ============================================================
# RED DE SALUD DE JUNÍN · grafo de referencia, tiempo de viaje y conectividad
#
# Capa BASE: se corre una vez al mes (cuando SUSALUD actualiza RENIPRESS), no cada
# vez que se actualizan los incendios. incendios_junin.py lee sus resultados.
#
# Pregunta: si un establecimiento pequeño recibe heridos, ¿a qué establecimiento de
# mayor nivel los refiere, cuánto se tarda en llegar y tiene señal para telemedicina?
#
# Cómo ejecutar:   python red_salud_junin.py
# Resultados:      salidas/red_salud_nodos.csv, salidas/red_salud_enlaces.csv,
#                  salidas/red_salud_junin.png
#
# AVISO: RENIPRESS trae la estructura administrativa (red y microrred), no los flujos
# oficiales de referencia y contrarreferencia, que define cada DIRESA o red asistencial.
# Este grafo es una APROXIMACIÓN con reglas explícitas (ver Proceso R2).
# ============================================================

# --- Paquetes ---
# Instalar (una sola vez por computadora):
#   uv pip install geopandas rasterio scikit-image scipy requests matplotlib
# Cargar (en cada ejecución): las líneas import de abajo.
import sys
import warnings
from pathlib import Path

import geopandas as gpd               # geopandas: tablas con geometría (puntos y líneas)
import matplotlib.pyplot as plt       # matplotlib: dibujo estático de la red para el README
import numpy as np                    # numpy: cálculos con arreglos de números
import pandas as pd                   # pandas: tablas, como Excel con código
import rasterio                       # rasterio: leer el raster de fricción (tiempo de viaje)
import requests                       # requests: descargar archivos
from rasterio.transform import rowcol  # rasterio: pasar de coordenadas a fila/columna del raster
from scipy.spatial import cKDTree     # scipy: buscar el centro poblado más cercano
from shapely.geometry import LineString  # shapely: dibujar cada enlace como una línea
from skimage.graph import MCP_Geometric  # scikit-image: camino de menor tiempo sobre un raster

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8")


# ============================================================
# PROCESO R0: Configuración
# Insumo:    nada (son decisiones del análisis)
# Resultado: rutas, direcciones web y reglas en un solo lugar
# ============================================================

CARPETA = Path(__file__).resolve().parent
DATOS = CARPETA / "datos"
SALIDAS = CARPETA / "salidas"
DATOS.mkdir(exist_ok=True)
SALIDAS.mkdir(exist_ok=True)

CRS_METROS = 32718                  # UTM zona 18 Sur, en metros
CAJA_JUNIN = (-77.0, -13.0, -73.0, -10.2)   # lon mín, lat mín, lon máx, lat máx (con margen)

# Establecimientos: RENIPRESS (SUSALUD), corte del 30-09-2026
URL_RENIPRESS = "https://www.datosabiertos.gob.pe/sites/default/files/RENIPRESS_30-09-2026.csv"
ARCHIVO_RENIPRESS = DATOS / "RENIPRESS_30-09-2026.csv"

# Fricción motorizada 2019 (Malaria Atlas Project; Weiss et al., 2020, Nature Medicine):
# minutos necesarios para recorrer 1 metro en cada celda de ~1 km, usando vías y terreno.
URL_FRICCION = ("https://data.malariaatlas.org/geoserver/Accessibility/ows?service=WCS&version=2.0.1"
                "&request=GetCoverage&coverageId=Accessibility__202001_Global_Motorized_Friction_Surface"
                "&subset=Lat(-13.0,-10.2)&subset=Long(-77.0,-73.0)&format=image/geotiff")
ARCHIVO_FRICCION = DATOS / "friccion_motorizada_junin.tif"

# Cobertura móvil por centro poblado (OSIPTEL, datos abiertos, 2025)
URL_OSIPTEL = ("https://www.datosabiertos.gob.pe/sites/default/files/Porcentaje%20de%20cobertura%20movil%20"
               "por%20centro%20poblado%20empresa%20operadora%20y%20tecnolog%C3%ADa.csv")
ARCHIVO_OSIPTEL = DATOS / "osiptel_cobertura_ccpp.csv"

# Reglas
ORDEN_CATEGORIAS = ["I-1", "I-2", "I-3", "I-4", "II-1", "II-2", "II-E", "III-1", "III-2", "III-E"]
INSTITUCIONES_DIRESA = ["GOBIERNO REGIONAL", "MINSA"]
# Establecimientos de una sola especialidad: no reciben quemados ni trauma, quedan fuera del árbol.
# Ojo: "IREN" no sirve como patrón porque también aparece en nombres como "CutivIRENI".
ESPECIALIZADOS = "MATERNO|NEOPL|SALUD MENTAL"
DIST_MAX_CCPP_M = 1_500                   # el centro poblado de OSIPTEL debe estar a <= 1.5 km
UMBRAL_4G = 50                            # 4G "garantizado" si cubre >= 50 % del área del centro poblado


def descargar(url, destino):
    """Descarga `url` en `destino` solo si el archivo todavía no existe."""
    if destino.exists():
        print(f"  ya estaba descargado: {destino.name}")
        return destino
    print(f"  descargando {destino.name} ...")
    respuesta = requests.get(url, timeout=900, headers={"User-Agent": "Mozilla/5.0"})
    respuesta.raise_for_status()
    destino.write_bytes(respuesta.content)
    return destino


# ============================================================
# PROCESO R1: Establecimientos de las dos redes públicas de Junín
# Insumo:    CSV RENIPRESS
# Resultado: nodos (un establecimiento por fila, con su sistema: DIRESA o EsSalud)
# ============================================================

descargar(URL_RENIPRESS, ARCHIVO_RENIPRESS)
renipress = pd.read_csv(ARCHIVO_RENIPRESS, sep=";", encoding="utf-8-sig", dtype=str)
renipress["lat"] = pd.to_numeric(renipress["NORTE"], errors="coerce")
renipress["lon"] = pd.to_numeric(renipress["ESTE"], errors="coerce")

# Filtros en cascada (se informa cuántos quedan en cada paso)
lon_min, lat_min, lon_max, lat_max = CAJA_JUNIN
filtros = [
    ("registros en RENIPRESS", pd.Series(True, index=renipress.index)),
    ("en Junín (UBIGEO 12...)", renipress["UBIGEO"].str.startswith("12")),
    ("estado ACTIVO", renipress["ESTADO"].eq("ACTIVO")),
    ("establecimiento de salud (no servicio de apoyo)", renipress["TIPO_ESTABLECIMIENTO"].str.startswith("ESTABLECIMIENTO")),
    ("categorizado (no '0')", renipress["CATEGORIA"].isin(ORDEN_CATEGORIAS)),
    ("de la DIRESA (Gob. Regional/MINSA) o EsSalud", renipress["INSTITUCION"].isin(INSTITUCIONES_DIRESA + ["ESSALUD"])),
]
conserva = pd.Series(True, index=renipress.index)
print("[R1] filtros de RENIPRESS:")
for descripcion, condicion in filtros:
    conserva &= condicion
    print(f"     {descripcion:50s} {conserva.sum():>6,}")
candidatos = renipress[conserva].copy()
con_coordenadas = candidatos["lat"].between(lat_min, lat_max) & candidatos["lon"].between(lon_min, lon_max)
print(f"     {'con coordenadas válidas (se pueden ubicar)':50s} {con_coordenadas.sum():>6,}")
print(f"     sin coordenadas (quedan fuera del grafo): {(~con_coordenadas).sum()}")

nodos = candidatos[con_coordenadas].reset_index(drop=True)
# Sistema: la red de la DIRESA (Gobierno Regional/MINSA) o la red asistencial de EsSalud
nodos["sistema"] = np.where(nodos["INSTITUCION"].eq("ESSALUD"), "EsSalud", "DIRESA Junín")
# Nivel numérico de la categoría (I-1 = 0 ... III-E = 9) para comparar "mayor" y "menor"
nodos["nivel"] = nodos["CATEGORIA"].map({c: i for i, c in enumerate(ORDEN_CATEGORIAS)})
print(f"     nodos por sistema: {nodos['sistema'].value_counts().to_dict()}")


# ============================================================
# PROCESO R2: Rol de cada establecimiento en la red
# Insumo:    nodos
# Resultado: columna "rol" (hospital regional, hospital, cabecera de microrred, establecimiento)
# ============================================================
# Reglas (explícitas para poder discutirlas):
#   - Hospital regional: categoría III, salvo los de una sola especialidad (materno, oncológico).
#   - Hospital: categoría II-1 o II-2.
#   - Cabecera de microrred (solo DIRESA): el de mayor categoría de su microrred; si empatan,
#     el más cercano al centro de la microrred.
#   - Establecimiento: todos los demás (puestos y centros de salud).
#   - Los especializados (II-E, III-E de una especialidad) quedan fuera del árbol.

nodos["especializado"] = (nodos["NOMBRE"].str.contains(ESPECIALIZADOS, case=False)
                          | nodos["CATEGORIA"].eq("II-E"))
nodos["rol"] = np.select(
    [nodos["especializado"],
     nodos["CATEGORIA"].str.startswith("III"),
     nodos["CATEGORIA"].isin(["II-1", "II-2"])],
    ["Especializado", "Hospital regional", "Hospital"], default="Establecimiento")

# Microrred válida: la que no dice "NO PERTENECE"
nodos["microrred_valida"] = (nodos["sistema"].eq("DIRESA Junín")
                             & ~nodos["MICRORED"].fillna("NO PERTENECE").str.contains("NO PERTENECE"))

# Cabecera: dentro de cada (red, microrred), el de mayor nivel; desempate por cercanía al centro
puntos_m = gpd.GeoSeries(gpd.points_from_xy(nodos["lon"], nodos["lat"]), crs=4326).to_crs(epsg=CRS_METROS)
nodos["x"], nodos["y"] = puntos_m.x.values, puntos_m.y.values
for (red, microrred), grupo in nodos[nodos["microrred_valida"] & nodos["rol"].eq("Establecimiento")].groupby(["RED", "MICRORED"]):
    mejores = grupo[grupo["nivel"] == grupo["nivel"].max()]
    centro_x, centro_y = grupo["x"].mean(), grupo["y"].mean()
    elegido = ((mejores["x"] - centro_x) ** 2 + (mejores["y"] - centro_y) ** 2).idxmin()
    nodos.loc[elegido, "rol"] = "Cabecera de microrred"

print(f"[R2] roles: {nodos.groupby(['sistema', 'rol']).size().to_dict()}")


# ============================================================
# PROCESO R3: Tiempo de viaje con el raster de fricción
# Insumo:    raster de fricción motorizada (minutos por metro) + nodos
# Resultado: una "superficie de tiempo" desde cada nodo que puede recibir referencias
# ============================================================
# Idea: el raster dice cuántos minutos cuesta cruzar cada celda (rápido en carretera, lento
# en selva o puna). El algoritmo de "costo mínimo" (MCP) busca, desde un hospital, el camino
# más rápido hacia todas las celdas. Es como un mapa de isócronas.

descargar(URL_FRICCION, ARCHIVO_FRICCION)
with rasterio.open(ARCHIVO_FRICCION) as fuente:
    friccion = fuente.read(1).astype("float64")
    transformacion = fuente.transform
    nodata_friccion = fuente.nodata
    lat_centro = (fuente.bounds.top + fuente.bounds.bottom) / 2
    alto_celda_m = fuente.res[1] * 110_600                               # 1° de latitud ≈ 110.6 km
    ancho_celda_m = fuente.res[0] * 111_320 * np.cos(np.radians(lat_centro))  # se achica con la latitud

# Celdas sin dato o con valor no válido: se vuelven muy lentas (no se cruzan)
invalidas = (friccion == nodata_friccion) | ~np.isfinite(friccion) | (friccion <= 0)
friccion[invalidas] = friccion[~invalidas].max() * 10
print(f"[R3] fricción: {friccion.shape} celdas de {ancho_celda_m:,.0f} × {alto_celda_m:,.0f} m "
      f"| celdas inválidas: {invalidas.sum()}")

# rowcol(): pasa la longitud y latitud de cada nodo a fila y columna del raster
filas, columnas = rowcol(transformacion, nodos["lon"].values, nodos["lat"].values)
nodos["fila"], nodos["columna"] = np.array(filas), np.array(columnas)

# MCP_Geometric(): prepara el cálculo; sampling = tamaño de la celda en metros (alto, ancho)
#   así el costo acumulado sale directamente en MINUTOS
calculadora = MCP_Geometric(friccion, sampling=(alto_celda_m, ancho_celda_m))
superficies = {}   # código IPRESS -> raster de minutos hasta ese nodo


def minutos_hasta(codigo_destino, fila, columna):
    """Minutos de viaje desde la celda (fila, columna) hasta el nodo destino."""
    if codigo_destino not in superficies:
        destino = nodos.loc[nodos["COD_IPRESS"] == codigo_destino].iloc[0]
        # find_costs(): desde el destino, el tiempo mínimo hacia todas las celdas
        #   Antes: raster de fricción | Después: raster de minutos hasta el destino
        # .copy() es obligatorio: find_costs() reutiliza su memoria interna y, sin copia,
        # cada nueva superficie sobrescribiría las anteriores.
        minutos, _ = calculadora.find_costs([(destino["fila"], destino["columna"])])
        superficies[codigo_destino] = minutos.copy()
    return float(superficies[codigo_destino][fila, columna])


def mas_cercano(fila_nodo, candidatos):
    """Entre los candidatos, el código del que se llega más rápido y los minutos."""
    tiempos = {c: minutos_hasta(c, fila_nodo["fila"], fila_nodo["columna"]) for c in candidatos["COD_IPRESS"]}
    codigo = min(tiempos, key=tiempos.get)
    return codigo, tiempos[codigo]


# ============================================================
# PROCESO R4: Enlaces de referencia (el grafo)
# Insumo:    nodos con rol + superficies de tiempo
# Resultado: columna "padre" (a quién refiere cada establecimiento) y minutos hasta él
# ============================================================
# DIRESA:  establecimiento -> cabecera de su microrred -> hospital de su red -> hospital regional
#          (si su red no tiene hospital, va al hospital más rápido de alcanzar)
# EsSalud: establecimiento -> hospital EsSalud más rápido de alcanzar -> hospital regional EsSalud

nodos["padre"] = None
nodos["min_a_padre"] = np.nan

for sistema in ["DIRESA Junín", "EsSalud"]:
    del_sistema = nodos[nodos["sistema"].eq(sistema)]
    regionales = del_sistema[del_sistema["rol"].eq("Hospital regional")]
    hospitales = del_sistema[del_sistema["rol"].eq("Hospital")]
    for i, nodo in del_sistema.iterrows():
        padre = None
        if nodo["rol"] == "Hospital":
            padre, minutos = mas_cercano(nodo, regionales)
        elif nodo["rol"] == "Cabecera de microrred":
            de_su_red = hospitales[hospitales["RED"] == nodo["RED"]]
            padre, minutos = mas_cercano(nodo, de_su_red if len(de_su_red) else pd.concat([hospitales, regionales]))
        elif nodo["rol"] == "Establecimiento":
            cabecera = del_sistema[del_sistema["rol"].eq("Cabecera de microrred")
                                   & del_sistema["RED"].eq(nodo["RED"])
                                   & del_sistema["MICRORED"].eq(nodo["MICRORED"])]
            if nodo["microrred_valida"] and len(cabecera):
                padre = cabecera["COD_IPRESS"].iat[0]
                minutos = minutos_hasta(padre, nodo["fila"], nodo["columna"])
            else:
                padre, minutos = mas_cercano(nodo, pd.concat([hospitales, regionales]))
        if padre is not None:
            nodos.loc[i, ["padre", "min_a_padre"]] = [padre, minutos]

# --- Hospital de referencia: el primer hospital (o regional) subiendo por el árbol ---
rol_de = nodos.set_index("COD_IPRESS")["rol"]
padre_de = nodos.set_index("COD_IPRESS")["padre"]


def hospital_de_referencia(codigo):
    """Sube por el árbol hasta encontrar un hospital o un hospital regional."""
    actual = padre_de.get(codigo)
    while actual is not None and rol_de.get(actual) not in ("Hospital", "Hospital regional"):
        actual = padre_de.get(actual)
    return actual


nodos["hospital_ref"] = [hospital_de_referencia(c) if r not in ("Hospital regional", "Especializado") else None
                         for c, r in zip(nodos["COD_IPRESS"], nodos["rol"])]
nodos["min_a_hospital"] = [minutos_hasta(h, f, c) if h else np.nan
                           for h, f, c in zip(nodos["hospital_ref"], nodos["fila"], nodos["columna"])]
nombre_de = nodos.set_index("COD_IPRESS")["NOMBRE"]
nodos["nombre_padre"] = nodos["padre"].map(nombre_de)
nodos["nombre_hospital_ref"] = nodos["hospital_ref"].map(nombre_de)
print(f"[R4] enlaces: {nodos['padre'].notna().sum()} | superficies de tiempo calculadas: {len(superficies)}")
print(f"     minutos al hospital de referencia (establecimientos): "
      f"{nodos.loc[nodos['rol'].eq('Establecimiento'), 'min_a_hospital'].describe()[['50%', 'max']].round(0).to_dict()}")


# ============================================================
# PROCESO R5: Conectividad para telemedicina (OSIPTEL)
# Insumo:    CSV de cobertura móvil por centro poblado + nodos
# Resultado: cobertura 4G garantizada del centro poblado más cercano a cada establecimiento
# ============================================================
# OSIPTEL solo lista centros poblados con cobertura declarada. Si no hay uno a <= 1.5 km,
# se asume "sin cobertura declarada". 4G es la condición mínima para una teleconsulta con video.

descargar(URL_OSIPTEL, ARCHIVO_OSIPTEL)
osiptel = pd.read_csv(ARCHIVO_OSIPTEL, sep=";", encoding="latin-1", dtype=str)
osiptel.columns = osiptel.columns.str.strip()          # los nombres traen espacios sobrantes
osiptel = osiptel[osiptel["Ubigeo"].str.startswith("12")].reset_index(drop=True)


def maximo_por_tecnologia(tabla, tecnologia):
    """Mayor % de cobertura garantizada (CG) entre todos los operadores para una tecnología."""
    columnas = [c for c in tabla.columns if f"_{tecnologia}_CG" in c and "CAR" not in c]
    # .str.replace("%", ""): "55%" -> "55"; pd.to_numeric: texto -> número
    valores = tabla[columnas].apply(lambda s: pd.to_numeric(s.str.replace("%", "").str.strip(), errors="coerce"))
    return valores.max(axis=1).fillna(0)


osiptel["cg_4g"] = maximo_por_tecnologia(osiptel, "4G")
osiptel["cg_voz"] = np.maximum(maximo_por_tecnologia(osiptel, "2G"), maximo_por_tecnologia(osiptel, "3G"))
ccpp_m = gpd.GeoSeries(gpd.points_from_xy(pd.to_numeric(osiptel["Longitud"]), pd.to_numeric(osiptel["Latitud"])),
                       crs=4326).to_crs(epsg=CRS_METROS)

# cKDTree(...).query(): para cada establecimiento, el centro poblado de OSIPTEL más cercano
distancia, indice = cKDTree(np.column_stack([ccpp_m.x, ccpp_m.y])).query(np.column_stack([nodos["x"], nodos["y"]]))
cercano = distancia <= DIST_MAX_CCPP_M
nodos["ccpp_osiptel"] = np.where(cercano, osiptel["CentroPoblado"].values[indice], "")
nodos["cobertura_4g_pct"] = np.where(cercano, osiptel["cg_4g"].values[indice], 0)
nodos["cobertura_voz_pct"] = np.where(cercano, osiptel["cg_voz"].values[indice], 0)
nodos["conectividad"] = np.select(
    [nodos["cobertura_4g_pct"] >= UMBRAL_4G, nodos["cobertura_4g_pct"] > 0, nodos["cobertura_voz_pct"] > 0],
    ["4G garantizado", "4G parcial", "Solo voz (2G/3G)"], default="Sin cobertura declarada")
print(f"[R5] conectividad: {nodos['conectividad'].value_counts().to_dict()}")


# ============================================================
# PROCESO R6: Guardar nodos y enlaces
# Insumo:    nodos
# Resultado: salidas/red_salud_nodos.csv y salidas/red_salud_enlaces.csv
# ============================================================

columnas_nodos = ["COD_IPRESS", "NOMBRE", "sistema", "INSTITUCION", "CATEGORIA", "nivel", "rol", "RED", "MICRORED",
                  "PROVINCIA", "DISTRITO", "UBIGEO", "lat", "lon", "padre", "nombre_padre", "min_a_padre",
                  "hospital_ref", "nombre_hospital_ref", "min_a_hospital", "ccpp_osiptel", "cobertura_4g_pct",
                  "cobertura_voz_pct", "conectividad"]
tabla_nodos = nodos[columnas_nodos].copy()
tabla_nodos[["min_a_padre", "min_a_hospital"]] = tabla_nodos[["min_a_padre", "min_a_hospital"]].round(0)
tabla_nodos.to_csv(SALIDAS / "red_salud_nodos.csv", index=False, encoding="utf-8-sig")

# Enlaces: una fila por par (establecimiento -> padre), con su línea para el mapa
coordenadas = nodos.set_index("COD_IPRESS")[["lon", "lat"]]
enlaces = nodos[nodos["padre"].notna()][["COD_IPRESS", "padre", "sistema", "rol", "min_a_padre"]].copy()
enlaces["geometry"] = [LineString([coordenadas.loc[a].values, coordenadas.loc[b].values])
                       for a, b in zip(enlaces["COD_IPRESS"], enlaces["padre"])]
enlaces = gpd.GeoDataFrame(enlaces, geometry="geometry", crs=4326)
enlaces.drop(columns="geometry").to_csv(SALIDAS / "red_salud_enlaces.csv", index=False, encoding="utf-8-sig")

resumen = (nodos[nodos["rol"].eq("Establecimiento")].groupby(["sistema", "nombre_hospital_ref"])
           .agg(establecimientos=("COD_IPRESS", "size"),
                minutos_mediana=("min_a_hospital", "median"),
                sin_4g=("conectividad", lambda s: (~s.eq("4G garantizado")).sum()))
           .sort_values("establecimientos", ascending=False))
print("[R6] establecimientos por hospital de referencia:")
print(resumen.round(0).to_string())


# ============================================================
# PROCESO R7: Dibujo de la red (imagen para el README)
# Insumo:    nodos y enlaces
# Resultado: salidas/red_salud_junin.png
# ============================================================

# Un solo color para la salud; el sistema se distingue por la forma (línea y marcador)
ESTILO_SISTEMA = {"DIRESA Junín": ("-", "s"), "EsSalud": ("--", "o")}
fig, ax = plt.subplots(figsize=(8, 8.5))
for sistema, (linea, marcador) in ESTILO_SISTEMA.items():
    enlaces[enlaces["sistema"] == sistema].plot(ax=ax, color="#008f8f", linewidth=0.6, alpha=0.6, linestyle=linea)
    del_sistema = nodos[nodos["sistema"] == sistema]
    tamano = del_sistema["rol"].map({"Establecimiento": 5, "Cabecera de microrred": 18,
                                     "Hospital": 55, "Hospital regional": 110, "Especializado": 30})
    ax.scatter(del_sistema["lon"], del_sistema["lat"], s=tamano, color="#008f8f", marker=marcador, edgecolor="white",
               linewidth=0.4, label=sistema, zorder=3)
ax.set_title("Red de referencia aproximada de los establecimientos públicos de Junín\n"
             "(puntos más grandes = mayor nivel; las líneas unen cada establecimiento con su referencia)",
             fontsize=10)
ax.legend(frameon=False, loc="lower left")
ax.set_aspect(1 / np.cos(np.radians(lat_centro)))
ax.axis("off")
plt.tight_layout()
plt.savefig(SALIDAS / "red_salud_junin.png", dpi=150)
print(f"[R7] imagen: {SALIDAS / 'red_salud_junin.png'}")
