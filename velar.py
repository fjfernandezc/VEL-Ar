'''
Modelo Velocidades Argentinas
Version 2
2022
'''

import numpy as np
import argparse
import json

# constantes
pi    = 3.141592653589793
torad = pi / 180.0

# WGS84
a  = 6378137.0
b  = 6.356752314245179e+06
e2 = 0.006694379990141

# epoca del terremoto de Maule y de Illapel
tmaule   = 2.0101589e+03
tillapel = 2.0157096E+03

encabezado = '''
		-----------------------------------------
		 Instituto Geografico Nacional (IGN-Ar) 
		     The Ohio State University (OSU)      
		      Modelo Velocidades Argentinas
		              (VEL-Ar v2.0)                 
		-----------------------------------------
'''

columnas = '{:^42s} {:^8s} {:^8s} {:^42s} {:^42s}'.format('coordenadas iniciales', 'ep inicial', 'ep final', 'coordenadas finales', 'coordenadas finales p07')

# texto de ayuda
descripcion = '''
{}

Programa para convertir coordenadas de una epoca a otra
utilizando el Modelo Velocidades Argentinas (VEL-Ar) desde
la terminal o desde un archivo de texto con el siguiente
formato:

latitud longitud altura epoca_incial epoca_final
latitud longitud altura epoca_incial epoca_final
...     ...      ...    ...          ...
latitud longitud altura epoca_incial epoca_final

Los valores deben estar separados por espacios. Alternativamente
se pueden utilizar coordenadas cartesianas geocentricas.
'''.format(encabezado)


def ecef2lla(x, y, z):

	'''
	funcion para convertir coordenadas cartesianas
	geocentricas a coordenadas geodesicas
	'''

	ep = np.sqrt((np.power(a, 2) - np.power(b, 2)) / np.power(b, 2))
	p  = np.sqrt(np.power(x, 2) + np.power(y, 2))
	th = np.arctan2(a * z, b * p)
	lon = np.arctan2(y, x)
	lat = np.arctan2((z + np.power(ep, 2) * b * np.power(np.sin(th), 3)), (p - e2 * a * np.power(np.cos(th), 3)))
	N = a / np.sqrt(1 - e2 * np.power(np.sin(lat), 2))
	h = p / np.cos(lat) - N

	# devuelve el valor de latitud en el intervalo [0,2*pi)
	lon = np.fmod(lon, 2 * pi) / torad
	lat = lat / torad

	return lat, lon, h


def lla2ecef(lat, lon, alt):

	'''
	funcion para convertir coordenadas geodesicas a
	coordenadas cartesianas geocentricas
	'''

	N = a / np.sqrt(1 - e2 * np.power(np.sin(lat * torad), 2))

	x = (N + alt) * np.cos(lat * torad) * np.cos(lon * torad)
	y = (N + alt) * np.cos(lat * torad) * np.sin(lon * torad)
	z = ((1 - e2) * N + alt) * np.sin(lat * torad)

	return x, y, z


def lg2ct(north, east, up, lat, lon):

	'''
	funcion para convertir desplazamientos desde un sistema
	topocentrico a un sistema cartesiano geocentrico
	'''

	R = np.zeros((3, 3))
	C = np.zeros((3, 1))

	R[0, 0] = -np.sin(lat * torad) * np.cos(lon * torad)
	R[1, 0] = -np.sin(lat * torad) * np.sin(lon * torad)
	R[2, 0] =  np.cos(lat * torad)
	R[0, 1] = -np.sin(lon * torad)
	R[1, 1] =  np.cos(lon * torad)
	R[2, 1] =  0
	R[0, 2] =  np.cos(lat * torad) * np.cos(lon * torad)
	R[1, 2] =  np.cos(lat * torad) * np.sin(lon * torad)
	R[2, 2] =  np.sin(lat * torad)

	C[0, 0] = north
	C[1, 0] = east
	C[2, 0] = up

	L = R.dot(C)

	dx = L[0, 0]
	dy = L[1, 0]
	dz = L[2, 0]

	return dx, dy, dz


def isecef(lat, lon, h):

	'''
	funcion para determinar si un juego de coordenadas es
	geodesico o cartesiano geocentrico
	'''

	R = np.sqrt(np.power(lat, 2) + np.power(lon, 2) + np.power(h, 2))

	if (R > 6e+06):
		return True
	else:
		return False


def lle2utm(lat, lon, lcm):

	'''
	funcion para convertir coordenadas geodesicas a planas UTM
	'''

	# utilizo WGS84
	f = 1 - np.sqrt(1 - e2)

	lat = lat * torad
	lon = lon * torad
	lcm = lcm * torad

	ko = 0.9996    # factor de escala
	No = 1e7       # falso norte
	Eo = 0         # falso este

	lam = lon - lcm
	  
	RN = a / np.power(1-e2 * np.power(np.sin(lat), 2), 0.5)
	RM = a * (1 - e2) / np.power(1 - e2 * np.power(np.sin(lat), 2), 1.5)
	h2 = e2 * np.power(np.cos(lat), 2) / (1 - e2)
	t = np.tan(lat)
	n = f / (2 - f)

	#----- Helmert (1880) expansion & simplification of Bessel series (faster)
	A0 = 1 + np.power(n, 2) / 4 + np.power(n, 4) / 64
	A2 = 3 / 2 * (n - np.power(n, 3) / 8)
	A4 = 15 / 16 * (np.power(n, 2) - np.power(n, 4) / 4)
	A6 = 35 / 48 * np.power(n, 3)
	A8 = 315 / 512 * np.power(n, 4)
	S = a / (1 + n) * (A0 * lat - A2 * np.sin(2 * lat) + A4 * np.sin(4 * lat) - A6 * np.sin(6 * lat) + A8 * np.sin(8 * lat))


	E1 = lam * np.cos(lat)
	E2 = np.power(lam, 3) * np.power(np.cos(lat), 3) / 6 * (1 - np.power(t, 2) + h2)
	E3 = np.power(lam, 5) * np.power(np.cos(lat), 5) / 120 * (5 - 18 * np.power(t, 2) + np.power(t, 4) + 14 * h2 - 58 * np.power(t, 2) * h2 + 13 * np.power(h2, 2) + 4 * np.power(h2, 3) - 64 * np.power(t, 2) * np.power(h2, 2) - 24 * np.power(t, 2) * np.power(h2, 3))
	E4 = np.power(lam, 7) * np.power(np.cos(lat), 7) / 5040 * (61 - 479 * np.power(t, 2) + 179 * np.power(t, 4) - np.power(t, 6))
	E = Eo + ko * RN * (E1 + E2 + E3 + E4)

	N1 = S / RN
	N2 = np.power(lam, 2) / 2 * np.sin(lat) * np.cos(lat)
	N3 = np.power(lam, 4) / 24 * np.sin(lat) * np.power(np.cos(lat), 3) * (5 - np.power(t, 2) + 9 * h2 + 4 * np.power(h2, 2))
	N4 = np.power(lam, 6) / 720 * np.sin(lat) * np.power(np.cos(lat), 5) * (61 - 58 * np.power(t, 2) + np.power(t, 4) +  270 * h2 - 330 * np.power(t, 2) * h2 + 445 * np.power(h2, 2) + 324 * np.power(h2, 3) - 680 * np.power(t, 2) * np.power(h2, 2) +  88 * np.power(h2, 4) - 600 * np.power(t, 2) * np.power(h2, 3) - 192 * np.power(t, 2) * np.power(h2, 4))
	N5 = np.power(lam, 8) / 40320 * np.sin(lat) * np.power(np.cos(lat), 7) * (1385 - 311 * np.power(t, 2) + 543 * np.power(t, 4) - np.power(t, 6))
	N = No + ko * RN * (N1 + N2 + N3 + N4 + N5)

	return np.array([N / 1000 - 6200, E / 1000])


def helmert7(x, y, z, marco_inicial='P07b'):

	'''
	Transformacion de 7 parametros de coordenadas geocentricas
	x: coordenada cartesiana X (m)
	y: coordenada cartesiana Y (m)
	z: coordenada cartesiana Z (m)
	return: coordenadas cartesianas transformadas
	'''

	# vector x y z original
	xyz_before = np.array([[x],
						   [y],
						   [z]])

	'''
	Parametros de transformacion de P07b a P07:
	T (mm) : -1.03 -1.00 -0.62
	R (ppb): -0.72  0.44 -0.05
	S (ppb):  1.8101
	'''

	tx = -1.03 / 1000
	ty = -1.00 / 1000
	tz = -0.62 / 1000
	rx = -0.72
	ry =  0.44
	rz = -0.05
	sc = 1.8101

	# vector de parametros
	parameters = np.array([[tx],
						   [ty],
						   [tz],
						   [sc],
						   [rx],
						   [ry],
						   [rz]])
	
	if marco_inicial == 'P07':
		parameters *= -1

	# matriz de diseno
	A = np.array([[1.0, 0.0, 0.0, x*1e-9,      0.0,  -z*1e-9,   y*1e-9],
				  [0.0, 1.0, 0.0, y*1e-9,   z*1e-9,      0.0,  -x*1e-9],
				  [0.0, 0.0, 1.0, z*1e-9,  -y*1e-9,   x*1e-9,      0.0]])

	# vector x y z transformado
	xyz_after = A.dot(parameters) + xyz_before

	xtrans = float(xyz_after[0])
	ytrans = float(xyz_after[1])
	ztrans = float(xyz_after[2])

	return xtrans, ytrans, ztrans


def get_variable(lat, lon, archivo, max_dist=None):

	'''
	funcion para interpolar las grillas del modelo.
	el formato de la grilla debe ser:

	latitud longitud coeficiente_norte coeficiente_este
	
	los valores deben estar separados por espacios
	'''

	# radio terrestre
	R = 6371.0

	# carga los archivos del modelo
	grilla = np.loadtxt(archivo)

	# pasa los valores a radianes
	grilla_rad = grilla * torad

	# busco los puntos mas cercanos a los parametros ingresados
	nfilas = grilla.shape[0]
	pos = np.ones((nfilas, 2))
	pos[:, 0] = pos[:, 0] * lat * torad
	pos[:, 1] = pos[:, 1] * lon * torad
	
	# distancia sobre circulo mayor utilizando la formula estable de haversine
	# https://en.wikipedia.org/wiki/Haversine_formula
	dist = np.zeros((nfilas, 1))
	d_phi = (grilla_rad[:, 0] - pos[:, 0]) / 2
	d_lam = (grilla_rad[:, 1] - pos[:, 1]) / 2
	dist = 2 * R * np.arcsin(np.sqrt(np.power(np.sin(d_phi), 2) + np.cos(grilla_rad[:, 0]) * np.cos(pos[:, 0]) * np.power(np.sin(d_lam), 2)))

	# ordeno las distancias entre el punto en estudio y los de las grillas de menor a mayor
	indices = np.argsort(dist)
	inter = np.concatenate((dist[indices].reshape((nfilas, 1)), grilla[indices, :]), axis=1)

	if max_dist is not None:
		# si hay una distancia maxima mas alla de la cual no se interpola
		if inter[0, 0] > max_dist:
			# si la distancia del punto mas cercano de la grilla al punto en estudio es mayor
			# a la distancia maxima especificada, no se interpola
			return np.array([0.0, 0.0])

	# me quedo con 4 elementos mas cercanos
	points = 4
	inter = inter[0: points, :]

	# para cada uno de los puntos mas cercanos convierto las coordenadas geodesicas a planas utm
	A = np.ones((points, 3))
	for i in range(points):
		tNE = lle2utm(inter[i, 1], inter[i, 2], lon)
		A[i, 1] = tNE[0]
		A[i, 2] = tNE[1]

	# convierto las coordenadas del punto en estudio a planas utm
	tNE = lle2utm(lat, lon ,lon)
	a = np.array([1.0, tNE[0], tNE[1]])

	# coeficientes de la grilla para la componente norte
	L1 = inter[:, 3]
	# coeficientes de la grilla para la componente este
	L2 = inter[:, 4]

	# resuelve un sistema lineal de ecuaciones para cada componente
	X1 = np.linalg.solve((A.T).dot(A), (A.T).dot(L1))
	X2 = np.linalg.solve((A.T).dot(A), (A.T).dot(L2))

	# devuelve el resultado de la interpolacion en norte y este
	res = np.array([a.dot(X1), a.dot(X2)])

	return res


def heaviside(ti, tf, ts):

	'''
	funcion escalon para determinar el sentido del salto cosismico
	'''

	if ((tf - ts) > 0.0) and ((ti - ts) < 0.0):
		# la fecha objetivo es despues del sismo y la fecha de partida es antes del sismo
		return 1.0
	elif ((tf - ts) < 0.0) and ((ti - ts) > 0.0):
		# la fecha objetivo es antes del sismo y la fecha de partida es despues del sismo
		return -1.0
	else:
		# ambas fechas son anteriores o posteriores al sismo
		return 0.0


def delta_log(ti, tf, ts):

	'''
	funcion para determinar el intervalo de tiempo que se ingresa
	dentro del logaritmo para calcular la componente postsismica
	'''

	if (ti - ts) >= 0:
		# la fecha de partida es despues del sismo
		return ti - ts
	elif (ti - ts) < 0:
		# la fecha de partida es anterior al sismo
		return tf - ts


def velar(lat, lon, alt, ti, tf):

	'''
	funcion principal para determinar los parametros del modelo para el punto
	en analisis y computar las coordenadas en la epoca objetivo
	'''

	# diccionario para guardar los resultados
	res = {}
	res['datos_entrada'] = {
		'coordenadas'  : [lat, lon, alt],
		'epoca_inicial': ti,
		'epoca_final'  : tf
	}

	# si la fecha inicial es 2006.632 aplica una transformacion de 7 parametros P07 -> P07b
	if ti == 2006.632:
		# si las coordenadas de entrada no estan en ecef hay que convertirlas
		if not isecef(lat, lon, alt):
			x_p07, y_p07, z_p07 = lla2ecef(lat, lon, alt)
			x_p07b, y_p07b, z_p07b = helmert7(x_p07, y_p07, z_p07, marco_inicial='P07')
			lat, lon, alt = ecef2lla(x_p07b, y_p07b, z_p07b)
		else:
			x_p07, y_p07, z_p07 = lat, lon, alt
			x_p07b, y_p07b, z_p07b = helmert7(x_p07, y_p07, z_p07, marco_inicial='P07')
			lat, lon, alt = x_p07b, y_p07b, z_p07b	

	# si las coordenadas especificadas son cartesianas geocentricas
	# hay que convertirlas a geodesicas para interpolar
	xyz = False
	if isecef(lat, lon, alt):
		lat, lon, alt = ecef2lla(lat, lon, alt)
		xyz = True

	# ------------------- sumamos la componente intersismica siempre ------------------------------

	# interpolo
	intersismica = get_variable(lat, lon, 'grids/vel-ar-lin.txt')
	# sumo los desplazamientos en norte y este
	desplazamientos = intersismica * (tf - ti)

	# agrego al vector de resultados
	if xyz:
		x, y, z = lg2ct(intersismica[0], intersismica[1], 0, lat, lon)
		res['componente_intersismica'] = [x, y, z]
	else:
		res['componente_intersismica'] = [intersismica[0], intersismica[1], 0.0]

	# ------------------- sumamos la influencia del sismo de maule --------------------------------

	# interpolo
	if lat < -25.0 and lat > -45.0 and lon < -55.0 and lon > -75.0:
		# limites de la grilla para el salto cosismico de maule
		cosismica = get_variable(lat, lon, 'grids/vel-ar-cos-maule.txt')
	else:
		cosismica = np.array([0.0, 0.0])
	# si el punto mas cercano esta a mas de 13.5 km no interpola
	postsismica = get_variable(lat, lon, 'grids/vel-ar-log-maule.txt', max_dist=13.5)
	
	# agrego al vector de resultados
	if xyz:
		x, y, z = lg2ct(cosismica[0], cosismica[1], 0, lat, lon)
		res['salto_cosismico_maule'] = [x, y, z]
		x, y, z = lg2ct(postsismica[0], postsismica[1], 0, lat, lon)
		res['componente_postsismica_maule'] = [x, y, z]
	else:
		res['salto_cosismico_maule'] = [cosismica[0], cosismica[1], 0.0]
		res['componente_postsismica_maule'] = [postsismica[0], postsismica[1], 0.0]

	# h determina si se suma (vamos hacia adelante), se resta (vamos hacia atras) 
	# o no pasamos por el sismo
	h = heaviside(ti, tf, tmaule)
	
	# componente cosismica de maule
	desplazamientos += h * cosismica

	# componente postsismica cuando pasamos por el sismo (si ambas fechas son anteriores
	# o posteriores h es 0 entonces no se suma)
	if h != 0.0:
		desplazamientos += h * postsismica * np.log10(1 + (delta_log(ti, tf, tmaule) / 0.5))

	# si las dos fechas son posteriores al sismo sumo la componente postsismica
	if (h == 0.0) and (ti > tmaule):
		desplazamientos += postsismica * np.log10(1 + ((tf - tmaule) / 0.5)) - postsismica * np.log10(1 + ((ti - tmaule) / 0.5))

	# ------------------- sumamos la influencia del sismo de illapel ------------------------------

	# interpolo
	cosismica   = get_variable(lat, lon, 'grids/vel-ar-cos-illapel.txt')
	# si el punto mas cercano esta a mas de 12 km no interpola
	postsismica = get_variable(lat, lon, 'grids/vel-ar-log-illapel.txt', max_dist=12.0)

	# agrego al vector de resultados
	if xyz:
		x, y, z = lg2ct(cosismica[0], cosismica[1], 0, lat, lon)
		res['salto_cosismico_illapel'] = [x, y, z]
		x, y, z = lg2ct(postsismica[0], postsismica[1], 0, lat, lon)
		res['componente_postsismica_illapel'] = [x, y, z]
	else:
		res['salto_cosismico_illapel'] = [cosismica[0], cosismica[1], 0.0]
		res['componente_postsismica_illapel'] = [postsismica[0], postsismica[1], 0.0]

	# h determina si se suma (vamos hacia adelante), se resta (vamos hacia atras) 
	# o no pasamos por el sismo
	h = heaviside(ti, tf, tillapel)
	
	# componente cosismica de illapel
	desplazamientos += h * cosismica

	# componente postsismica cuando pasamos por el sismo (si ambas fechas son anteriores
	# o posteriores h es 0 entonces no se suma)
	if h != 0.0:
		desplazamientos += h * postsismica * np.log10(1 + (delta_log(ti, tf, tillapel) / 0.5))

	# si las dos fechas son posteriores al sismo sumo la componente postsismica
	if (h == 0.0) and (ti > tillapel):
		desplazamientos += postsismica * np.log10(1 + ((tf - tillapel) / 0.5)) - postsismica * np.log10(1 + ((ti - tillapel) / 0.5))

	# ------------------- obtenemos las coordenadas finales ---------------------------------------

	# convierto los valores de norte este a lat lon
	v = a / np.sqrt(1 - e2 * np.power(np.sin(lat * torad) , 2))
	r = v * (1 - e2) / (1 - e2 * np.power(np.sin(lat * torad), 2))
	dlat = (desplazamientos[0] / (r + alt)) / torad
	dlon = (desplazamientos[1] / np.cos(lat * torad) / (v + alt)) / torad

	# agrego al vector de resultados las coordenadas transformadas
	if xyz:
		x, y, z = lla2ecef(lat + dlat, lon + dlon, alt)
		res['coordenadas_transformadas'] = [x, y, z]
		x, y, z = lg2ct(desplazamientos[0], desplazamientos[1], 0, lat, lon)
		res['desplazamientos_aplicados'] = [x, y, z]
	else:
		res['coordenadas_transformadas'] = [lat + dlat, lon + dlon, alt]
		res['desplazamientos_aplicados'] = [desplazamientos[0], desplazamientos[1], 0.0]
	
	# si la fecha final es 2006.632 aplica una transformacion de 7 parametros P07b -> P07
	if tf == 2006.632:
		if xyz:
			x_p07b, y_p07b, z_p07b = res['coordenadas_transformadas']
			x_p07 , y_p07 , z_p07  = helmert7(x_p07b, y_p07b, z_p07b)
			res['coordenadas_transformadas_p07'] = [x_p07, y_p07, z_p07]
		else:
			lat_p07b, lon_p07b, alt_p07b = res['coordenadas_transformadas']
			x_p07b  , y_p07b  , z_p07b   = lla2ecef(lat_p07b, lon_p07b, alt_p07b)
			x_p07   , y_p07   , z_p07    = helmert7(x_p07b, y_p07b, z_p07b)
			lat_p07 , lon_p07 , alt_p07  = ecef2lla(x_p07, y_p07, z_p07)
			res['coordenadas_transformadas_p07'] = [lat_p07, lon_p07, alt_p07]

	return res


def imprime_resultados(res):

	xyz = isecef(res['datos_entrada']['coordenadas'][0], res['datos_entrada']['coordenadas'][1], res['datos_entrada']['coordenadas'][2])
	if xyz:
		if 'coordenadas_transformadas_p07' in res.keys():
			formato = '{:>+13.4f} {:>+13.4f} {:>+13.4f} {:>+9.4f} {:>+9.4f} {:>+13.4f} {:>+13.4f} {:>+13.4f} {:>+13.4f} {:>+13.4f} {:>+13.4f}\n'
		else:
			formato = '{:>+13.4f} {:>+13.4f} {:>+13.4f} {:>+9.4f} {:>+9.4f} {:>+13.4f} {:>+13.4f} {:>+13.4f}\n'
	else:
		if 'coordenadas_transformadas_p07' in res.keys():
			formato = '{:>+12.10f} {:>+12.10f} {:>+9.4f} {:>+9.4f} {:>+9.4f} {:>+12.10f} {:>+12.10f} {:>+9.4f} {:>+12.10f} {:>+12.10f} {:>+9.4f}\n'
		else:
			formato = '{:>+12.10f} {:>+12.10f} {:>+9.4f} {:>+9.4f} {:>+9.4f} {:>+12.10f} {:>+12.10f} {:>+9.4f}\n'

	if 'coordenadas_transformadas_p07' in res.keys():
		return formato.format(*res['datos_entrada']['coordenadas'], res['datos_entrada']['epoca_inicial'], res['datos_entrada']['epoca_final'], *res['coordenadas_transformadas'], *res['coordenadas_transformadas_p07'])
	else:
		return formato.format(*res['datos_entrada']['coordenadas'], res['datos_entrada']['epoca_inicial'], res['datos_entrada']['epoca_final'], *res['coordenadas_transformadas'])


def main():


	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=descripcion)

	grupo1 = parser.add_argument_group('para transformacion de coordenadas desde terminal')
	grupo2 = parser.add_argument_group('para transformacion de coordenadas desde un archivo')

	grupo1.add_argument('-crd', type=float, nargs=3, metavar='crd', help='especificar las coordenadas a transformar')
	grupo1.add_argument('-ti', type=float, nargs=1, metavar='epoca', help='especificar la epoca inicial')
	grupo1.add_argument('-tf', type=float, nargs=1, metavar='epoca', help='especificar la epoca final')
	grupo2.add_argument('-f', type=str, nargs=1, metavar='archivo', help='especificar archivo de coordenadas')
	grupo2.add_argument('-o', type=str, nargs=1, metavar='archivo', help='especificar archivo de resultados')
	grupo2.add_argument('-json', action='store_true', help='guardar los resultados en un json')

	args = parser.parse_args()

	# planteamos dos opciones

	if args.crd:

		# 1) convertir las coordenadas de un unico punto a partir de los parametros pasados por linea de comandos
		if (not args.ti) or (not args.tf):
			parser.print_help()
			return

		# aplico velar
		res = velar(args.crd[0], args.crd[1], args.crd[2], args.ti[0], args.tf[0])
		
		# imprimo los resultados
		print(encabezado)
		print(columnas)
		print(imprime_resultados(res), end='')

	elif args.f:

		# 2) convertir las coordenadas de uno o varios puntos especificados en un archivo de texto y guardar
		#    los resultados en un archivo de texto (solo los datos de entrada y las coordenadas de salida) o
		#    en un json (todos los datos, incluidos los saltos y coeficientes aplicado en cada punto)

		if not args.o:
			parser.print_help()
			return

		print(encabezado)
		print(columnas)

		try:
			puntos = np.loadtxt(args.f[0], dtype=np.float64)
		except OSError:
			exit('El archivo de puntos {} no existe'.format(args.f[0]))
		except ValueError:
			exit('Verificar el formato del archivo de puntos {}'.format(args.f[0]))

		if len(puntos.shape) == 1:
			# el archivo contiene un solo punto
			puntos = puntos.reshape((-1, 5))
			
		resultados = {}
		# aplico velar para cada punto
		for i in range(puntos.shape[0]):
			lat = puntos[i, 0]
			lon = puntos[i, 1]
			alt = puntos[i, 2]
			ti  = puntos[i, 3]
			tf  = puntos[i, 4]
			res = velar(lat, lon, alt, ti, tf)
			resultados['punto {}'.format(i + 1)] = res
			print(imprime_resultados(res), end='')

		if args.json:
			# guardo en un json con todos los parametros
			with open(args.o[0], 'w') as f:
				json.dump(resultados, f, indent=4)
		else:
			# guardo en un txt con los datos de entrada y las coordenadas unicamente
			with open(args.o[0], 'w') as f:
				for punto in resultados:
					f.write(imprime_resultados(resultados[punto]))

	else:
		parser.print_help()



if __name__ == '__main__':

	main()

	
