# VEL-Ar

Durante los días 27 de febrero de 2010 y 16 de septiembre de 2015 ocurrieron sismos de magnitud 8.8 en la región de Maule (República de Chile) y de 8.3 en la región de Coquimbo (República de Chile) respectivamente, que provocaron un desplazamiento co-sísmico de la corteza terrestre entre las latitudes 25º S y 40º S, con valores superiores a los cinco metros en la zona del epicentro (para el caso de Maule) y de aproximadamente dos centímetros en la costa Atlántica. Como consecuencia de estos sismos, las coordenadas de los puntos próximos a los epicentros dejaron de evolucionar linealmente, dificultando en el territorio nacional la determinación de coordenadas precisas en el marco de referencia [POSGAR07] (definido para la época 2006.632) en las zonas más afectadas por los sismos (principalmente en las provincias de San Juan, Mendoza y Neuquén).

El Modelo de Predicción de Trayectorias (MPT) de la República Argentina, denominado VEL-Ar, permite la transformación temporal de coordenadas del marco de referencia operacional del IGN, basado actualmente en IGS14 (denominado [POSGAR07b]), y facilita el acceso de los usuarios al marco de referencia oficial POSGAR07 para la época 2006.632.

Los MPT surgen con la necesidad de obtener predicciones de las posiciones de vértices geodésicos pasivos (es decir, aquellos sin monitoreo GNSS continuo) para evaluar la deformación relativa de dichos puntos respecto de un juego de coordenadas conocidas en la época convencional de POSGAR07 (2006.632). La aplicación de VEL-Ar y POSGAR07b como herramientas para el acceso al marco de referencia oficial POSGAR07 convierten a la República Argentina en el primer país de las Américas en implementar un marco de referencia cinemático (MRC). A diferencia de los marcos de referencia estáticos (en los cuales las coordenadas no cambian con el tiempo), los MRC observan y modelan fenómenos geofísicos tales como la deformación por interacción entre bordes de placas, los saltos instantáneos debido a sismos (cosísmicos), la deformación posterior a los sismos (postsísmica), y la deformación vertical y horizontal debido a cambios en la hidrósfera y criósfera. Como resultado de la aplicación de un MRC, VEL-Ar y POSGAR07b aseguran el acceso a POSGAR07 con una incertidumbre inferior a 0.05 m en casi la totalidad del territorio nacional.

VEL-Ar es un desarrollo conjunto del IGN, Ohio State University (EEUU) y The University of Memphis (EEUU), y es calculado utilizando los parámetros obtenidos de los Modelos Extendidos de Trayectorias (MET) de las series de tiempo GNSS de [RAMSAC]. La última versión disponible de VEL-Ar (2022) es denominada VEL-Ar v2.0, y contempla los efectos horizontales de los siguientes fenómenos geofísicos:

1. las velocidades intersísmicas de la placa sudamericana (vinculadas a IGS14);
2. el salto cosísmico producto del sismo del 27 de febrero de 2010 en Maule, Chile;
3. las variaciones no lineales producidas por efectos postsísmicos debido al sismo de Maule;
4. el salto cosísmico producto del sismo del 16 de septiembre de 2015 en Illapel, Chile;
5. las variaciones no lineales producidas por efectos postsísmicos debido al sismo de Illapel.

## Dependencias

- Python3 (testeado en Python 3.8)
- NumPy (testeado en Numpy 1.22.4)


[POSGAR07]: <https://www.ign.gob.ar/NuestrasActividades/Geodesia/Posgar07>
[POSGAR07b]: <https://www.ign.gob.ar/NuestrasActividades/Geodesia/Posgar07>
[RAMSAC]: <https://www.ign.gob.ar/NuestrasActividades/Geodesia/Ramsac>
