# README

Histogramas modelo para probar instalación

* Cada histograma esta normalizado de forma tal que recoge el valor esperado de partículas en un detector volumétrico en un área de $1~ \rm m^2$ durante un tiempo de $1~ \rm h$.
* Hay dos casos que representan diferentes rigideces magnéticas.
	* h_alta es un ejemplo de ruta con alta tasa de dosis
	* h_baja es un ejemplo de ruta con baja tasa de dosis
* La estructura de los archivos es la siguiente.
	* Cada fila corresponde a una energía en el rango 100 keV - 1 TeV.
	* Las energías estan distribuídas de forma logarítmica con 20 divisiones por decada.
	* Cada columna tiene el valor esperado de partículas por cada bin de energía.
	* Las columnas estan identificadas siguiendo la nomenclatura GPS de GEANT4. Los nombres son autocontenidos en casi todos los casos. Los más comunes son:
		* gamma
		* e+
		* e-
		* mu+
		* mu-
		* pi0
		* pi+
		* pi-
		* neutron
		* proton
		* anti_proton

## Encabezado

| E\_GeV | gamma | e- | e+ | neutron | proton |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 0.0001 | 2269855 | 0 | 0 | 0 | 0 |
| 0.0001122018 | 2470446 | 0 | 0 | 0 | 0 |
| 0.0001258925 | 2686895 | 0 | 0 | 0 | 0 |
| 0.0001412538 | 2896712 | 0 | 0 | 0 | 0 |
| 0.0001584893 | 3113886 | 0 | 0 | 0 | 0 |
| 0.0001778279 | 3326657 | 0 | 0 | 0 | 0 |
| 0.0001995262 | 3524066 | 0 | 0 | 0 | 0 |
| 0.0002238721 | 3714898 | 11552 | 31 | 0 | 0 |
| 0.0002511886 | 3890205 | 22213 | 62 | 0 | 0 |
| 0.0002818383 | 4059691 | 27201 | 137 | 0 | 0 |
| 0.0003162278 | 4197511 | 33290 | 200 | 0 | 0 |
| 0.0003548134 | 4322831 | 40270 | 387 | 0 | 0 |
| 0.0003981072 | 4422880 | 48705 | 620 | 0 | 0 |
| 0.0004466836 | 4506310 | 58550 | 1021 | 0 | 0 |
| 0.0005011872 | 4845805 | 69704 | 1739 | 0 | 0 |
| 0.0005623413 | 4557295 | 82321 | 2618 | 0 | 0 |
| 0.0006309573 | 4596543 | 96051 | 3783 | 0 | 0 |
| 0.0007079458 | 4627534 | 109140 | 5312 | 0 | 0 |
| 0.0007943282 | 4651700 | 124875 | 7499 | 0 | 0 |
| 0.0008912509 | 4673587 | 140158 | 10308 | 0 | 0 |
| 0.001 | 4700841 | 155509 | 13614 | 0 | 0 |

## Uso
Los histogramas se copian en la carpeta Ejemplo/histogramas y se usan en la función 

```python
sampleFromHistogram(inFile='inFile',outFile='outFile',sampleSize=0.01)
```

para generar la muestra adaptada a la geometría del modelo y el submuestreo seleccionado.
