# Metadatos del Dataset

Este dataset proviene del estudio sobre **2018-MetabotypingPaper** y contiene datos de metabolitos obtenidos a partir de muestras de sujetos. A continuación, se describen los metadatos y las variables que componen este conjunto de datos.

## Descripción General

El conjunto de datos contiene información sobre los metabolitos de 39 muestras de sujetos, así como metadatos asociados a cada muestra, como la cirugía realizada, edad, género, grupo de pertenencia, y varios otros parámetros relacionados con la salud. El dataset está compuesto por dos archivos:

1. `DataValues_S013.csv`: Datos de los metabolitos.
2. `DataInfo_S013.csv`: Metadatos asociados a cada muestra.

## Variables en los Metadatos

A continuación, se detallan las variables presentes en los metadatos del dataset summarized_experiment:

| Variable       | Descripción                                                      |
|----------------|------------------------------------------------------------------|
| `SUBJECTS`     | Identificador único del sujeto (número de identificación).       |
| `SURGERY`      | Tipo de cirugía realizada al sujeto (ByPass, tubular).           |
| `AGE`          | Edad (en años).                                                  |
| `GENDER`       | Género del sujeto (Masculino, Femenino).                         |
| `GROUP`        | Grupo al que pertenece el sujeto en el estudio (Grupo1, Grupo2). |
| `MEDMM_T0`     | Datos iniciales con valores categóricos.                         |
| `MEDCOL_T0`    | Datos iniciales con valores categóricos.                         |
| `MEDINF_T0`    | Datos iniciales con valores categóricos.                         |
| `MEDHTA_T0`    | Datos iniciales con valores categóricos.                         |

## Descripción de los Datos

### `DataValues_S013.csv`

Este archivo contiene los datos de los metabolitos para cada uno de los 39 sujetos. Los metabolitos están organizados en filas y las muestras en columnas. Las primeras 9 filas son utilizadas para los metadatos.

- **Dimensiones**: 686 metabolitos (filas) y 39 muestras (columnas).
- **Formato**: Los valores están almacenados en formato numérico.

### `DataInfo_S013.csv`

Este archivo contiene los metadatos de los sujetos, como se describe anteriormente. Es utilizado para asociar información adicional a cada muestra en el conjunto de datos de metabolitos.

- **Dimensiones**: 39 muestras (filas) y 9 variables (columnas).
- **Formato**: Los valores son principalmente categóricos y numéricos (este último en SUBJECTS, AGE, MEDMM_T0, MEDCOL_T0, MEDINF_T0 y MEDHTA_T0.