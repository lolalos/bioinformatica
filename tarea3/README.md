# Universidad Nacional de San Antonio Abad del Cusco

## Asignatura: Bioinformática
## Práctica: Archivos FASTQ

## Profesora: María del Pilar Venegas Vergara 
## Alumno: Efraín Vitorino Marín Cod: 160337

## Propósito:
El propósito de esta práctica es proporcionar una comprensión profunda de los datos de secuencias de nucleótidos en formato FASTQ, un formato ampliamente utilizado en bioinformática para almacenar datos de secuenciación. Además, se busca familiarizar al estudiante con herramientas bioinformáticas como FastQC y Galaxy, que son esenciales para la evaluación y el análisis de la calidad de los datos secuenciados. Finalmente, se pretende que el estudiante implemente una función que permita determinar la calidad de los nucleótidos en archivos FASTQ, reforzando así sus habilidades prácticas y teóricas.

## Trabajo
### Actividad 1: Responder a las siguientes interrogantes

1. **¿Cuáles son las secciones de un formato FASTA?**  
    El formato FASTA consta de dos secciones principales:  
    - Una línea de encabezado que comienza con el carácter `>` seguida de una descripción opcional.  
    - Las secuencias de nucleótidos o aminoácidos representadas en líneas subsecuentes.  
    Este formato es utilizado principalmente para almacenar secuencias biológicas de manera simple y legible.

2. **¿Cuáles son las secciones de un formato FASTQ?**  
    El formato FASTQ tiene cuatro secciones por cada entrada:  
    - Una línea de encabezado que comienza con `@` seguida de un identificador único.  
    - La secuencia de nucleótidos.  
    - Una línea separadora que comienza con `+`, opcionalmente seguida de información adicional.  
    - Una línea de valores de calidad en formato ASCII que corresponde a cada nucleótido de la secuencia.  
    Este formato combina información de secuencia y calidad, lo que lo hace ideal para datos de secuenciación masiva.

3. **¿Qué información guarda un archivo FASTQ?**  
    Un archivo FASTQ almacena:  
    - Las secuencias de nucleótidos obtenidas de la secuenciación.  
    - La calidad asociada a cada nucleótido, representada en un formato codificado (Phred score).  
    Esto permite evaluar la precisión de las lecturas y realizar análisis posteriores con datos confiables.

4. **¿Qué es Galaxy? Menciona algunas aplicaciones que tiene Galaxy.**  
    Galaxy es una plataforma bioinformática basada en la web que permite realizar análisis de datos biológicos sin necesidad de conocimientos avanzados de programación. Algunas de sus aplicaciones incluyen:  
    - Análisis de datos de secuenciación masiva (NGS).  
    - Procesamiento de datos genómicos, transcriptómicos y proteómicos.  
    - Visualización de datos biológicos.  
    - Automatización de flujos de trabajo bioinformáticos.  
    Galaxy es ampliamente utilizado por su accesibilidad y capacidad para integrar múltiples herramientas en un solo entorno.

5. **¿Qué puede realizar la aplicación FastQC?**  
    FastQC es una herramienta diseñada para evaluar la calidad de los datos de secuenciación. Sus funciones incluyen:  
    - Generar informes detallados sobre la calidad de las lecturas.  
    - Identificar problemas comunes como sesgos en la secuenciación o baja calidad en ciertas regiones.  
    - Proporcionar gráficos y estadísticas para interpretar los datos.  
    Es una herramienta esencial para garantizar que los datos sean adecuados para análisis posteriores.

6. **¿Qué puede realizar la aplicación Cutadapt?**  
    Cutadapt es una herramienta utilizada para:  
    - Recortar adaptadores y secuencias no deseadas de las lecturas.  
    - Eliminar regiones de baja calidad en los extremos de las secuencias.  
    - Preparar los datos para análisis bioinformáticos posteriores.  
    Su uso es crucial para limpiar los datos de secuenciación y mejorar la precisión de los resultados.
### Actividad 2: Implementar los procesos de identificación de secuencias, transcripción de ADN y traducción de ARN en Python

En esta actividad, se busca trabajar con datos biológicos en formato FASTQ, que es ampliamente utilizado para almacenar secuencias de ADN junto con información de calidad. Además, se utilizarán herramientas bioinformáticas para analizar y procesar estas secuencias.

#### Código en Python para procesar un archivo FASTQ:
El siguiente código utiliza la biblioteca `Biopython`, una herramienta poderosa para trabajar con datos biológicos. Este script permite leer un archivo FASTQ, extraer las secuencias y sus calidades, y mostrarlas en la consola.

```python
from Bio import SeqIO

def calcular_calidad_fastq(archivo_fastq):
    """
    Función para leer un archivo FASTQ y mostrar la información de las secuencias y su calidad.
    
    Parámetros:
    archivo_fastq (str): Ruta al archivo FASTQ que se desea procesar.
    """
    with open(archivo_fastq, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            print(f"ID: {record.id}")
            print(f"Secuencia: {record.seq}")
            print(f"Calidad: {record.letter_annotations['phred_quality']}")

# Ruta del archivo FASTQ
archivo_fastq = "ruta/a/tu/archivo.fastq"
calcular_calidad_fastq(archivo_fastq)
```

**Explicación del código:**
1. **Importación de módulos:** Se importa `SeqIO` de `Biopython`, que permite manejar formatos de secuencias como FASTQ.
2. **Función `calcular_calidad_fastq`:** 
   - Abre el archivo FASTQ en modo lectura.
   - Itera sobre cada registro en el archivo utilizando `SeqIO.parse`.
   - Extrae y muestra el identificador (`record.id`), la secuencia (`record.seq`) y la calidad de las bases en formato Phred (`record.letter_annotations['phred_quality']`).
3. **Uso del script:** Se define la ruta del archivo FASTQ y se llama a la función para procesarlo.

#### Pasos adicionales:
Para complementar el análisis de las secuencias, se realizarán las siguientes actividades:

1. **Crear una cuenta en Galaxy:**  
   Galaxy es una plataforma bioinformática en línea que permite realizar análisis de datos genómicos sin necesidad de instalar software localmente. Puedes registrarte en:  
   [https://usegalaxy.org/](https://usegalaxy.org/)

2. **Descargar secuencias de ADN:**  
   Se trabajará con las siguientes muestras en formato FASTQ, disponibles en la base de datos NCBI SRA:  
   - SRR32856753  
   - SRR1039500  
   - SRR1553426_1  
   - SRR11825937_1  
   - ERR164407  
   Puedes descargarlas desde:  
   [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)

3. **Análisis de las secuencias:**  
   Responder la pregunta:  
   **¿Qué información registran las secuencias SRR32856753, SRR1039500, SRR1553426_1, SRR11825937_1, ERR164407?**  
   *(Dejar espacio para completar manualmente en el informe).*

4. **Evaluación de calidad con FastQC:**  
   - Utilizar la herramienta **FastQC** en Galaxy para analizar las muestras FASTQ descargadas.
   - Adjuntar el diagrama de bigotes “Per base sequence quality” en el informe. Este gráfico muestra la calidad de las bases a lo largo de las secuencias, lo que ayuda a identificar regiones de baja calidad.  
     ![Diagrama de bigotes](ruta/a/imagen.png)

5. **Recorte de secuencias con Cutadapt:**  
   - Usar la herramienta **Cutadapt** en Galaxy para recortar las regiones de baja calidad en las secuencias.
   - Capturar y adjuntar el resultado del recorte en el informe. Esto asegura que las secuencias procesadas sean de alta calidad para análisis posteriores.  
     ![Resultado del recorte](ruta/a/imagen_recorte.png)

#### Argumentación:
El uso de herramientas como Python y Galaxy permite automatizar y simplificar el análisis de datos genómicos, que de otra manera sería tedioso y propenso a errores. Python, con bibliotecas como `Biopython`, es ideal para manejar grandes volúmenes de datos biológicos, mientras que Galaxy proporciona una interfaz gráfica amigable para realizar análisis avanzados sin necesidad de conocimientos profundos de programación.

Además, el análisis de calidad de las secuencias es crucial para garantizar la fiabilidad de los resultados en estudios genómicos. Herramientas como FastQC y Cutadapt son estándares en la comunidad científica para este propósito.

Este flujo de trabajo asegura que las secuencias procesadas sean de alta calidad y listas para aplicaciones posteriores, como la identificación de genes, la transcripción de ADN a ARN y la traducción de ARN a proteínas.

#### Código en Python para procesar un archivo FASTQ:
```python
def leer_fastq(ruta_archivo):
    secuencias = []
    calidades = []
    with open(ruta_archivo, 'r', encoding='utf-8') as archivo:
        while True:
            encabezado = archivo.readline().strip()
            if not encabezado:
                break
            secuencia = archivo.readline().strip()
            archivo.readline()  # Línea de separador
            calidad = archivo.readline().strip()
            secuencias.append(secuencia)
            calidades.append(calidad)
    return secuencias, calidades

def convertir_calidad_phred(calidad):
    return [ord(c) - 33 for c in calidad]

def procesar_fastq(ruta_archivo):
    secuencias, calidades = leer_fastq(ruta_archivo)
    ruta_salida = ruta_archivo.replace('.fastq', '_resultado.txt')

    with open(ruta_salida, 'w', encoding='utf-8') as archivo_salida:
        for i, (secuencia, calidad) in enumerate(zip(secuencias, calidades)):
            valores_calidad = convertir_calidad_phred(calidad)
            resultado = (
                f"🔬 Secuencia {i+1}: {secuencia}\n"
                f"📊 Calidad numérica: {valores_calidad}\n"
                f"📉 Calidad promedio: {sum(valores_calidad) / len(valores_calidad):.2f}\n"
                + '-' * 50 + '\n'
            )
            print(resultado)  # Mostrar en consola
            archivo_salida.write(resultado)  # Escribir en el archivo

# Ruta del archivo FASTQ
ruta_archivo_fastq = r"C:\Users\LOVIT\OneDrive\Documentos\GitHub\bioinformatica\Tarea 3\SRR32856753.fastq"

# Ejecutar el procesamiento
import sys
sys.stdout.reconfigure(encoding='utf-8')

procesar_fastq(ruta_archivo_fastq)

```
## resultado 
```python
from Bio import SeqIO

def calcular_calidad_fastq(archivo_fastq):
    with open(archivo_fastq, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            print(f"ID: {record.id}")
            print(f"Secuencia: {record.seq}")
            print(f"Calidad: {record.letter_annotations['phred_quality']}")

# Ruta del archivo FASTQ
archivo_fastq = "ruta/a/tu/archivo.fastq"
calcular_calidad_fastq(archivo_fastq)
```


#### Pasos adicionales:
1. Crear una cuenta en la plataforma Galaxy: [https://usegalaxy.org/](https://usegalaxy.org/)

2. Descargar secuencias de ADN:  
   - SRR32856753  
   - SRR1039500  
   - SRR1553426_1  
   - SRR11825937_1  
   - ERR164407  
   en formato FASTQ desde NCBI: [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)

3. **¿Qué información registran las secuencias SRR32856753, SRR1039500, SRR1553426_1, SRR11825937_1, ERR164407?**  
   (Dejar espacio para completar manualmente).

4. Utilizar FastQC de Galaxy para procesar las muestras FASTQ descargadas. Adjuntar el diagrama de bigotes “Per base sequence quality” en el informe y realizar su interpretación.  
   ![Diagrama de bigotes](ruta/a/imagen.png)

5. Luego, utilizar la aplicación Cutadapt de Galaxy para recortar las secuencias que tienen poca calidad. Capturar el resultado del recorte y adjuntar en el informe.  
   ![Resultado del recorte](ruta/a/imagen_recorte.png)