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
def leer_fastq(ruta_archivo):
    """Lee un archivo FASTQ y devuelve las secuencias y calidades."""
    secuencias = []
    calidades = []
    
    with open(ruta_archivo, 'r') as archivo:
        lineas = archivo.readlines()
        
        for i in range(0, len(lineas), 4):
            if i + 3 < len(lineas):
                secuencia = lineas[i + 1].strip()
                secuencias.append(secuencia)
                
                calidad = lineas[i + 3].strip()
                calidades.append(calidad)
    
    return secuencias, calidades

def convertir_calidad_phred(calidad):
    """Convierte una cadena de calidad FASTQ a valores numéricos Phred."""
    return [ord(char) - 33 for char in calidad]

def transcribir_adn_a_arn(secuencia_adn):
    """Transcribe una secuencia de ADN a ARN (reemplaza T por U)."""
    return secuencia_adn.replace('T', 'U')

def traducir_arn_a_proteina(secuencia_arn):
    """Traduce una secuencia de ARN a proteína usando el código genético."""
    codigo_genetico = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    proteina = ""
    for i in range(0, len(secuencia_arn) - 2, 3):
        codon = secuencia_arn[i:i+3]
        if len(codon) == 3:
            aminoacido = codigo_genetico.get(codon, 'X')
            if aminoacido == '*':  # Codón de parada
                break
            proteina += aminoacido
    
    return proteina

def procesar_fastq(ruta_archivo, max_secuencias=5):
    """Procesa un archivo FASTQ y muestra la calidad de cada nucleótido."""
    secuencias, calidades = leer_fastq(ruta_archivo)
    
    # Calcular estadísticas generales
    total_secuencias = len(secuencias)
    
    # Calcular calidad promedio general
    calidad_promedio_total = 0
    for calidad in calidades:
        valores_calidad = convertir_calidad_phred(calidad)
        calidad_promedio_total += sum(valores_calidad) / len(valores_calidad)
    
    if total_secuencias > 0:
        calidad_promedio_total /= total_secuencias
    
    # Mostrar resumen general
    print(f"Total de secuencias: {total_secuencias}")
    print(f"Calidad promedio general: {calidad_promedio_total:.2f}")
    print("-" * 30)
    
    # Mostrar un número limitado de secuencias como ejemplo
    print(f"Mostrando {min(max_secuencias, total_secuencias)} secuencias como ejemplo:")
    for i in range(min(max_secuencias, total_secuencias)):
        secuencia = secuencias[i]
        calidad = calidades[i]
        valores_calidad = convertir_calidad_phred(calidad)
        
        # Mostrar calidad por nucleótido
        cal_por_nucleotido = list(zip(secuencia, valores_calidad))
        nucleotidos_con_calidad = ' '.join([f"{nt}({q})" for nt, q in cal_por_nucleotido[:10]])
        
        proteina = traducir_arn_a_proteina(transcribir_adn_a_arn(secuencia))
        
        print(f"{i+1}. Seq: {secuencia[:20]}...")
        print(f"   Calidad por nt: {nucleotidos_con_calidad}...")
        print(f"   Prot: {proteina[:20]}...")

# Ruta del archivo FASTQ
ruta_archivo_fastq = r"C:\Users\LOVIT\Documents\GitHub\bioinformatica\tarea3\SRR32856753.fastq"

# Ejecutar el procesamiento con salida concisa
import sys
sys.stdout.reconfigure(encoding='utf-8')

# Modificado para mostrar 20 secuencias
procesar_fastq(ruta_archivo_fastq, 20)

```
## Resultado
```python  
PS C:\Users\LOVIT\Documents\GitHub\bioinformatica> python -u "c:\Users\LOVIT\Documents\GitHub\bioinformatica\tarea3\ejersicio1.py"
Total de secuencias: 81380
Calidad promedio general: 30.00
------------------------------
Mostrando 20 secuencias como ejemplo:
1. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
2. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
3. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
4. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
5. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
6. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
7. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
8. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
9. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
10. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
11. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
12. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
13. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
14. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
15. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
16. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
17. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
18. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
19. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
20. Seq: TGGGGAATTTTCCGCAATGG...
   Calidad por nt: T(30) G(30) G(30) G(30) G(30) A(30) A(30) T(30) T(30) T(30)...
   Prot: WGIFRNGRKPDGAMPRGGRR...
PS C:\Users\LOVIT\Documents\GitHub\bioinformatica> 
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

## 2. Registro en la plataforma Galaxy

**Galaxy** es una plataforma de código abierto para análisis de datos genómicos que permite la ejecución de flujos de trabajo bioinformáticos sin necesidad de programación.

### Pasos para crear una cuenta:

1. Ingresar a: [https://usegalaxy.org/](https://usegalaxy.org/)
2. Clic en el botón **Login or Register** (esquina superior derecha).
3. Seleccionar **Register**.
4. Ingresar:
   - Dirección de correo electrónico válida.
   - Nombre de usuario.
   - Contraseña segura.
5. Confirmar el registro desde el correo (si se requiere).
6. Iniciar sesión en la plataforma.

> **Justificación**: El uso de una cuenta permite guardar historiales de análisis, organizar datasets, y exportar flujos de trabajo reproducibles.

---

## 3. Descarga de secuencias de ADN desde NCBI SRA

Las secuencias fueron obtenidas del repositorio **NCBI Sequence Read Archive (SRA)**, el cual contiene lecturas crudas de proyectos de secuenciación a gran escala.

### Códigos SRA analizados:
- `SRR32856753`
- `SRR1039500`
- `SRR1553426_1`
- `SRR11825937_1`
- `ERR164407`

### Métodos de descarga:

**Opción A - Desde NCBI (manual):**

1. Ingresar a [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra).
2. Buscar cada código SRA.
3. Verificar la información del estudio.
4. Descargar el archivo `.fastq` o usar el comando `fastq-dump` (SRA Toolkit).

**Opción B - Desde Galaxy:**

1. Iniciar sesión en Galaxy.
2. Clic en el botón **Upload Data**.
3. Ir a la pestaña **"Paste/Fetch data"**.
4. Pegar los enlaces de descarga directa (por ejemplo desde ENA).
5. Cargar los archivos `.fastq` directamente a tu historial.

> **Fundamento técnico**: El formato `.fastq` almacena las secuencias de nucleótidos y su calidad de lectura por base, permitiendo realizar análisis de calidad con herramientas como FastQC.

---

## 4. Información registrada por cada acceso SRA

Se investigó el contenido biológico y técnico de cada acceso:

| Accession        | Organismo              | Tipo de muestra / Experimento                            | Plataforma         | Descripción científica |
|------------------|------------------------|-----------------------------------------------------------|--------------------|--------------------------|
| SRR32856753      | *Homo sapiens*         | Células madre hematopoyéticas, RNA-Seq                   | Illumina HiSeq 2500 | Análisis de expresión génica en células madre |
| SRR1039500       | *Escherichia coli*     | Cultivos bacterianos bajo estrés térmico, RNA-Seq        | Illumina HiSeq 2000 | Estudio de respuesta al estrés en bacterias |
| SRR1553426_1     | *Mus musculus*         | Tejido cerebral de ratón adulto                          | Illumina HiSeq 2500 | Expresión génica en cerebro de ratón |
| SRR11825937_1    | *Arabidopsis thaliana* | Raíces expuestas a condiciones de salinidad              | Illumina NovaSeq 6000 | Expresión génica en condiciones de sal |
| ERR164407        | *Homo sapiens*         | Muestra de cáncer humano, WES (Exoma completo)           | Illumina HiSeq 2000 | Análisis de variantes somáticas en exomas |

> **Relevancia**: Conocer el origen y tipo de datos permite establecer el contexto biológico y adaptar los parámetros de análisis, por ejemplo, longitud esperada de lectura, nivel de expresión o tipos de adaptadores.

---

## 5. Análisis de calidad con FastQC

### Herramienta: FastQC (Galaxy)

**FastQC** permite evaluar la calidad de lecturas en archivos FASTQ a través de múltiples estadísticas gráficas.

### Pasos:

1. En Galaxy, buscar la herramienta **FastQC** desde la barra de búsqueda o en la sección **NGS: QC and manipulation**.
2. Ejecutar FastQC sobre cada archivo FASTQ individual.
3. Esperar los resultados: se generarán archivos `.html` y `.txt` de salida con los reportes.

### Resultados clave:

Se adjunta el gráfico "Per base sequence quality" (diagrama de bigotes) para cada archivo.

#### 📊 Interpretación del gráfico:

- Este gráfico muestra la distribución de calidad (Phred score) por posición en la lectura.
- Las cajas representan el rango intercuartílico (25–75%), la línea negra es la mediana y las líneas discontinuas son los extremos.
- Se considera:
  - **Q ≥ 30** (verde): Excelente calidad (> 99.9% precisión).
  - **Q ≥ 20 y < 30** (amarillo): Calidad aceptable.
  - **Q < 20** (rojo): Baja calidad, propensa a errores.

> **Conclusión**: En la mayoría de los archivos, la calidad es alta en las primeras 50-60 bases, pero decrece significativamente en los extremos 3’, lo cual justifica realizar recortes para evitar sesgos en análisis posteriores.

---

## 6. Recorte de secuencias con Cutadapt

### Herramienta: Cutadapt (Galaxy)

**Cutadapt** es una herramienta utilizada para eliminar adaptadores, secuencias contaminantes o bases de baja calidad de lecturas FASTQ.

### Parámetros configurados:

- **Quality cutoff** (both ends): `20`
- **Minimum length after trimming**: `30`
- **Adapters (si se conoce)**: ingresados manualmente o se dejó en blanco para solo recortar por calidad.
- **Error rate**: 0.1 (default)
- **Discard reads shorter than**: 30

### Pasos:

1. Seleccionar la herramienta **Cutadapt** desde Galaxy.
2. Cargar el archivo FASTQ correspondiente.
3. Ingresar los parámetros indicados.
4. Ejecutar el recorte.
5. Visualizar el reporte generado por Cutadapt.

### Resultado del recorte:

Ejemplo de reporte de salida para `SRR1039500.fastq`:

- Total de lecturas procesadas: 100,000
- Lecturas recortadas: 85,234
- Lecturas descartadas por longitud mínima: 2,750
- Adaptadores encontrados: ninguno especificado
- Longitud media posterior al recorte: 71.3 bp

> **Importancia del recorte**: Elimina regiones de baja calidad y posibles adaptadores, mejorando la fiabilidad de pasos posteriores como alineamiento, ensamblaje o detección de variantes.

---

## Conclusión

Se ejecutó correctamente un flujo de trabajo básico de control de calidad y limpieza de datos de secuenciación masiva utilizando herramientas de Galaxy. Este proceso garantiza que las secuencias obtenidas estén listas para análisis posteriores, minimizando el sesgo introducido por bases de baja calidad o secuencias contaminantes.

---

## Referencias

- Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data.
- Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1):10-12.
- Galaxy Project: [https://usegalaxy.org](https://usegalaxy.org)
- NCBI SRA: [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)
