# UNIVERSIDAD NACIONAL DE SAN ANTONIO ABAD DEL CUSCO

## ASIGNATURA: BIOINFORMÁTICA
## PRÁCTICA: TRANSCRIPCIÓN, TRADUCCIÓN
## profesora : MARIA DEL PILAR VENEGAS VERGARA
## Alumno: Efrain Vitorino Marin

#  análisis detallado de los puntos clave de esta práctica:

1.  **Ejecución y Fundamentación de Actividades:**

    *   **Acciones Solicitadas:** Se requiere la ejecución precisa de cada paso en las actividades de transcripción y traducción. Esto implica aplicar herramientas bioinformáticas y algoritmos específicos para simular y analizar estos procesos biológicos fundamentales.
    *   **Justificación Bioinformática:** Cada acción debe estar intrínsecamente ligada a su relevancia bioinformática. Por ejemplo, al simular la transcripción, es crucial comprender cómo las secuencias promotoras influyen en la eficiencia de la transcripción y cómo las diferentes isoformas de ARNm pueden surgir del splicing alternativo. En la traducción, se debe analizar el impacto de las modificaciones postraduccionales en la estructura y función de las proteínas.
    *   **Interrogantes:** Las respuestas a las interrogantes deben demostrar una comprensión profunda de los mecanismos moleculares subyacentes y cómo las herramientas bioinformáticas nos permiten modelarlos y analizarlos. Se espera que se utilicen bases de datos biológicas (NCBI, UniProt) para validar y enriquecer las respuestas.

2.  **Elaboración de Informe Técnico-Científico:**

    *   **Metodología:** El informe debe detallar la metodología utilizada, incluyendo las herramientas de software (e.g., EMBOSS, Biopython), las bases de datos consultadas y los parámetros de configuración empleados. Se debe justificar la elección de cada herramienta y parámetro en función de los objetivos de la práctica.
    *   **Resultados:** Los resultados deben presentarse de manera clara y concisa, utilizando tablas, gráficos y visualizaciones generadas a partir de los datos obtenidos. Se deben incluir métricas de evaluación relevantes, como la precisión de la predicción de sitios de inicio de la transcripción o la exactitud de la traducción in silico.
    *   **Análisis:** El análisis debe interpretar los resultados a la luz de la teoría biológica y los hallazgos previos en la literatura científica. Se deben discutir las limitaciones de los métodos utilizados y proponer posibles mejoras o experimentos complementarios. Es fundamental demostrar una comprensión crítica de los resultados y su significado biológico.
    *   **Estructura:** El informe debe seguir una estructura lógica y coherente, incluyendo una introducción, una sección de materiales y métodos, una sección de resultados, una sección de discusión y una conclusión. Se deben citar las fuentes bibliográficas utilizadas de acuerdo con un estilo de citación estándar (e.g., APA, MLA).

3.  **Entrega a través del Aula Virtual:**

    *   **Administración Centralizada:** La entrega a través del aula virtual facilita la gestión y el seguimiento de los trabajos. Permite una evaluación eficiente y una retroalimentación oportuna.
    *   **Estándares Académicos:** El proceso de entrega y revisión se adhiere a los estándares académicos establecidos por la universidad. Esto garantiza la calidad y la integridad de la evaluación.
    *   **Retroalimentación:** La retroalimentación proporcionada a través del aula virtual tiene como objetivo mejorar la comprensión y el desempeño del estudiante. Se espera que los estudiantes utilicen esta retroalimentación para reflexionar sobre su aprendizaje y mejorar sus habilidades bioinformáticas.
## Ejercicio 1: Transcripción de ADN a ARN

1. **Desarrollo de un aplicativo en Python para la transcripción de ADN a ARN:**

    * **Objetivo:** Crear un programa en Python que simule el proceso biológico de la transcripción, convirtiendo una secuencia de ADN (en formato FASTA) a una secuencia de ARN (en formato FASTA).
    
    * **Fundamento biológico:** La transcripción es un proceso celular esencial donde la información genética codificada en el ADN se copia a una molécula de ARN. Durante este proceso, la enzima ARN polimerasa sintetiza una cadena de ARN complementaria utilizando una hebra de ADN como plantilla, reemplazando las timinas (T) por uracilos (U).
    
    * **Requisitos técnicos:**
        - Lectura de archivos en formato FASTA (cabecera iniciando con ">" seguida de la secuencia)
        - Validación de la secuencia de entrada (verificar que solo contenga nucleótidos válidos)
        - Implementación del algoritmo de transcripción (sustitución de T por U)
        - Generación del archivo de salida en formato FASTA con la secuencia de ARN resultante
        
    * **Implementación sugerida:**
        - Uso de la biblioteca Biopython para el manejo de archivos FASTA
        - Definición de funciones para la lectura, validación, transcripción y escritura
        - Manejo de excepciones para casos de secuencias inválidas o problemas con los archivos
        
    * **Resultados esperados:**
        - Archivo FASTA con la secuencia de ARN transcrita correctamente
        - Informe sobre el proceso (longitud de la secuencia, cambios realizados)

 ## python
```python
def transcribe_dna_to_rna(input_fasta):
    """
    Transcribe a DNA sequence to RNA from a FASTA file and prints the result.

    Args:
        input_fasta (str): Path to the input FASTA file containing DNA sequences.
    """
    with open(input_fasta, 'r') as infile:
        for line in infile:
            if line.startswith('>'):  # Header line
                print(line.strip())  # Print header
            else:  # DNA sequence line
                rna_sequence = line.strip().replace('T', 'U')
                print(rna_sequence)  # Print transcribed sequence

if __name__ == "__main__":
    # Input files
    files = [
        r"C:\Users\LOVIT\OneDrive\Documentos\GitHub\bioinformatica\tarea2\transcripcion_traduccion_01.fasta",
        r"C:\Users\LOVIT\OneDrive\Documentos\GitHub\bioinformatica\tarea2\transcripcion_traduccion_02.fasta"
    ]

    for input_fasta in files:
        try:
            print(f"\nProcessing file: {input_fasta}")
            # Call the function to perform transcription and display the result
            transcribe_dna_to_rna(input_fasta)
        except FileNotFoundError:
            print(f"Error: Input file not found: {input_fasta}. Please check the path.")
        except Exception as e:
            print(f"Unexpected error while processing {input_fasta}: {e}")
```

## Resultado 
```python
>Macaco_japones
UCAAGACCUAGCAGUGGGACAGCCAGACAGACGGCACGAUGGCACUGAGCUCCCAGAUCUGGGCCACUUGCCUUCUCCUCCUUCUCCUCCUCGCCAGCCUGACCAGUGGCUCCGUUUUCCCACAACAGACGGGACAACUUGCAGAGCUGCAACCUCAGGACAGAGCUGGAGCCAGGGACAGCUGGACGCCCAUGCUCCAGAGGCGAAGGAGGCGAGACACCCACUUCCCCAUCUGCAUUUUCUGCUGCGGCUGCUGUCAUCGAUCAAAGUGUGGGAUGUGCUGCAGGACGUAGAACCUUCCUGCCCUGCCCCCAUCCCCUCCCUUCCUUAUUUAUUCCUGCUGCCCCAGAACACAGGUCUUGGAAUAAAACGGCUGAUUCUUUUGUUUUCC
>Macaco_cangrejero
UCAAGACCUAGCAGUGGGACAGCCAGACAGACGGCACGAUGGCACUGAGCUCCCAGAUCUGGGCCACUUGCCUCCUCCUCCUUCUCCUCCUCGCCAGCCUGACCAGUGGCUCCGUUUUCCCACAACAGACGGGACAACUUGCAGAGCUGCAACCUCAGGACAGAGCUGGAGCCAGGGCCAGCUGGACGCCCAUGCUCCAGAGGCGAAGGAGGCGAGACACCCACUUCCCCAUCUGCAUUUUCUGCUGCGGCUGCUGUCAUCGAUCAAAGUGUGGGAUGUGCUGCAGGACGUAGAACCUUCCUGCCCUGCCCCCAUCCCCUCCCUUCCUUAUUUAUUCCUGCUGCCCCAGAACACAGGUCUUGGAAUAAAACGGCUGAUUCUUUUGUUUUCC
>Gorila
UCAAGACCCAGCAGUGGGACAGCCAGACAGACGGCACGAUGGCACUGAGCUCCCAGAUCUGGGCCGCUUGCCUCCUGCUCCUCCUCCUCCUCGCCAGCCUGACCAGUGGCUCUGUUUUCCCACAACAGACGGGACAACUUGCAGAGCUGCAACCCCAGGACAGAGCUGGAGCCAGGGCCGGCUGGACGCCCAUGCUCCAGAGGCGAAGGAGGCGAGACACCCACUUCCCCAUCUGCAUUUUCUGCUGCGGCUGCUGUCAUCGAUCAAAGUGUGGGAUGUGCUGCAAGACGUAGAACCUACCUGCCCUGCCCCCGUCCCCUCCCUUCCUUAUUUAUUCCUGCUGCCCCAGAACAUAGGUCUUGGAAUAAAAUGGCUGGUUCUUUUGUUUUCC
>Titi_comun
UCAAGACCCAGCAGUGGGACAGCCGGACAGACAGCACGAUGGCACUGAGCUCCCAGAUCUGGGCCGUCUGCCUCUUCCUCCUCCUCCUCCUCACCAGCCUGACCAGUGGCUUUGUUUUCCCACAACAGGCAGGACAACUCACCGAGCUUCAACCCCAGGACAGAGCUGGAGCCAGGGCCAGCUGGAUGCCCAUGAUCCAGAGGCGAAGAAGGCGAGACACCCACUUCCCCAUCUGCAUUUUCUGCUGCGGCUGCUGUCGUCAAUCAAACUGUGGGAUGUGCUGCAAGACGUAGAACCUUCCUGCCCUGCCCCCGUCUCCUCCCUUCCUUAUUUAUUCCUGCUGCCCCAGAACACAGGUCUGGGAAUAAAAUGGCUGGUUCUUUUGUUUUCC
>Tota
UCAAGACCUAGCAGUGGGACAGCCAGACAGACGGCACGAUGGCACUGAGCUCCCAGAUCUGGGCCACUUGCCUCUUCCUCCUUCUCCUCCUCGCCAGCCUGACCAGUGGCUCCGUUUUCCCACAACAGACGGGACAACUUGCAGAGCUGCAACCUCAGGACAGAGCUGGAGCCAGGGCCAGCUGGACGCCCAUGCUCCAGAGGCGAAGGAGGCGAGACACCCACUUCCCCAUCUGCAUUUUCUGCUGUGGCUGCUGUCAUCGAUCAAAGUGUGGGAUGUGCUGCAGGACGUAGAAUCUUCCUGCCCUGCCCCCAUCUCCUCCCUUCCUUAUUUAUUCCUGCUGCCCCAGAAUACAGGUCUUGGAAUAAAACGGCUGAUUCUUUUGUUUUCC
>Mono_arana
UCAAGACCCAGCAGUGGGACAGCCGGACAGACAGCACGAUGGCACUGAGCUCCCAGAUCUGGGCCGUCUGCUUCUUCCUCCUCCUCCUCCUCACCAGCCUGACCAGUGGCUUUGUUUUCCCACAACAGAGAGGACAACUCACUGAGCUUCAACCCCAGGACAGAGCUGGAGCCAGGGCCGGCUGGACGCCCAUGAUCCAGAGGCGAAGGAGGCGAGACACCCACUUCCCCAUCUGCAUUUUCUGCUGCGGCUGCUGUCGUCAACCAAACUGUGGGAUGUGCUGCAAGACGUAGAACCUUCCCGCCCCGCCCCCAUCCCUCCCUUCCUUAUUUAUUCCUGCUGCCCCAGAACACAGGUCUGGGAAUAAAAUGGCUGGUUCUUUUGUUUUCC
>Gibon_cresta
UCAAGACCCAGCAGUGGGACAGCCAGACAGACGGCACGAUGGCACUGAGCUCCCAGAUCUGGGCCGCUUGCCUCCUCCUCCUCCUCCUCCUCGCCAGCCUGACCAGUGGCUCCGUUUUCCCACAACAGACGGGACAACUUGCAGAGCUGCAAGGACAGCUGGACAGAGCUGGAGCCAGGGCCGGCUGGACGCCCAUGCUCCAGAGGCGAAGGAGGCGAGACACCCACUUCCCCAUCUGCAUUUUCUGCUGCGGCUGCUGUCAUCGAUCAAAGUGUGGGAUGUGCUGCAAGACGUAGAACCUUCCUGCCCUGCCCCCGUCCCCUCCCUUCCUUAUUUAUUCCUGCUGCCCCAGAACACAGGUCUUGGAAUAAAAUGGCUGGUUCUUUUGUUUUCC
>Gibon_manos
UCAAGACCCAGCAGUGGGACAGCCAGACAGACGGCACGAUGGCACUGAGCUCCCAGAUCUGGGCUGCUUGCCUCCUCCUCCUCCUCCUCCUCGCCAGCCUGACCAGUGGCUCCGUUUUCCCACAACAGAUGGGACAACUUGCAGAGCUGCAAGGACAGCUGGACAGAGCUGGAGCCAGGGCCGGCUGGACGCCCAUGCUCCAGAGGCGAAGGAGGCGAGACACCCACUUCCCCAUCUGCAUUUUCUGCUGCGGCUGCUGUCAUCGAUCAAAGUGUGGGAUGUGCUGCAAGACGUAGAACCUUCCUGCCCUGCCCCCGUCCCCUCCCUUCCUUAUUUAUUCCUGCUGCCCCAGAACACAGGUCUUGGAAUAAAAUGGCUGGUUCUUUUGUUUUCC
>Hombre
UCAAGACCCAGCAGUGGGACAGCCAGACAGACGGCACGAUGGCACUGAGCUCCCAGAUCUGGGCCGCUUGCCUCCUGCUCCUCCUCCUCCUCGCCAGCCUGACCAGUGGCUCUGUUUUCCCACAACAGACGGGACAACUUGCAGAGCUGCAACCCCAGGACAGAGCUGGAGCCAGGGCCAGCUGGAUGCCCAUGUUCCAGAGGCGAAGGAGGCGAGACACCCACUUCCCCAUCUGCAUUUUCUGCUGCGGCUGCUGUCAUCGAUCAAAGUGUGGGAUGUGCUGCAAGACGUAGAACCUACCUGCCCUGCCCCCGUCCCCUCCCUUCCUUAUUUAUUCCUGCUGCCCCAGAACAUAGGUCUUGGAAUAAAAUGGCUGGUUCUUUUGUUUUCC
```

## Ejercicio 2: Traducción de ARN a proteína

2. **Desarrollo de un aplicativo en Python para la traducción de ARN a proteína:**

    * **Objetivo:** Crear un programa en Python que simule el proceso biológico de la traducción, convirtiendo una secuencia de ARN (en formato FASTA) a una secuencia de aminoácidos (en formato FASTA).
    
    * **Fundamento biológico:** La traducción es el proceso mediante el cual la información genética contenida en el ARN mensajero (ARNm) se utiliza para sintetizar proteínas. Los ribosomas leen el ARNm en grupos de tres nucleótidos (codones), y cada codón especifica un aminoácido particular según el código genético universal, con algunos codones funcionando como señales de inicio (AUG) y de terminación (UAA, UAG, UGA).
    
    * **Requisitos técnicos:**
        - Lectura de archivos en formato FASTA con secuencias de ARN
        - Validación de la secuencia de entrada (verificar que solo contenga nucleótidos válidos de ARN)
        - Implementación del algoritmo de traducción (conversión de codones a aminoácidos)
        - Detección adecuada de codones de inicio y terminación
        - Generación del archivo de salida en formato FASTA con la secuencia de aminoácidos resultante
        
    * **Implementación sugerida:**
        - Uso de la biblioteca Biopython para el manejo de archivos FASTA y traducción
        - Definición de un diccionario de código genético para mapear codones a aminoácidos
        - Manejo de marcos de lectura y posibles múltiples productos de traducción
        - Tratamiento adecuado de codones de inicio y terminación
        
    * **Resultados esperados:**
        - Archivo FASTA con la secuencia de aminoácidos correctamente traducida
        - Informe sobre el proceso (longitud de la secuencia, posibles ORFs identificados)
 ## python
 ```python
def read_fasta(input_file):
    sequences = []
    header = None
    sequence = ""
    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequences.append((header, sequence))
                header = line  # nueva cabecera encontrada
                sequence = ""
            else:
                sequence += line.upper().replace(" ", "")
        if header is not None:
            sequences.append((header, sequence))
    return sequences

def translate_rna(rna_seq):
    # Tabla de codones para ARN estándar.
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    protein = ""
    # Recorrer la secuencia en grupos de 3 nucleótidos
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        amino_acid = codon_table.get(codon, 'X')  # X para codones desconocidos
        protein += amino_acid
        if amino_acid == '*':  # detener la traducción en codón de terminación
            break
    return protein

def write_fasta(header, sequence, output_file):
    with open(output_file, "a") as f:
        f.write(header + "\n")
        f.write(sequence + "\n")

def main():
    input_path = r"C:\Users\LOVIT\OneDrive\Documentos\GitHub\bioinformatica\tarea2\transcripcion_traduccion_01.fasta"
    output_path = r"C:\Users\LOVIT\OneDrive\Documentos\GitHub\bioinformatica\tarea2\transcripcion_traduccion_02.fasta"
    
    sequences = read_fasta(input_path)
    all_proteins = ""
    print("Secuencias traducidas:")
    for header, rna_seq in sequences:
        protein_seq = translate_rna(rna_seq)
        result_header = header if header.startswith(">") else ">Traducción"
        all_proteins += result_header + "\n"
        for i in range(0, len(protein_seq), 60):
            all_proteins += protein_seq[i:i+60] + "\n"
        print(result_header)
        print(protein_seq)
    
    with open(output_path, "w") as f:
        f.write(all_proteins)
    result_header = header if header.startswith(">") else ">Traducción"
    write_fasta(result_header, protein_seq, output_path)
    
    print("Secuencia traducida:")
    print(result_header)
    print(protein_seq)

if __name__ == "__main__":
    main()
```
## Resultado 
 ```python
Secuencias traducidas:
>exon18_MYH16_chimpance
EQXNKXXTTXXSTAPXXXRXXXPXEXKQX
>exon18_MYH16_humano
EQXNKXXTTXXSRTPXXPXXXPQXXXAX
Secuencia traducida:
>exon18_MYH16_humano
EQXNKXXTTXXSRTPXXPXXXPQXXXAX
```
