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
