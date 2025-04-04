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

    for i, (secuencia, calidad) in enumerate(zip(secuencias, calidades)):
        valores_calidad = convertir_calidad_phred(calidad)
        resultado = (
            f"🔬 Secuencia {i+1}: {secuencia}\n"
            f"📊 Calidad numérica: {valores_calidad}\n"
            f"📉 Calidad promedio: {sum(valores_calidad) / len(valores_calidad):.2f}\n"
            + '-' * 50 + '\n'
        )
        print(resultado)  # Mostrar en consola

# Ruta del archivo FASTQ
ruta_archivo_fastq = r"C:\Users\LOVIT\OneDrive\Documentos\GitHub\bioinformatica\Tarea 3\SRR32856753.fastq"

# Ejecutar el procesamiento
import sys
sys.stdout.reconfigure(encoding='utf-8')

procesar_fastq(ruta_archivo_fastq)
