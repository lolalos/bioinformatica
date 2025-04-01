def transcribe_dna_to_rna(input_fasta):
    """
    Transcribe a DNA sequence to RNA from a FASTA file and prints the result.

    Args:
        input_fasta (str): Path to the input FASTA file containing DNA sequences.
    """
    with open(input_fasta, 'r') as infile:
        for line in infile:
            if line.startswith('>'):  # Header line
                print(line.strip())  # Imprimir encabezado
            else:  # DNA sequence line
                rna_sequence = line.strip().replace('T', 'U')
                print(rna_sequence)  # Imprimir secuencia transcrita

if __name__ == "__main__":
    # Archivos de entrada
    files = [
        "D:\\bioinformatica\\tarea2\\transcripcion_traduccion_01.fasta",
        "D:\\bioinformatica\\tarea2\\transcripcion_traduccion_02.fasta"
    ]

    for input_fasta in files:
        try:
            print(f"\nProcesando archivo: {input_fasta}")
            # Llamar a la función para realizar la transcripción y mostrar el resultado
            transcribe_dna_to_rna(input_fasta)
        except FileNotFoundError:
            print(f"Error: No se encontró el archivo de entrada: {input_fasta}. Verifique la ruta.")
        except Exception as e:
            print(f"Error inesperado al procesar {input_fasta}: {e}")