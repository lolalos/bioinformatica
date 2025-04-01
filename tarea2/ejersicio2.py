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