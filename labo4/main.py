from Bio import SeqIO
import Levenshtein

def detectar_exones(gen_seq, exones_ref):
    exones_detectados = []
    last_pos = 0
    gen_seq = gen_seq.upper()
    
    for exon in exones_ref:
        exon = exon.upper()
        best_match = ''
        best_score = 0
        best_pos = 0

        for i in range(last_pos, len(gen_seq) - len(exon) + 1):
            subseq = gen_seq[i:i+len(exon)]
            score = Levenshtein.ratio(subseq, exon)

            if score > best_score:
                best_score = score
                best_match = subseq
                best_pos = i

            if score == 1.0:
                break

        if best_score > 0.8:
            exones_detectados.append(best_match)
            last_pos = best_pos + len(exon)
        else:
            exones_detectados.append('')
    return exones_detectados

def leer_exones_referencia(path):
    exones_por_gen = {}
    with open(path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
        if len(lines) % 2 != 0:
            raise ValueError("El archivo de exones de referencia no tiene líneas en pares (gen y exones).")
        for i in range(0, len(lines), 2):
            gen_id = lines[i]
            exones = [exon.strip() for exon in lines[i + 1].split(",")]
            exones_por_gen[gen_id] = exones
    return exones_por_gen

def procesar_fasta(archivo_fasta, archivo_exones_ref, archivo_salida):
    exones_por_gen = leer_exones_referencia(archivo_exones_ref)

    with open(archivo_salida, "w") as out_fasta:
        for gen in SeqIO.parse(archivo_fasta, "fasta"):
            gen_id = gen.id.split()[0]
            A = str(gen.seq).upper()
            R = exones_por_gen.get(gen_id, [])
            S_parts = detectar_exones(A, R)
            S = ''.join(S_parts)

            out_fasta.write(f">{gen.id}\n{S}\n")

            print(f"Gen: {gen.id}")
            print(f"A: {A}")
            print(f"R: {', '.join(R)}")
            print(f"S: {S}")
            print("-" * 50)

# Archivos
archivo_fasta = "ei.fasta"
archivo_exones_ref = "ei_exon_ref.fasta"
archivo_salida = "exones_empalmados.fasta"

procesar_fasta(archivo_fasta, archivo_exones_ref, archivo_salida)
print("✅ Proceso completado. Resultados guardados en exones_empalmados.fasta.")
