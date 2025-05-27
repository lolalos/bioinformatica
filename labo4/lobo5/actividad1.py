from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# Importar la matriz BLOSUM62 de forma robusta
try:
    from Bio.Align import substitution_matrices
    blosum62 = substitution_matrices.load("BLOSUM62")
except Exception:
    print("No se pudo importar la matriz BLOSUM62. Instala Biopython correctamente.")
    raise

# Secuencias de ADN
seq1 = "AACGTTTCCAGTCCAAATAGCTAGGC"
seq2 = "AGTCGAAAT"
seq3 = "GAGCAGCTGAACAAGCTGATGACCACCCTCCATAGCACCGCACCCCATTTTGTCCGCTGTATTATTCCCCAATGAGTTTAAGCAATCGG"
seq4 = "ACCGCTCCATAGCCGCACCCCATTTTGTCCGCTGTTTTA"

# Secuencias de proteína
prot5 = (
    "MTPMRKINPLMKLINHSFIDLPTPSSNISAWWNFGSLLGACLILQITTGLFLAMHYSPDASFAFS"
    "SIAHITRDVNYGWIIRYLHANGASMFFICLFLHIGRGLYYGSFLYSETNWIGIILLLATMATA"
    "FMGYVLPWGQMSFWGATVITNLLSAIPYIGTDLVQWIWGGYSVDSPTLTRFFTFHFILPFIIAA"
    "LATLHLLFLHETGSNNPLGITSHSDKITFHPYYTIKDALGLLLFLSLMTLTLFSPDLLGDPDNYTLANPLA"
)
prot6 = (
    "ISAWWNFGSLLGACMILQITTGLFLAMHYYPDASTAFSSIALATRDVNYGWIIRYLHANGASMF"
    "FICLFLHIGRGLYYGSFLYSETNWIGIIFQMSTATAFMGYVLPWGQMSFWGATVITNLLSAIPY"
    "IGTDLVQWIWGGYSVANPLA"
)

pairs = [
    ("Seq1 vs Seq2 (DNA)", seq1, seq2, "dna"),
    ("Seq3 vs Seq4 (DNA)", seq3, seq4, "dna"),
    ("Prot5 vs Prot6 (Protein)", prot5, prot6, "protein")
]

for title, a, b, kind in pairs:
    print(title)
    if kind == "dna":
        # Needleman-Wunsch (global)
        g_alignments = pairwise2.align.globalms(a, b, 2, -1, -0.5, -0.1)
        print("Global (Needleman-Wunsch):")
        if g_alignments:
            print(format_alignment(*g_alignments[0]))
        else:
            print("No alignment found.")

        # Smith-Waterman (local)
        l_alignments = pairwise2.align.localms(a, b, 2, -1, -0.5, -0.1)
        print("Local (Smith-Waterman):")
        if l_alignments:
            print(format_alignment(*l_alignments[0]))
        else:
            print("No alignment found.")

    else:  # kind == "protein"
        # Protein: BLOSUM62, gap open = -10, gap extend = -0.5
        g_alignments = pairwise2.align.globalds(a, b, blosum62, -10, -0.5)
        print("Global (Needleman-Wunsch):")
        if g_alignments:
            print(format_alignment(*g_alignments[0]))
        else:
            print("No alignment found.")

        l_alignments = pairwise2.align.localds(a, b, blosum62, -10, -0.5)
        print("Local (Smith-Waterman):")
        if l_alignments:
            print(format_alignment(*l_alignments[0]))
        else:
            print("No alignment found.")
    print()