from Bio import pairwise2
from Bio.pairwise2 import format_alignment
try:
    from Bio.Align.substitution_matrices import load
    blosum62 = load("BLOSUM62")
except ImportError:
    # Fallback for older Biopython versions
    from Bio.SubsMat.MatrixInfo import blosum62

# GitHub Copilot


# Secuencias de ADN
seq1 = "AACGTTTCCAGTCCAAATAGCTAGGC"
seq2 = "AGTCGAAAT"
seq3 = ("GAGCAGCTGAACAAGCTGATGACCACCCTCCATAGCACCGCACCCCATTTTGTCCGCTGT"
    "ATTATCCCCAATGAGTTTAAGCAATCGG")
seq4 = "ACCGCTCCATAGCCGCACCCCATTTTGTCCGCTGTTTA"

# Secuencias de proteína
prot5 = ("MTPMRKINPLMKLINHSFIDLPTPSNISAWWNFGSLLGACLILQITTGLFLAMHYSPDASTAF"
     "SSIAHITRDVNYGWIIRYLHANGASMFFICLFLHIGRGLYYGSFLYSETWNIGIILLLATMATA"
     "FMGYVLPWGQMSFWGATVITNLLSAIPYIGTDLVQWIWGGYSVDSPTLTRFFTFHFILPFIIAAL"
     "ATLHLLFLHETGSNNPLGITSHSDKITFHPYYTIKDALGLLLFLLSLMTLTLFSPDLLGDPDNYTLANPLA")
prot6 = ("ISAWWNFGSLLGACMILQITTGLFLAMHYYPDASTAFSSIALATRDVNYGWIIRYLHANGASMFF"
     "ICLFLHIGRGLYYGSFLYSETWNIGIIFLQMSTATAFMGYVLPWGQMSFWGATVITNLLSAIPYIG"
     "TDLVQWIWGGYSVANPLA")

pairs = [
    ("Seq1 vs Seq2 (ADN)", seq1, seq2, "dna"),
    ("Seq3 vs Seq4 (ADN)", seq3, seq4, "dna"),
]

for title, a, b, kind in pairs:
    print(title)
    if kind == "dna":
        # Needleman-Wunsch (global) con conteo simple match/mismatch
        g = pairwise2.align.globalxx(a, b)
        print("Global (Needleman–Wunsch):")
        print(format_alignment(*g[0]))
        # Smith-Waterman (local)
        l = pairwise2.align.localxx(a, b)
        print("Local (Smith–Waterman):")
        print(format_alignment(*l[0]))
    else:
        # Proteína: BLOSUM62, gap open=-10, gap extend=-0.5
        g = pairwise2.align.globalds(a, b, blosum62, -10, -0.5)
        print("Global (Needleman–Wunsch):")
        print(format_alignment(*g[0]))
        l = pairwise2.align.localds(a, b, blosum62, -10, -0.5)
        print("Local (Smith–Waterman):")
        print(format_alignment(*l[0]))
    print()