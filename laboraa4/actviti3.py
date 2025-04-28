from typing import List, Tuple
from Levenshtein import distance as levenshtein_distance
from collections import defaultdict

def boyer_moore_preprocessing(pattern: str) -> Tuple[dict, list]:
    m = len(pattern)
    bad_char = defaultdict(lambda: -1)
    for i, char in enumerate(pattern):
        bad_char[char] = i

    # Good suffix rule preprocessing (standard implementation)
    good_suffix = [0] * (m + 1)
    suff = [0] * (m + 1)
    suff[m] = m
    g = m
    f = 0
    for i in range(m - 1, -1, -1):
        if i > 0 and pattern[i] == pattern[m - 1]:
            suff[i] = suff[i + 1] + 1
        else:
            suff[i] = 0
    for i in range(m + 1):
        good_suffix[i] = m
    j = 0
    for i in range(m - 1, -1, -1):
        if suff[i] == i + 1:
            for j in range(m - i - 1):
                if good_suffix[j] == m:
                    good_suffix[j] = m - i - 1
    for i in range(m - 1):
        good_suffix[m - suff[i]] = m - i - 1

    return bad_char, good_suffix


def boyer_moore_search(text: str, pattern: str) -> List[int]:
    n, m = len(text), len(pattern)
    if m == 0 or n == 0 or m > n:
        return []

    bad_char, good_suffix = boyer_moore_preprocessing(pattern)
    positions = []
    s = 0

    while s <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[s + j]:
            j -= 1

        if j < 0:
            positions.append(s)
            s += good_suffix[0]
        else:
            bc_shift = max(1, j - bad_char[text[s + j]])
            gs_shift = good_suffix[j + 1]
            s += max(bc_shift, gs_shift)

    return positions

def hamming_similarity(s1: str, s2: str) -> float:
    if len(s1) != len(s2):
        return 0.0
    if len(s1) == 0:
        return 1.0
    mismatches = sum(c1 != c2 for c1, c2 in zip(s1, s2))
    return 1.0 - (mismatches / len(s1))

def levenshtein_similarity(s1: str, s2: str) -> float:
    max_len = max(len(s1), len(s2))
    if max_len == 0:
        return 1.0
    dist = levenshtein_distance(s1, s2)
    return 1.0 - (dist / max_len)

def ngram_similarity(s1: str, s2: str, n: int = 2) -> float:
    if n <= 0:
        raise ValueError("n debe ser positivo para n-gramas")

    def get_ngrams(s: str, n_gram: int) -> set:
        if len(s) < n_gram:
            return set()
        return set(s[i:i+n_gram] for i in range(len(s) - n_gram + 1))

    set1 = get_ngrams(s1, n)
    set2 = get_ngrams(s2, n)

    if not set1 and not set2:
        return 1.0
    if not set1 or not set2:
        return 0.0

    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))

    return intersection / union if union != 0 else 0.0


def find_mutants(C: str, P: str, S: float, method: str = 'levenshtein', n_gram_size: int = 2) -> List[List]:
    len_p = len(P)
    len_c = len(C)

    if len_p == 0 or len_c == 0 or len_p > len_c:
        return []

    exact_match_positions = boyer_moore_search(C, P)
    results = []
    processed_indices = set()

    for pos in exact_match_positions:
        results.append([pos, P, 1.0])
        processed_indices.add(pos)

    for i in range(len_c - len_p + 1):
        if i in processed_indices:
            continue

        kmer = C[i : i + len_p]

        sim = 0.0
        if method == 'hamming':
            sim = hamming_similarity(P, kmer)
        elif method == 'levenshtein':
            sim = levenshtein_similarity(P, kmer)
        elif method == 'ngram':
            sim = ngram_similarity(P, kmer, n=n_gram_size)
        else:
            raise ValueError("Método de similitud no válido. Usar 'hamming', 'levenshtein' o 'ngram'.")

        if sim >= S:
            results.append([i, kmer, round(sim, 4)])
            processed_indices.add(i)

    results.sort(key=lambda x: x[0])

    return results

def main():
    C_ejemplo = "CGCCCGAATCCAGAACGCATTCCCCTGGCCTCCATTCTGGAA CGGTACGGACGTCAATCAAAT".replace(" ", "")
    P_ejemplo = "ATTCTGGA"
    S_ejemplo = 0.625

    print(f"Secuencia C: {C_ejemplo}")
    print(f"Patrón P: {P_ejemplo}")
    print(f"Umbral S: {S_ejemplo}")
    print("-" * 30)

    print("Resultados usando Hamming:")
    mutantes_hamming = find_mutants(C_ejemplo, P_ejemplo, S_ejemplo, method='hamming')
    print(mutantes_hamming)
    print("-" * 30)

    print("Resultados usando Levenshtein:")
    mutantes_levenshtein = find_mutants(C_ejemplo, P_ejemplo, S_ejemplo, method='levenshtein')
    print(mutantes_levenshtein)
    print("-" * 30)

    print("Resultados usando N-gramas (n=2):")
    mutantes_ngram2 = find_mutants(C_ejemplo, P_ejemplo, S_ejemplo, method='ngram', n_gram_size=2)
    print(mutantes_ngram2)
    print("-" * 30)

    print("Resultados usando N-gramas (n=3):")
    mutantes_ngram3 = find_mutants(C_ejemplo, P_ejemplo, S_ejemplo, method='ngram', n_gram_size=3)
    print(mutantes_ngram3)
    print("-" * 30)

if __name__ == "__main__":
    main()
