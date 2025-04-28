# actividad3.py

import time
import random
import matplotlib.pyplot as plt

# ----------------------
# Algoritmo Fuerza Bruta
# ----------------------
def fuerza_bruta(C, P):
    posiciones = []
    for i in range(len(C) - len(P) + 1):
        if C[i:i+len(P)] == P:
            posiciones.append(i)
    return posiciones

# ----------------------
# Algoritmo Boyer-Moore
# ----------------------
def preprocess_boyer_moore(P):
    bad_char = [-1]*256
    for i in range(len(P)):
        bad_char[ord(P[i])] = i
    return bad_char

def boyer_moore(C, P):
    posiciones = []
    bad_char = preprocess_boyer_moore(P)
    m = len(P)
    n = len(C)
    s = 0
    while(s <= n - m):
        j = m - 1
        while j >= 0 and P[j] == C[s+j]:
            j -= 1
        if j < 0:
            posiciones.append(s)
            s += (m - bad_char[ord(C[s+m])] if s+m < n else 1)
        else:
            s += max(1, j - bad_char[ord(C[s+j])])
    return posiciones

# ----------------------
# Función de búsqueda
# ----------------------
def buscar_patron(C, P, A):
    inicio = time.time()
    if A == 'fuerza_bruta':
        posiciones = fuerza_bruta(C, P)
    elif A == 'boyer_moore':
        posiciones = boyer_moore(C, P)
    else:
        raise ValueError("Algoritmo no válido. Use 'fuerza_bruta' o 'boyer_moore'.")
    fin = time.time()
    frecuencia = len(posiciones)
    tiempo = fin - inicio
    return frecuencia, posiciones, tiempo

# ----------------------
# Función para prueba masiva
# ----------------------
def prueba_tiempos():
    bases = ['A', 'C', 'G', 'T']
    tamanos = [1000, 5000, 10000, 20000, 50000]  # Tamaños de secuencias
    tiempos_fb = []  # Fuerza Bruta
    tiempos_bm = []  # Boyer Moore

    P = 'GCAT'  # Patrón fijo

    print(f"{'Tamaño':<10}{'Tiempo FB':<15}{'Tiempo BM':<15}")
    print('-'*40)

    for tam in tamanos:
        C = ''.join(random.choices(bases, k=tam))

        _, _, tiempo_fb = buscar_patron(C, P, 'fuerza_bruta')
        _, _, tiempo_bm = buscar_patron(C, P, 'boyer_moore')

        tiempos_fb.append(tiempo_fb)
        tiempos_bm.append(tiempo_bm)

        print(f"{tam:<10}{tiempo_fb:<15.6f}{tiempo_bm:<15.6f}")

    return tamanos, tiempos_fb, tiempos_bm

# ----------------------
# Función para graficar
# ----------------------
def graficar(tamanos, tiempos_fb, tiempos_bm):
    plt.figure(figsize=(10,6))
    plt.plot(tamanos, tiempos_fb, marker='o', label='Fuerza Bruta')
    plt.plot(tamanos, tiempos_bm, marker='s', label='Boyer Moore')
    plt.title('Comparación de Tiempos de Ejecución')
    plt.xlabel('Tamaño de Secuencia')
    plt.ylabel('Tiempo (segundos)')
    plt.legend()
    plt.grid(True)
    plt.savefig('grafico_comparacion.png')  # Guarda el gráfico
    plt.show()

# ----------------------
# Programa Principal
# ----------------------
if __name__ == "__main__":
    # PRUEBAS DE TIEMPOS
    tamanos, tiempos_fb, tiempos_bm = prueba_tiempos()

    # GRAFICAR
    graficar(tamanos, tiempos_fb, tiempos_bm)
