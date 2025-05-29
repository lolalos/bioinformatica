"""
Implementación completa de un ensamblador de genomas usando:
1. Grafos de De Bruijn
2. Algoritmo para encontrar caminos/ciclos eulerianos
3. Ensamblaje de la secuencia genómica final
"""

from collections import defaultdict, deque

class GenomeAssembler:
    def __init__(self, k):
        """
        Inicializa el ensamblador con el tamaño de k-mers a usar
        
        Args:
            k (int): Tamaño de los k-mers para construir el grafo
        """
        self.k = k  # Tamaño de los k-mers
        self.graph = defaultdict(list)  # Grafo de De Bruijn (diccionario de listas)
        self.in_degree = defaultdict(int)  # Contador de grados de entrada
        self.out_degree = defaultdict(int)  # Contador de grados de salida

    def build_de_bruijn_graph(self, reads):
        """
        Construye el grafo de De Bruijn a partir de las lecturas de ADN
        
        Args:
            reads (list): Lista de cadenas de ADN (lecturas)
        """
        for read in reads:
            # Dividir cada lectura en k-mers y crear aristas en el grafo
            for i in range(len(read) - self.k + 1):
                kmer = read[i:i+self.k]
                left = kmer[:-1]  # (k-1)-mer prefijo
                right = kmer[1:]   # (k-1)-mer sufijo
                
                # Añadir arista al grafo
                self.graph[left].append(right)
                # Actualizar grados
                self.out_degree[left] += 1
                self.in_degree[right] += 1

    def find_eulerian_path_or_cycle(self):
        """
        Encuentra un camino o ciclo euleriano en el grafo usando el algoritmo de Hierholzer
        
        Returns:
            list: Lista de nodos en el camino/ciclo euleriano
        """
        # Paso 1: Encontrar nodos inicial y final (para camino euleriano)
        start_node = None
        end_node = None
        
        # Obtener todos los nodos del grafo
        all_nodes = set(self.graph.keys()).union(
            set(x for neighbors in self.graph.values() for x in neighbors))
        
        # Buscar nodos con grados desbalanceados (para camino euleriano)
        for node in all_nodes:
            if self.out_degree.get(node, 0) - self.in_degree.get(node, 0) == 1:
                start_node = node
            elif self.in_degree.get(node, 0) - self.out_degree.get(node, 0) == 1:
                end_node = node
        
        # Si no hay nodos desbalanceados, es un ciclo euleriano
        if start_node is None:
            start_node = next(iter(self.graph.keys())) if self.graph else None
        
        # Paso 2: Algoritmo de Hierholzer para encontrar el camino/ciclo
        stack = []
        path = []
        current_path = [start_node]
        
        while current_path:
            current_node = current_path[-1]
            if self.graph[current_node]:
                # Hay aristas salientes sin visitar
                next_node = self.graph[current_node].pop()
                current_path.append(next_node)
            else:
                # No hay más aristas salientes, añadir al camino
                path.append(current_path.pop())
        
        # Invertir para obtener el orden correcto
        path = path[::-1]
        return path

    def assemble_genome(self, path):
        """
        Ensambla la secuencia genómica a partir del camino euleriano
        
        Args:
            path (list): Camino euleriano (lista de nodos)
            
        Returns:
            str: Secuencia de ADN ensamblada
        """
        if not path:
            return ""
        
        # El genoma comienza con el primer nodo
        genome = path[0]
        
        # Para cada nodo subsiguiente, añadir el último carácter
        for node in path[1:]:
            genome += node[-1]
            
        return genome

    def visualize_graph(self):
        """
        Muestra una representación visual del grafo de De Bruijn
        """
        print("\nGrafo de De Bruijn:")
        print("Formato: Nodo -> Vecino1, Vecino2, ...")
        for node, neighbors in sorted(self.graph.items()):
            print(f"{node} -> {', '.join(neighbors)}")


def run_assembly_case(L, S, K, case_num):
    """
    Ejecuta un caso completo de ensamblaje y muestra los resultados
    
    Args:
        L (list): Lista de lecturas de ADN
        S (int): Tamaño de solapamiento (no usado directamente en esta implementación)
        K (int): Tamaño de k-mers a usar
        case_num (int): Número de caso para mostrar en los resultados
    """
    print(f"\n{'='*50}")
    print(f"Caso de Prueba {case_num}")
    print(f"Lecturas: {L}")
    print(f"K (tamaño de k-mers): {K}")
    print(f"S (solapamiento): {S}")
    print(f"{'='*50}")
    
    # 1. Construir el grafo de De Bruijn
    assembler = GenomeAssembler(K)
    assembler.build_de_bruijn_graph(L)
    
    # 2. Mostrar el grafo
    assembler.visualize_graph()
    
    # 3. Encontrar camino/ciclo euleriano
    path = assembler.find_eulerian_path_or_cycle()
    
    # Determinar si es camino o ciclo
    if path and path[0] == path[-1]:
        print("\nTipo: CICLO euleriano")
    else:
        print("\nTipo: CAMINO euleriano")
    
    print("\nCamino/Ciclo Euleriano encontrado:")
    print(" -> ".join(path))
    
    # 4. Ensamblar el genoma
    genome = assembler.assemble_genome(path)
    print("\nGenoma ensamblado final:")
    print(genome)
    
    return genome


# =============================================================================
# Ejecución de los casos de prueba solicitados
# =============================================================================

if __name__ == "__main__":
    print("="*70)
    print("APLICACIÓN PARA ENSAMBLAJE DE GENOMAS")
    print("="*70)
    print("Este programa implementa:")
    print("- Construcción de grafos de De Bruijn")
    print("- Búsqueda de caminos/ciclos eulerianos")
    print("- Ensamblaje de secuencias genómicas")
    print("="*70)
    
    # Caso 1
    L1 = ['TAG', 'CTG', 'GCT', 'ACT', 'CTA', 'TGC']
    S1 = 1
    K1 = 3
    genome1 = run_assembly_case(L1, S1, K1, 1)
    
    # Caso 2
    L2 = ['TAC', 'CTA', 'TAG', 'ATA', 'CAT', 'ACA', 'CCT', 'GAC', 'CCC', 'ACC']
    S2 = 1
    K2 = 3
    genome2 = run_assembly_case(L2, S2, K2, 2)
    
    # Caso 3
    L3 = ['AACGTTTC', 'TTCCAGTC', 'GTCCAAAT', 'AATAGCTA', 'CTAGGCAT']
    S3 = 3
    K3 = 4
    genome3 = run_assembly_case(L3, S3, K3, 3)
    
    print("\n" + "="*70)
    print("Resumen de resultados:")
    print(f"Caso 1: Genoma ensamblado = {genome1}")
    print(f"Caso 2: Genoma ensamblado = {genome2}")
    print(f"Caso 3: Genoma ensamblado = {genome3}")
    print("="*70)