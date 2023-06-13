import numpy as np
import math
import time

TOL1 = 10**-3       #Tolerâncias usadas como critério de parada no método SOR
TOL2 = 10**-7

K = 10      #Número de iterações
P = 5       #Número de iterações adicionais

np.set_printoptions(precision=4, floatmode = 'fixed')   #Configura precisão de 4 dígitos decimais na impressão dos vetores e matrizes

def iniciaSistemaTridiagonal(n):
    '''Inicia e retorna a matriz com as diagonais não nulas da matriz A (matriz BT), vetor b e vetor x (iniciado com os mesmos valores do vetor b).'''
    matriz_BT = np.zeros((n,3))             
    for i in range(n):
        if i != 0: 
            matriz_BT[i,0] = -1 
            matriz_BT[i-1,2] = -2
        matriz_BT[i,1] = 8
    vetor_b = matriz_BT.sum(axis = 1)
    vetor_x = np.copy(vetor_b)              
    return matriz_BT, vetor_b, vetor_x

def iniciaSistemaPentadiagonal(n, q):
    '''Inicia e retorna a matriz com as diagonais não nulas da matriz A (matriz BP), vetor b e vetor x (iniciado com os mesmos valores do vetor b).'''
    matriz_BP = np.zeros((n,5))
    for i in range(n):
        if i != 0: 
            matriz_BP[i,1] = -1 
            matriz_BP[i-1,3] = -2
        if i <= n-q: 
            matriz_BP[i,4] = -2
            matriz_BP[i+q-1,0] = -1
        matriz_BP[i,2] = 8
    vetor_b = matriz_BP.sum(axis = 1)
    vetor_x= np.copy(vetor_b)
    return matriz_BP, vetor_b, vetor_x

def matrizA(n,q):
    '''Inicia e retorna a matriz A pentadiagonal dos coeficientes do sistema.'''
    matriz_A = np.zeros((n,n))
    for i in range(n):
        if i != 0: 
            matriz_A[i,i-1] = -1 
            matriz_A[n-i-1,n-i] = -2
        if i <= n-q: 
            matriz_A[i,i+q-1] = -2
            matriz_A[i+q-1,i] = -1
        matriz_A[i,i] = 8
    return matriz_A

def thomas(matriz_BT, vetor_b, vetor_x, n):
    '''
    Algoritmo de Thomas. Preenche o vetor_x com a solução do sistema.
    '''
    a = matriz_BT.T[0]
    b = matriz_BT.T[1]
    c = matriz_BT.T[2]
    d = vetor_b

    c[0] = c[0]/b[0]
    d[0] = d[0]/b[0]
    for i in range(1, n):
        c[i] = c[i]/(b[i]-(a[i]*c[i-1]))
        d[i] = (d[i]-(a[i]*d[i-1]))/(b[i]-(a[i]*c[i-1]))
    
    vetor_x[n-1] = d[n-1]
    for i in range(n-2,-1,-1):
        vetor_x[i] = d[i] - (c[i]*vetor_x[i+1])

def algoritmoThomas(n):
    matriz_BT, vetor_b, vetor_x  = iniciaSistemaTridiagonal(n)

    tempo_inicial = time.time()
    thomas(matriz_BT, vetor_b, vetor_x, n)
    tempo_final = time.time()

    print(vetor_x)
    print("Algoritmo de Thomas: ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms") 

def sorTri(matriz_BT, vetor_b, matriz_x, n, iteracoes, tolerancia = 0, w = 1):
    '''
    Método iterativo SOR para matriz tridiagonal. Modifica a matriz_x que consiste em dois vetores: 
    matriz_x[0] - vetor x da iteração mais recente (x^(k+1));
    metriz_x[1] - vetor x da última iteração (x^(k)).
    '''
    for j in range(iteracoes):
        matriz_x[1] = np.copy(matriz_x[0]) 
        for i in range(n):
            x_d1 = 0 if i==0 else matriz_x[0,i-1]
            x_d2 = 0 if i == n-1 else matriz_x[0,i+1]

            matriz_x[0,i] = (w/matriz_BT[i,1])*(vetor_b[i] - matriz_BT[i,0]*x_d1 - matriz_BT[i,2]*x_d2) + (1-w)*matriz_x[1,i]
        tol = (max(abs(matriz_x[0]-matriz_x[1])))/max(abs(matriz_x[0]))
        if tol <= tolerancia:
            break

def algoritmoTri(n, tol = 0, w = 0):
    '''
    Algoritmo responsável por coordenar a execução do método SOR para matriz tridiagonal. 
    Imprime o vetor solução no terminal e o tempo gasto na execução do método SOR.
    '''
    matriz_x = np.zeros((2,n))             
    matriz_BT, vetor_b, matriz_x[0] = iniciaSistemaTridiagonal(n)

    if w == 0:      #Se w não foi especificado calcula w ótimo, executando o "algoritmo 01"
        sorTri(matriz_BT,vetor_b, matriz_x, n, K)   
        e_k = max(abs((matriz_x[0]-matriz_x[1])))

        sorTri(matriz_BT,vetor_b, matriz_x, n, P)   
        e_kp = max(abs((matriz_x[0]-matriz_x[1])))

        w = 2/(1+math.sqrt(1-((e_kp/e_k)**(1/P))))

    tempo_inicial = time.time()
    sorTri(matriz_BT, vetor_b, matriz_x, n, K, tol, w) 
    tempo_final = time.time()
    
    print(matriz_x[0])
    print("Algoritmo SOR para matriz Tridiagonal: ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms")
    print("w =", w)

def sorPenta(matriz_BP, vetor_b, matriz_x, n, q, iteracoes, tolerancia = 0, w = 1):
    '''
    Método iterativo SOR para matriz pentadiagonal. Modifica a matriz_x que consiste em dois vetores: 
    matriz_x[0] - vetor x da iteração mais recente (x^(k+1));
    metriz_x[1] - vetor x da última iteração (x^(k)).
    '''
    for j in range(iteracoes):
        matriz_x[1] = np.copy(matriz_x[0]) 
        for i in range(n):
            x_d1 = 0 if i==0 else matriz_x[0,i-1]
            x_df1 = 0 if i-q+1 < 0 else matriz_x[0,i-q+1]
            x_d2 = 0 if i == n-1 else matriz_x[0,i+1]
            x_df2 = 0 if i+q-1 >= n else matriz_x[0,i+q-1]

            matriz_x[0,i] = (w/matriz_BP[i,2])*(vetor_b[i] - matriz_BP[i,0]*x_df1 - matriz_BP[i,1]*x_d1 - matriz_BP[i,3]*x_d2 - matriz_BP[i,4]*x_df2) + (1-w)*matriz_x[1,i]
        tol = (max(abs(matriz_x[0]-matriz_x[1])))/max(abs(matriz_x[0]))
        if tol <= tolerancia:
            break

def algoritmoPenta(n, q = 3, tol = 0, w = 0):
    '''
    Algoritmo responsável por coordenar a execução do método SOR para matriz pentadiagonal. 
    Imprime o vetor solução no terminal e o tempo gasto na execução do método SOR.
    '''
    matriz_x = np.zeros((2,n))
    matriz_BP, vetor_b, matriz_x[0] = iniciaSistemaPentadiagonal(n, q)
    
    if w == 0:      #Se w não foi especificado calcula w ótimo, executando o "algoritmo 01"
        sorPenta(matriz_BP,vetor_b, matriz_x, n, q, K)
        e_k = max(abs((matriz_x[0]-matriz_x[1])))

        sorPenta(matriz_BP,vetor_b, matriz_x, n, q, P)
        e_kp = max(abs((matriz_x[0]-matriz_x[1])))

        w = 2/(1+math.sqrt(1-((e_kp/e_k)**(1/P))))

    tempo_inicial = time.time()
    sorPenta(matriz_BP, vetor_b, matriz_x, n, q, K, tol, w) 
    tempo_final = time.time()

    print(matriz_x[0])
    print("Algoritmo SOR para matriz Pentadiagonal: ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms") 
    print("w = ", w)

def resolve(n, q):
    '''
    Inicia a matriz A e os vetores b e x do sistema pentadiagonal.
    Imprime no terminal o vetor solução e o tempo gasto na execução da função linalg.solve.
    '''
    matriz_A = matrizA(n, q)
    vetor_b = matriz_A.sum(axis = 1)
    vetor_x = np.zeros((n))

    tempo_inicial = time.time()
    vetor_x = np.linalg.solve(matriz_A,vetor_b)
    tempo_final = time.time()

    print(vetor_x)
    print("Algoritmo de resolução da biblioteca 'linalg': ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms") 

def main():
    n = 10      #Ordem do sistema

#-------------------------------- Matriz Tridiagonal --------------------------------#

    # algoritmoThomas(n)                            #Algoritmo de Thomas

    # algoritmoTri(n, tol = TOL1, w = 1.1)          #SOR

#------------------------------- Matriz Pentadiagonal -------------------------------#

    # resolve(n, q = 5)                                 #linalg.solve

    # algoritmoPenta(n, q = 5, tol = TOL1, w = 1.1)     #SOR

main()