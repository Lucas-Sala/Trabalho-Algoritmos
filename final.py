import numpy as np
import math
import time

TOL1 = 10**-3 
TOL2 = 10**-7

np.set_printoptions(precision=4, floatmode = 'fixed')

k = 10      #Número de iterações
p = 5       #Número de operações adicionais

def matrizA(n,q):
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

def matrizBT(n):
    matriz_BT = np.zeros((n,3))
    for i in range(n):
        if i != 0: 
            matriz_BT[i,0] = -1 
            matriz_BT[i-1,2] = -2
        matriz_BT[i,1] = 8
    return matriz_BT

def matrizBP(n, q):
    matriz_BP = np.zeros((n,5))
    for i in range(n):
        if i != 0: 
            matriz_BP[i,1] = -1 
            matriz_BP[i-1,3] = -2
        if i <= n-q: 
            matriz_BP[i,4] = -2
            matriz_BP[i+q-1,0] = -1
        matriz_BP[i,2] = 8
    return matriz_BP

def thomas(matriz_BT, vetor_b, vetor_x, n):
    a = np.copy(matriz_BT.T[0])
    b = np.copy(matriz_BT.T[1])
    c = np.copy(matriz_BT.T[2])
    d = np.copy(vetor_b)

    c[0] = c[0]/b[0]
    d[0] = d[0]/b[0]
    for i in range(1, n):
        c[i] = c[i]/(b[i]-(a[i]*c[i-1]))
        d[i] = (d[i]-(a[i]*d[i-1]))/(b[i]-(a[i]*c[i-1]))
    
    vetor_x[n-1] = d[n-1]
    for i in range(n-2,-1,-1):
        vetor_x[i] = d[i] - (c[i]*vetor_x[i+1])

def algoritmoThomas(n):
    matriz_BT = matrizBT(n)
    vetor_b = matriz_BT.sum(axis = 1)
    vetor_x = np.zeros(n)

    tempo_inicial = time.time()
    thomas(matriz_BT, vetor_b, vetor_x, n)
    tempo_final = time.time()

    print(vetor_x)
    print("Algoritmo de Thomas: ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms") 

def sorTri(matriz_BT, vetor_b, matriz_x, w, n, iteracoes):
    for j in range(iteracoes):
        matriz_x[1] = np.copy(matriz_x[0]) 
        for i in range(n):
            x_d1 = 0 if i==0 else matriz_x[0,i-1]
            x_d2 = 0 if i == n-1 else matriz_x[0,i+1]

            matriz_x[0,i] = (w/matriz_BT[i,1])*(vetor_b[i] - matriz_BT[i,0]*x_d1 - matriz_BT[i,2]*x_d2) + (1-w)*matriz_x[1,i]

def algoritmoTri(n):
    matriz_BT = matrizBT(n)              #definindo matriz_BT
    vetor_b = matriz_BT.sum(axis = 1)     #definindo vetor_b
    matriz_x = np.zeros((2,n))             #definindo vetor_x
    matriz_x[0] = np.copy(vetor_b)

    sorTri(matriz_BT,vetor_b, matriz_x, 1, n, k)   #k iterações com w = 1
    e_k = max(abs((matriz_x[0]-matriz_x[1])))

    sorTri(matriz_BT,vetor_b, matriz_x, 1, n, p)   #p iterações com w = 1
    e_kp = max(abs((matriz_x[0]-matriz_x[1])))

    w = 2/(1+math.sqrt(1-((e_kp/e_k)**(1/p))))

    tempo_inicial = time.time()
    sorTri(matriz_BT, vetor_b, matriz_x, w, n, k)    #SOR
    tempo_final = time.time()
    
    print(matriz_x[0])
    print("Algoritmo SOR para matriz Tridiagonal: ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms")
    print("w = ", w)

def wTri(n, w):
    matriz_BT = matrizBT(n)              #definindo matriz_BT
    vetor_b = matriz_BT.sum(axis = 1)     #definindo vetor_b
    matriz_x = np.zeros((2,n))             #definindo vetor_x
    matriz_x[0] = np.copy(vetor_b)

    tempo_inicial = time.time()
    sorTri(matriz_BT, vetor_b, matriz_x, w, n, k)    #SOR
    tempo_final = time.time()
    
    print(matriz_x[0])
    print("Algoritmo SOR para matriz Tridiagonal: ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms")
    print("w = ", w)

def sorPenta(matriz_BP, vetor_b, matriz_x, w, n, q, iteracoes):
    for j in range(iteracoes):
        matriz_x[1] = np.copy(matriz_x[0]) 
        for i in range(n):
            x_d1 = 0 if i==0 else matriz_x[0,i-1]
            x_df1 = 0 if i-q+1 < 0 else matriz_x[0,i-q+1]
            x_d2 = 0 if i == n-1 else matriz_x[0,i+1]
            x_df2 = 0 if i+q-1 >= n else matriz_x[0,i+q-1]

            matriz_x[0,i] = (w/matriz_BP[i,2])*(vetor_b[i] - matriz_BP[i,0]*x_df1 - matriz_BP[i,1]*x_d1 - matriz_BP[i,3]*x_d2 - matriz_BP[i,4]*x_df2) + (1-w)*matriz_x[1,i]

def algoritmoPenta(n, q):
    matriz_BP = matrizBP(n, q)
    vetor_b = matriz_BP.sum(axis = 1)
    matriz_x = np.zeros((2,n))
    matriz_x[0] = np.copy(vetor_b)
    
    sorPenta(matriz_BP,vetor_b, matriz_x, 1, n, q, k)   #k iterações com w = 1
    e_k = max(abs((matriz_x[0]-matriz_x[1])))

    sorPenta(matriz_BP,vetor_b, matriz_x, 1, n, q, p)   #p iterações com w = 1
    e_kp = max(abs((matriz_x[0]-matriz_x[1])))

    w = 2/(1+math.sqrt(1-((e_kp/e_k)**(1/p))))

    tempo_inicial = time.time()
    sorPenta(matriz_BP, vetor_b, matriz_x, w, n, q, k) 
    tempo_final = time.time()

    print(matriz_x[0])
    print("Algoritmo SOR para matriz Pentadiagonal: ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms") 
    print("w = ", w)

def wPenta(n, q, w):
    matriz_BP = matrizBP(n, q)
    vetor_b = matriz_BP.sum(axis = 1)
    matriz_x = np.zeros((2,n))
    matriz_x[0] = np.copy(vetor_b)
    
    tempo_inicial = time.time()
    sorPenta(matriz_BP, vetor_b, matriz_x, w, n, q, k) 
    tempo_final = time.time()

    print(matriz_x[0])
    print("Algoritmo SOR para matriz Pentadiagonal: ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms") 
    print("w = ", w)

def resolve(n, q):
    matriz_A = matrizA(n, q)
    vetor_b = matriz_A.sum(axis = 1)
    vetor_x = np.zeros((n))

    tempo_inicial = time.time()
    vetor_x = np.linalg.solve(matriz_A,vetor_b)
    tempo_final = time.time()

    print(vetor_x)
    print("Algoritmo de resolução da biblioteca 'linalg': ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms") 

def main():
    n = 10
    
#Tridiagonal
    #algoritmoThomas(n)
    #algoritmoTri(n)

#Pentadiagonal
    #q = 5
    #resolve(n, q)
    #algoritmoPenta(n, q)

#Análise para diferentes w
    #w = 1.1
    #wTri(n, w)

    #q = 5
    #wPenta(n, q, w)

main()