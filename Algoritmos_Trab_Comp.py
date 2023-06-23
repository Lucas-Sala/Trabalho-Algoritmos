import numpy as np
import math
import warnings
import time

warnings.filterwarnings("error")

TOL1 = 10**-3       #Tolerâncias usadas como critério de parada no método SOR
TOL2 = 10**-7

K = 10      #Número de iterações
P = 5       #Número de iterações adicionais

np.set_printoptions(precision = 4, floatmode = 'fixed')   #Configura precisão de 4 dígitos decimais na impressão dos vetores e matrizes

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
    '''
    Algoritmo que coordena a execução do algoritmo de Thomas. Cria os vetores/matriz do sistema, faz a chamada da função 'thomas',
    imprime no terminal o vetor solução e o tempo de execução da função 'thomas'.
    '''
    matriz_BT, vetor_b, vetor_x  = iniciaSistemaTridiagonal(n)

    tempo_inicial = time.time()
    thomas(matriz_BT, vetor_b, vetor_x, n)      #Chamada do método de Thomas
    tempo_final = time.time()

    print("\n")
    print("Algoritmo de Thomas: ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms") 
    return vetor_x

def sorTri(matriz_BT, vetor_b, matriz_x, n, iteracoes, tolerancia = 0, w = 1):
    '''
    Método iterativo SOR para matriz tridiagonal. Modifica a matriz_x que consiste em dois vetores: 
    matriz_x[0] - vetor x da iteração mais recente (x^(k+1));
    metriz_x[1] - vetor x da última iteração (x^(k)).
    '''
    for j in range(iteracoes):
        matriz_x[1] = np.copy(matriz_x[0]) 
        for i in range(n):
            x_d1 = 0 if i==0 else matriz_x[0,i-1]           #Captação do valor de x correspontente às diagonais
            x_d2 = 0 if i == n-1 else matriz_x[0,i+1]

                #    x[i] = (1 - w) * x[i] + (w / A[i,i]) * (b[i] - Soma)
            matriz_x[0,i] = (1-w)*matriz_x[1,i] + (w/matriz_BT[i,1])*(vetor_b[i] - matriz_BT[i,0]*x_d1 - matriz_BT[i,2]*x_d2)
        
        err_rel = (max(abs(matriz_x[0]-matriz_x[1])))/max(abs(matriz_x[0]))     #Erro relativo
        if err_rel <= tolerancia:
            print("Tolerância atingida na iteração", j+1)
            break

def algoritmoTri(n, tol, w = 0):
    '''
    Algoritmo responsável por coordenar a execução do método SOR para matriz tridiagonal. 
    Imprime no terminal, o vetor solução, o tempo gasto na execução do método SOR e o fator w utilizado.
    '''
    matriz_x = np.zeros((2,n))      #Matriz que armazena o vetor x de duas iterações [x^(k+1) e x^(k)]
    matriz_BT, vetor_b, matriz_x[0] = iniciaSistemaTridiagonal(n)

    if w == 0:      #Se w não foi especificado calcula w ótimo, executando o "algoritmo 01"
        sorTri(matriz_BT,vetor_b, matriz_x, n, K)   
        e_k = max(abs((matriz_x[0]-matriz_x[1])))

        sorTri(matriz_BT,vetor_b, matriz_x, n, P)   
        e_kp = max(abs((matriz_x[0]-matriz_x[1])))

        w = 2/(1+math.sqrt(1-((e_kp/e_k)**(1/P))))
    print("\nFator w =", w)
    try:
        tempo_inicial = time.time()
        sorTri(matriz_BT, vetor_b, matriz_x, n, 100000, tol, w)      #Chamada do método SOR
        tempo_final = time.time()
        
        print("Algoritmo SOR para matriz Tridiagonal: ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms")
        return matriz_x[0]
    except:
        print("Método SOR não convergiu!")      #Caso o método não seja convergente para o w dado, um overflow ocorre e esse trecho é acionado
        return " "

def sorPenta(matriz_BP, vetor_b, matriz_x, n, q, iteracoes, tolerancia = 0, w = 1):
    '''
    Método iterativo SOR para matriz pentadiagonal. Modifica a matriz_x que consiste em dois vetores: 
    matriz_x[0] - vetor x da iteração mais recente (x^(k+1));
    metriz_x[1] - vetor x da última iteração (x^(k)).
    '''
    for j in range(iteracoes):
        matriz_x[1] = np.copy(matriz_x[0]) 
        for i in range(n):
            x_d1 = 0 if i==0 else matriz_x[0,i-1]           #Captação do valor de x correspontente às diagonais
            x_df1 = 0 if i-q+1 < 0 else matriz_x[0,i-q+1]
            x_d2 = 0 if i == n-1 else matriz_x[0,i+1]
            x_df2 = 0 if i+q-1 >= n else matriz_x[0,i+q-1]
                #    x[i] = (1 - w) * x[i] + (w / A[i][i]) * (b[i] - Soma)
            matriz_x[0,i] = (1-w)*matriz_x[1,i] + (w/matriz_BP[i,2])*(vetor_b[i] - matriz_BP[i,0]*x_df1 - matriz_BP[i,1]*x_d1 - matriz_BP[i,3]*x_d2 - matriz_BP[i,4]*x_df2)
       
        err_rel = (max(abs(matriz_x[0]-matriz_x[1])))/max(abs(matriz_x[0]))     #Erro relativo
        if err_rel <= tolerancia:
            print("Tolerância atingida na iteração", j+1)
            break

def algoritmoPenta(n, q, tol, w = 0):
    '''
    Algoritmo responsável por coordenar a execução do método SOR para matriz pentadiagonal. 
    Imprime no terminal, o vetor solução, o tempo gasto na execução do método SOR e o fator w utilizado.
    '''
    matriz_x = np.zeros((2,n))      #Matriz que armazena o vetor x de duas iterações [x^(k+1) e x^(k)]
    matriz_BP, vetor_b, matriz_x[0] = iniciaSistemaPentadiagonal(n, q)
    if w == 0:      #Se w não foi especificado calcula w ótimo, executando o "algoritmo 01"
        sorPenta(matriz_BP,vetor_b, matriz_x, n, q, K)
        e_k = max(abs((matriz_x[0]-matriz_x[1])))

        sorPenta(matriz_BP,vetor_b, matriz_x, n, q, P)
        e_kp = max(abs((matriz_x[0]-matriz_x[1])))

        w = 2/(1+math.sqrt(1-((e_kp/e_k)**(1/P))))
    print("\nFator w =", w)
    try:
        tempo_inicial = time.time()
        sorPenta(matriz_BP, vetor_b, matriz_x, n, q, 100000, tol, w)     #Chamada do método SOR 
        tempo_final = time.time()

        print("Algoritmo SOR para matriz Pentadiagonal: ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms") 
        return matriz_x[0]
    except:
        print("Método SOR não convergiu!")      #Caso o método não seja convergente para o w dado, um overflow ocorre e esse trecho é acionado
        return " "

def resolve(n, q):
    '''
    Inicia a matriz A e os vetores b e x do sistema pentadiagonal.
    Imprime no terminal o vetor solução e o tempo gasto na execução da função linalg.solve.
    '''
    matriz_A = matrizA(n, q)
    vetor_b = matriz_A.sum(axis = 1)        #Iniciação do sistema
    vetor_x = np.zeros((n))

    tempo_inicial = time.time()
    vetor_x = np.linalg.solve(matriz_A,vetor_b)     #Chamada do método 'solve'
    tempo_final = time.time()

    print("\n")
    print("Algoritmo de resolução da biblioteca 'linalg': ","%.4f" % ((tempo_final - tempo_inicial)*10**3), "ms")
    return vetor_x

def main():
    #   n  -  Ordem do sistema
    #   q  -  Distância entre as diagonais principal e flutuante
    # tol  -  Tolerância usada como critério de parada
    #   w  -  Fator de relaxação (Opcional)

##-------------------------------- Matriz Tridiagonal --------------------------------##
    # algoritmoThomas(n)                    #Algoritmo de Thomas
    # algoritmoTri(n, tol, w)               #SOR 
    
##------------------------------- Matriz Pentadiagonal -------------------------------##
    # resolve(n, q)                         #linalg.solve
    # algoritmoPenta(n, q, tol, w)          #SOR

    algoritmoThomas(200)

    algoritmoTri(500, TOL1, 1.3)

    resolve(100, 15)

    print(algoritmoPenta(10, 5, TOL2))

main()
