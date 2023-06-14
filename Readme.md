# Guia de como executar as funções

## Configurações iniciais
O primeiro passo é configurar algumas informações, mudando os valores atribuídos para o número de iterações (K) e (P) e a precisão (precision=4) dos números flutuantes impressos no vetor/matriz, caso desejado:

``` 
K = 10      #Número de iterações
P = 5       #Número de iterações adicionais

np.set_printoptions(precision=4, floatmode = 'fixed')  #Configura precisão de 4 dígitos decimais na impressão dos vetores e matrizes
``` 
## Determine a ordem do sistema

Na função main, a variável n dita a ordem do sistema, altere-o para o valor desejado. 
``` 
n = 10
``` 

## As funções

Determinada a ordem do sistema, é possivel que a matriz dos coeficientes seja tridiagonal ou pentadiagonal. Para o primeiro caso, as funções `algoritmoThomas` e `sorTri` solucionam o sistema. Para matriz pentadiagonal, a solução é obtida através das funções `resolve` e `sorPenta`.

## Executando as funções

Para executar as funções, basta descomentá-las na main e preencher os parâmetros com os valores desejados.

**algoritmoThomas(n)** *Algoritmo de Thomas*
>parâmetros: 
    - n - Ordem do sistema 

**algoritmoTri(n, tol = TOL1, w = 1.1)**    *SOR para matriz Tridiagonal*
>parâmetros: 
- n - Ordem do sistema 
- tol - Tolerância usada como critério de parada *(opcional)*
    Padrão: tol = 0
- w - Fator de relaxação *(opcional)*
    Padrão: w = 0

**resolve(n, q = 5)**    *Solução obtida pela função linalg.solve*
>parâmetros: 
- n - Ordem do sistema      
- q - Distância entra as diagonais principal e flutuante *(opcional)*
    Padrão: q = 3                     

**algoritmoPenta(n, q = 5, tol = TOL1, w = 1.1)** *SOR para matriz Pentadiagonal*
>parâmetros
>    n - Ordem do sistema 
>    q - Distância entra as diagonais principal e flutuante *(opcional)*
        Padrão: q = 3
    tol - Tolerância usada como critério de parada *(opcional)*
        Padrão: tol = 0
    w - Fator de relaxação *(opcional)*
        Padrão: w = 0