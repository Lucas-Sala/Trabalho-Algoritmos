# Guia para execução as funções

## Configurações iniciais
O primeiro passo é configurar algumas informações, mudando os valores atribuídos para o número de iterações (K) e (P) e a precisão (precision = 4) dos números flutuantes impressos no vetor/matriz, caso desejado:

```Python
K = 10      #Número de iterações
P = 5       #Número de iterações adicionais

np.set_printoptions(precision = 4, floatmode = 'fixed')  #Configura precisão de 4 dígitos na impressão dos vetores e matrizes
``` 
## As funções

É possivel que a matriz dos coeficientes seja tridiagonal ou pentadiagonal. Para o primeiro caso, as funções `algoritmoThomas` e `sorTri` solucionam o sistema. Para matriz pentadiagonal, a solução é obtida através das funções `resolve` e `sorPenta`.

**Algoritmo de Thomas**

`algoritmoThomas(n)` 
```Python
parâmetros: 
    n   Ordem do sistema
```

**SOR para matriz Tridiagonal**

`algoritmoTri(n, tol, w)`
```Python
parâmetros: 
    n     Ordem do sistema
    tol   Tolerância usada como critério de parada
    w     Fator de relaxação 'opcional'
```

**Solução obtida pela função linalg.solve**

`resolve(n, q)`
```Python
parâmetros: 
    n   Ordem do sistema      
    q   Distância entre a diagonal principal e diagonal flutuante
```            

**SOR para matriz Pentadiagonal**

`algoritmoPenta(n, q, tol, w)`
```Python
parâmetros:
    n     Ordem do sistema 
    q     Distância entre a diagonal principal e diagonal flutuante
    tol   Tolerância usada como critério de parada
    w     Fator de relaxação 'opcional'
```

## Executando as funções

Para executar as funções, basta descomentá-las (retirar o caractere #) na main e preencher os parâmetros com os valores desejados, obedecendo algumas restrições.

>tol = TOL1 ou tol = TOL2

>0 < w < 1 (sub-relaxação) ou 1 < w < 2 (sobre-relaxação)

>3 < q <= n

As funções imprimem o tempo gasto na execução do método (em milissegundos) e retornam o vetor solução (nos casos em que sua determinação é possível). Para exibí-lo no terminal, basta incluir a chamada da função dentro da função print. Ex: `print(algoritmoTri(7,TOL1,1.2))`

## Exemplos de chamadas (execuções)

```Python
#-------------------------------- Matriz Tridiagonal --------------------------------#
    algoritmoThomas(10)                             #Algoritmo de Thomas

    algoritmoTri(n = 10, tol = TOL1)                #SOR
    algoritmoTri(n = 10, tol = TOL1, w = 1.1)    
    print(algoritmoTri(n = 10, tol = TOL2)) 

#------------------------------- Matriz Pentadiagonal -------------------------------#
    resolve(n = 10, q = 5)                                 #linalg.solve
              
    algoritmoPenta(n = 10, q = 5, tol = TOL1)              #SOR
    algoritmoPenta(n = 10, q = 5, tol = TOL1, w = 1.1)
    print(algoritmoPenta(n = 10, q = 5, tol = TOL2))
```

Observe que é possível fazer chamada das funções `algoritmoTri` e `algoritmoPenta` sem passar o parâmetro *w*. Caso o parâmetro *w* não for especificado, um valor ótimo para *w* será calculado.