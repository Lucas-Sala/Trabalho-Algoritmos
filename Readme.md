# Guia para execução as funções

## Configurações iniciais
O primeiro passo é configurar algumas informações, mudando os valores atribuídos para o número de iterações (K) e (P) e a precisão (precision=4) dos números flutuantes impressos no vetor/matriz, caso desejado:

```Python
K = 10      #Número de iterações
P = 5       #Número de iterações adicionais

np.set_printoptions(precision=4, floatmode = 'fixed')  #Configura precisão de 4 dígitos na impressão dos vetores e matrizes
``` 
## As funções

É possivel que a matriz dos coeficientes seja tridiagonal ou pentadiagonal. Para o primeiro caso, as funções `algoritmoThomas` e `sorTri` solucionam o sistema. Para matriz pentadiagonal, a solução é obtida através das funções `resolve` e `sorPenta`.

## Executando as funções

Para executar as funções, basta descomentá-las (retirar o caractere #) na main e preencher os parâmetros com os valores desejados, obedecendo algumas restrições.

>tol = TOL1 ou tol = TOL2

>0 < w < 1 (sub-relaxação) ou 1 < w < 2 (sobre-relaxação)

>3 < q <= n


**Algoritmo de Thomas**

`algoritmoThomas(n)` 
```
parâmetros: 
    n - Ordem do sistema
```

**SOR para matriz Tridiagonal**

`algoritmoTri(n, tol, w)`
```
parâmetros: 
    n - Ordem do sistema
    tol - Tolerância usada como critério de parada (opcional)
        Padrão: tol = 0
    w - Fator de relaxação (opcional)
        Padrão: w = 0
```

**Solução obtida pela função linalg.solve**

`resolve(n, q)`
```
parâmetros: 
    n - Ordem do sistema      
    q - Distância entra as diagonais principal e flutuante (opcional)
        Padrão: q = 3
```            

**SOR para matriz Pentadiagonal**

`algoritmoPenta(n, q, tol, w)`
```
parâmetros:
    n - Ordem do sistema 
    q - Distância entra as diagonais principal e flutuante (opcional)
        Padrão: q = 3
    tol - Tolerância usada como critério de parada (opcional)
        Padrão: tol = 0
    w - Fator de relaxação (opcional)
        Padrão: w = 0
```

## Exemplos de chamadas (execuções)

```Python
#-------------------------------- Matriz Tridiagonal --------------------------------#
    algoritmoThomas(10)                             #Algoritmo de Thomas

    algoritmoTri(n = 10)                            #SOR
    algoritmoTri(n = 10, tol = TOL1)
    algoritmoTri(n = 10, w = 1.1)   
    algoritmoTri(n = 10, tol = TOL1, w = 1.1)     

#------------------------------- Matriz Pentadiagonal -------------------------------#
    resolve(n = 10, q = 5)                                 #linalg.solve

    algoritmoPenta(n = 10, q = 5)                           #SOR
    algoritmoPenta(n = 10, q = 5, tol = TOL1)      
    algoritmoPenta(n = 10, q = 5, w = 1.1)
    algoritmoPenta(n = 10, q = 5, tol = TOL1, w = 1.1)
```

Observe que é possível fazer chamada das funções `algoritmoTri` e `algoritmo` sem passar os parâmetros *tol* e *w*. Caso o parâmetro *tol* não seja passado, o método SOR executará K iterações. Se *w* não for especificado, um valor ótimo para *w* será calculado.