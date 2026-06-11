# CÁLCULOS NA MÃO — DETERMINANTES, POLINÔMIOS, RAÍZES E TABELAS

Cola para raciocínio no papel antes de programar.

Este arquivo foca no que pode aparecer como parte manuscrita ou como preparação para montar o código.

---

## 1. Como identificar o tipo de problema

```text
Pedir determinante                -> usar fórmula 2x2, 3x3 ou eliminação
Pedir polinômio característico     -> montar det(A(λ)) = 0
Pedir raízes de função f(x)=0      -> bisseção, falsa posição, Newton, secante
Pedir raízes de polinômio          -> Bhaskara, fatoração, tabela de coeficientes, Bairstow
Pedir estabilidade por tabela      -> tabela de Routh-Hurwitz
Pedir solução de sistema linear    -> Gauss, LU, Cramer, Thomas, Cholesky, Gauss-Seidel
```

---

## 2. Determinante 2x2

Para:

```text
| a  b |
| c  d |
```

Fórmula:

```text
D = ad - bc
```

Exemplo:

```text
| 3   2 |
| -1  2 |
```

```text
D = 3*2 - 2*(-1)
D = 6 + 2
D = 8
```

---

## 3. Determinante 3x3 por cofatores

Para:

```text
| a11  a12  a13 |
| a21  a22  a23 |
| a31  a32  a33 |
```

Expansão pela primeira linha:

```text
D = a11*(a22*a33 - a23*a32)
  - a12*(a21*a33 - a23*a31)
  + a13*(a21*a32 - a22*a31)
```

Exemplo:

```text
| 1  2  3 |
| 0  1  4 |
| 5  6  0 |
```

```text
D = 1*(1*0 - 4*6) - 2*(0*0 - 4*5) + 3*(0*6 - 1*5)
D = 1*(0 - 24) - 2*(0 - 20) + 3*(0 - 5)
D = -24 + 40 - 15
D = 1
```

---

## 4. Determinante 3x3 pela regra de Sarrus

Só serve para matriz 3x3.

Para:

```text
| a  b  c |
| d  e  f |
| g  h  i |
```

Fórmula:

```text
D = aei + bfg + cdh - ceg - bdi - afh
```

Exemplo:

```text
| 1  2  3 |
| 0  1  4 |
| 5  6  0 |
```

```text
D = 1*1*0 + 2*4*5 + 3*0*6 - 3*1*5 - 2*0*0 - 1*4*6
D = 0 + 40 + 0 - 15 - 0 - 24
D = 1
```

---

## 5. Determinante por eliminação

Use para matrizes maiores.

Ideia:

```text
Transformar A em matriz triangular superior.
O determinante vira o produto da diagonal.
```

Atenção:

```text
Trocar duas linhas muda o sinal do determinante.
Multiplicar uma linha por k multiplica o determinante por k.
Somar múltiplo de uma linha em outra não muda o determinante.
```

Exemplo:

```text
A =
| 2  1  1 |
| 4 -6  0 |
|-2  7  2 |
```

Eliminar abaixo do pivô 2:

```text
L2 <- L2 - 2L1 = [0, -8, -2]
L3 <- L3 + L1  = [0,  8,  3]
```

Matriz:

```text
| 2   1   1 |
| 0  -8  -2 |
| 0   8   3 |
```

Eliminar abaixo do pivô -8:

```text
L3 <- L3 + L2 = [0, 0, 1]
```

Triangular:

```text
| 2   1   1 |
| 0  -8  -2 |
| 0   0   1 |
```

Determinante:

```text
D = 2*(-8)*1 = -16
```

---

## 6. Regra de Cramer

Para sistema:

```text
A x = b
```

A solução é:

```text
x1 = D1/D
x2 = D2/D
x3 = D3/D
```

Onde:

```text
D  = det(A)
D1 = det(A trocando coluna 1 por b)
D2 = det(A trocando coluna 2 por b)
D3 = det(A trocando coluna 3 por b)
```

Use Cramer só para sistema pequeno, principalmente 2x2 ou 3x3.

---

## 7. Polinômio característico de sistema massa-mola-amortecedor

Equação matricial:

```text
M x¨ + C x˙ + K x = 0
```

Assume-se solução modal:

```text
x(t) = φ e^(λt)
```

Derivadas:

```text
x˙(t) = λ φ e^(λt)
x¨(t) = λ² φ e^(λt)
```

Substituindo:

```text
(M λ² + C λ + K) φ = 0
```

Para solução não trivial:

```text
det(M λ² + C λ + K) = 0
```

Esse é o polinômio característico.

Se o sistema tem N graus de liberdade:

```text
grau do polinômio = 2N
```

---

## 8. Exemplo de polinômio característico 2x2

Suponha:

```text
A(λ) =
| 2λ² + 3λ + 6      -λ - 2 |
| -λ - 2             λ² + λ + 2 |
```

Para matriz 2x2:

```text
det(A) = a11*a22 - a12*a21
```

Então:

```text
P(λ) = (2λ² + 3λ + 6)(λ² + λ + 2) - (-λ - 2)(-λ - 2)
```

Primeiro produto:

```text
(2λ² + 3λ + 6)(λ² + λ + 2)

2λ²*λ² = 2λ⁴
2λ²*λ  = 2λ³
2λ²*2  = 4λ²

3λ*λ² = 3λ³
3λ*λ  = 3λ²
3λ*2  = 6λ

6*λ² = 6λ²
6*λ  = 6λ
6*2  = 12
```

Somando:

```text
2λ⁴ + 5λ³ + 13λ² + 12λ + 12
```

Segundo produto:

```text
(-λ - 2)(-λ - 2) = (λ + 2)² = λ² + 4λ + 4
```

Logo:

```text
P(λ) = 2λ⁴ + 5λ³ + 13λ² + 12λ + 12 - (λ² + 4λ + 4)

P(λ) = 2λ⁴ + 5λ³ + 12λ² + 8λ + 8
```

Tabela de coeficientes:

```text
grau:       4   3   2   1   0
coef:       2   5   12  8   8
```

Em forma de lista para Python:

```python
coef = [2, 5, 12, 8, 8]
```

---

## 9. Como montar a tabela de coeficientes do polinômio

Para:

```text
P(x) = a4 x⁴ + a3 x³ + a2 x² + a1 x + a0
```

Tabela:

```text
grau        4    3    2    1    0
coef        a4   a3   a2   a1   a0
```

Exemplo:

```text
P(x) = 2x⁴ + 5x³ + 12x² + 8x + 8
```

Tabela:

```text
grau        4    3    2    1    0
coef        2    5    12   8    8
```

Lista para código:

```python
coef = [2, 5, 12, 8, 8]
```

---

## 10. Como avaliar polinômio na mão

Exemplo:

```text
P(x) = 2x⁴ + 5x³ + 12x² + 8x + 8
```

Para x = 1:

```text
P(1) = 2(1)⁴ + 5(1)³ + 12(1)² + 8(1) + 8
P(1) = 2 + 5 + 12 + 8 + 8
P(1) = 35
```

Para x = -1:

```text
P(-1) = 2(-1)⁴ + 5(-1)³ + 12(-1)² + 8(-1) + 8
P(-1) = 2 - 5 + 12 - 8 + 8
P(-1) = 9
```

---

## 11. Tabela de Horner

Serve para avaliar polinômio mais rápido.

Exemplo:

```text
P(x) = 2x⁴ + 5x³ + 12x² + 8x + 8
```

Coeficientes:

```text
2   5   12   8   8
```

Avaliar em x = 1:

```text
        2    5    12    8    8
x=1          2     7    19   27
        -------------------------
        2    7    19   27   35
```

Resultado:

```text
P(1) = 35
```

Como fazer:

```text
1. Desce o primeiro coeficiente.
2. Multiplica pelo x escolhido.
3. Soma com o próximo coeficiente.
4. Repete até o final.
```

---

## 12. Divisão sintética por Horner

Serve para dividir o polinômio por:

```text
(x - r)
```

Exemplo:

```text
P(x) = x³ - 6x² + 11x - 6
```

Coeficientes:

```text
1   -6   11   -6
```

Testar raiz r = 1.

Tabela:

```text
        1   -6    11   -6
r=1         1    -5     6
        -------------------
        1   -5     6    0
```

Resto zero, então x = 1 é raiz.

Quociente:

```text
x² - 5x + 6
```

Agora fatorando:

```text
x² - 5x + 6 = (x - 2)(x - 3)
```

Raízes:

```text
x = 1, 2, 3
```

---

## 13. Fórmula de Bhaskara

Para:

```text
ax² + bx + c = 0
```

Discriminante:

```text
Δ = b² - 4ac
```

Raízes:

```text
x = (-b ± sqrt(Δ))/(2a)
```

Casos:

```text
Δ > 0  -> duas raízes reais
Δ = 0  -> uma raiz real dupla
Δ < 0  -> duas raízes complexas conjugadas
```

Se Δ < 0:

```text
sqrt(Δ) = i sqrt(|Δ|)
```

Exemplo:

```text
x² - 5x + 6 = 0
```

```text
a = 1
b = -5
c = 6

Δ = (-5)² - 4(1)(6)
Δ = 25 - 24
Δ = 1
```

```text
x = (5 ± 1)/2
x1 = 3
x2 = 2
```

---

## 14. Método da bisseção na mão

Usado para:

```text
f(x) = 0
```

Exige intervalo [a,b] com troca de sinal:

```text
f(a)*f(b) < 0
```

Fórmula do ponto médio:

```text
xr = (a+b)/2
```

Tabela:

```text
iteração | a | b | xr | f(a) | f(xr) | novo intervalo
```

Regra:

```text
Se f(a)*f(xr) < 0, a raiz está entre a e xr -> b = xr
Se f(a)*f(xr) > 0, a raiz está entre xr e b -> a = xr
Se f(xr) = 0, achou a raiz exata
```

Exemplo:

```text
f(x) = x³ - x - 2
```

Intervalo:

```text
a = 1, b = 2
```

```text
f(1) = 1 - 1 - 2 = -2
f(2) = 8 - 2 - 2 = 4
```

Tem troca de sinal.

Iteração 1:

```text
xr = (1+2)/2 = 1.5
f(1.5) = 1.5³ - 1.5 - 2 = 3.375 - 3.5 = -0.125
```

Como:

```text
f(1)*f(1.5) = (-2)(-0.125) > 0
```

Novo intervalo:

```text
[1.5, 2]
```

---

## 15. Método da falsa posição na mão

Também exige:

```text
f(a)*f(b) < 0
```

Fórmula:

```text
xr = b - f(b)*(a-b)/(f(a)-f(b))
```

Forma equivalente:

```text
xr = a - f(a)*(b-a)/(f(b)-f(a))
```

Tabela:

```text
iteração | a | b | f(a) | f(b) | xr | f(xr)
```

Regra de atualização é igual à bisseção:

```text
Se f(a)*f(xr) < 0 -> b = xr
Se f(a)*f(xr) > 0 -> a = xr
```

Exemplo:

```text
f(x) = x³ - x - 2
a = 1, b = 2
f(a) = -2, f(b) = 4
```

```text
xr = 2 - 4*(1-2)/(-2-4)
xr = 2 - 4*(-1)/(-6)
xr = 2 - 4/6
xr = 1.3333
```

```text
f(1.3333) ≈ 1.3333³ - 1.3333 - 2
f(1.3333) ≈ -0.963
```

Como f(a) e f(xr) têm mesmo sinal:

```text
novo intervalo = [1.3333, 2]
```

---

## 16. Método de Newton-Raphson na mão

Fórmula:

```text
x_{i+1} = x_i - f(x_i)/f'(x_i)
```

Tabela:

```text
i | xi | f(xi) | f'(xi) | x_{i+1} | erro
```

Erro aproximado percentual:

```text
ea = |(x_novo - x_antigo)/x_novo| * 100
```

Exemplo:

```text
f(x) = x³ - x - 2
f'(x) = 3x² - 1
x0 = 1.5
```

Iteração 1:

```text
f(1.5) = 3.375 - 1.5 - 2 = -0.125
f'(1.5) = 3(1.5²) - 1 = 6.75 - 1 = 5.75
```

```text
x1 = 1.5 - (-0.125)/5.75
x1 = 1.521739
```

---

## 17. Método da secante na mão

Não precisa derivada.

Fórmula:

```text
x_{i+1} = x_i - f(x_i)*(x_{i-1} - x_i)/(f(x_{i-1}) - f(x_i))
```

Tabela:

```text
i | x_{i-1} | x_i | f(x_{i-1}) | f(x_i) | x_{i+1}
```

Exemplo:

```text
f(x) = x³ - x - 2
x0 = 1
x1 = 2
```

```text
f(1) = -2
f(2) = 4
```

```text
x2 = 2 - 4*(1 - 2)/(-2 - 4)
x2 = 1.3333
```

---

## 18. Tabela de Routh-Hurwitz

Serve para avaliar estabilidade sem calcular diretamente todas as raízes.

Para um polinômio:

```text
a_n s^n + a_{n-1} s^{n-1} + ... + a_1 s + a_0
```

Condição principal:

```text
Se todos os termos da primeira coluna da tabela de Routh tiverem o mesmo sinal,
o sistema não tem raízes com parte real positiva.
```

Se houver mudança de sinal na primeira coluna:

```text
número de mudanças de sinal = número de raízes no semiplano direito
```

---

## 19. Routh-Hurwitz para polinômio de 4ª ordem

Para:

```text
P(s) = a4 s⁴ + a3 s³ + a2 s² + a1 s + a0
```

Monte:

```text
s⁴ | a4    a2    a0
s³ | a3    a1    0
s² | b1    b2    0
s¹ | c1    0     0
s⁰ | d1    0     0
```

Onde:

```text
b1 = (a3*a2 - a4*a1)/a3
b2 = (a3*a0 - a4*0)/a3 = a0

c1 = (b1*a1 - a3*b2)/b1

d1 = a0
```

A primeira coluna é:

```text
a4, a3, b1, c1, d1
```

---

## 20. Exemplo de tabela de Routh-Hurwitz

Polinômio:

```text
P(s) = 2s⁴ + 5s³ + 12s² + 8s + 8
```

Coeficientes:

```text
a4 = 2
a3 = 5
a2 = 12
a1 = 8
a0 = 8
```

Tabela inicial:

```text
s⁴ | 2    12    8
s³ | 5     8    0
```

Linha s²:

```text
b1 = (5*12 - 2*8)/5 = (60 - 16)/5 = 44/5 = 8.8
b2 = 8
```

Linha s¹:

```text
c1 = (8.8*8 - 5*8)/8.8
c1 = (70.4 - 40)/8.8
c1 = 30.4/8.8
c1 = 3.4545
```

Linha s⁰:

```text
d1 = 8
```

Tabela:

```text
s⁴ | 2       12    8
s³ | 5        8    0
s² | 8.8      8    0
s¹ | 3.4545   0    0
s⁰ | 8        0    0
```

Primeira coluna:

```text
2, 5, 8.8, 3.4545, 8
```

Todos positivos.

Conclusão:

```text
Não há mudança de sinal na primeira coluna.
Logo, não há raízes com parte real positiva.
```

---

## 21. Bairstow — ideia da tabela b e c

Bairstow é usado para raízes de polinômios de grau alto.

A ideia é dividir o polinômio por um fator quadrático:

```text
x² - r x - s
```

Ao achar r e s, as raízes desse fator são:

```text
x = (r ± sqrt(r² + 4s))/2
```

Se o termo dentro da raiz for negativo:

```text
x = r/2 ± i*sqrt(|r²+4s|)/2
```

---

## 22. Bairstow — tabela de recorrência

Considere:

```text
P(x) = a_n x^n + a_{n-1} x^{n-1} + ... + a_1 x + a_0
```

Coeficientes em ordem decrescente:

```text
a_n, a_{n-1}, ..., a_1, a_0
```

Para usar uma forma simples de tabela, considere os índices da esquerda para a direita.

Recorrência para b:

```text
b_n     = a_n
b_{n-1} = a_{n-1} + r*b_n
b_i     = a_i + r*b_{i+1} + s*b_{i+2}
```

Recorrência para c:

```text
c_n     = b_n
c_{n-1} = b_{n-1} + r*c_n
c_i     = b_i + r*c_{i+1} + s*c_{i+2}
```

Depois calcula-se correções para r e s por sistema 2x2.

Na prova, se for pedir código, use o kit. Se pedir raciocínio, mostre que:

```text
1. Monta os coeficientes.
2. Escolhe r e s iniciais.
3. Calcula tabela b.
4. Calcula tabela c.
5. Corrige r e s.
6. Repete até erro pequeno.
7. Extrai raízes do fator quadrático.
8. Reduz o polinômio e repete.
```

---

## 23. Exemplo simples de raízes de polinômio por fatoração

```text
P(x) = x³ - 6x² + 11x - 6
```

Testar raízes inteiras possíveis:

```text
±1, ±2, ±3, ±6
```

Testar x = 1:

```text
P(1) = 1 - 6 + 11 - 6 = 0
```

Então:

```text
x = 1 é raiz.
```

Dividir por (x - 1) usando Horner:

```text
        1   -6    11   -6
1           1    -5     6
        -------------------
        1   -5     6    0
```

Quociente:

```text
x² - 5x + 6
```

Fatorando:

```text
x² - 5x + 6 = (x - 2)(x - 3)
```

Raízes:

```text
x = 1, 2, 3
```

---

## 24. Como saber se raiz complexa vem em par conjugado

Se o polinômio tem coeficientes reais:

```text
raiz complexa:     λ = σ + iω
outra raiz:        λ = σ - iω
```

Exemplo:

```text
λ = -2 + 3i
λ = -2 - 3i
```

---

## 25. Interpretação física das raízes em sistemas dinâmicos

Para sistemas massa-mola-amortecedor:

```text
λ = σ ± iωd
```

Interpretação:

```text
σ < 0    -> resposta decai no tempo -> estável
σ = 0    -> oscilação permanente -> marginal
σ > 0    -> resposta cresce -> instável
```

Se:

```text
λ real negativo
```

então:

```text
decaimento exponencial sem oscilação
```

Se:

```text
λ real positivo
```

então:

```text
crescimento exponencial / instabilidade
```

---

## 26. Checklist para problema manuscrito

Antes de calcular:

```text
1. O problema pede determinante, raiz, polinômio ou estabilidade?
2. Se for sistema dinâmico, montar Mλ² + Cλ + K.
3. Se for 2x2, usar det = a11*a22 - a12*a21.
4. Expandir com calma e agrupar potências.
5. Montar tabela de coeficientes.
6. Se pedir estabilidade, montar Routh-Hurwitz.
7. Se pedir raízes numéricas, escolher método:
   - Bisseção/falsa posição se tiver intervalo.
   - Newton se tiver derivada.
   - Secante se não quiser derivada.
   - Bairstow se for polinômio de grau alto.
```

---

## 27. Checklist para transformar em código

Depois de fazer no papel:

```text
1. Defina os coeficientes ou a função.
2. Escolha o método do kit.
3. Coloque tolerância.
4. Rode.
5. Imprima resultado.
6. Se necessário, faça gráfico.
```

Exemplo:

```python
from kit_prova_calculo_numerico import bissecao

def f(x):
    return x**3 - x - 2

raiz = bissecao(f, 1, 2)
print("raiz =", raiz)
print("f(raiz) =", f(raiz))
```

---

## 28. Fórmulas rápidas

```text
Determinante 2x2:
D = ad - bc

Determinante 3x3:
D = a11(a22a33-a23a32) - a12(a21a33-a23a31) + a13(a21a32-a22a31)

Polinômio característico:
P(λ) = det(Mλ² + Cλ + K)

Bisseção:
xr = (a+b)/2

Falsa posição:
xr = b - f(b)(a-b)/(f(a)-f(b))

Newton:
x_{i+1} = x_i - f(x_i)/f'(x_i)

Secante:
x_{i+1} = x_i - f(x_i)(x_{i-1}-x_i)/(f(x_{i-1})-f(x_i))

Erro aproximado:
ea = |(novo - antigo)/novo| * 100

Bhaskara:
x = (-b ± sqrt(b² - 4ac))/(2a)

Bairstow:
fator quadrático = x² - r x - s
raízes = (r ± sqrt(r² + 4s))/2
```
