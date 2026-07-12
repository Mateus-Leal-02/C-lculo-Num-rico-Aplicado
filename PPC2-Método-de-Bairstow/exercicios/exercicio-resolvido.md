# Exercício resolvido — uma iteração de Bairstow

Considere:

```math
P(x)=x^4-5x^2+4
```

com chute inicial:

```math
r_0=0,1,\qquad s_0=0,9.
```

Para o fator `x²-rx-s`, a divisão sintética produz:

```text
b = [1; 0,1; -4,09; -0,319; 0,2871]
```

Assim:

```math
b_1=-0,319,\qquad b_0=0,2871.
```

A segunda recorrência fornece:

```text
c = [1; 0,2; -3,17; -0,456]
```

Usando a notação do método:

```math
c_1=-0,456,\qquad c_2=-3,17,\qquad c_3=0,2.
```

O determinante é:

```math
\det=c_2^2-c_1c_3=10,1401.
```

As correções são:

```math
\Delta r=\frac{b_0c_3-b_1c_2}{\det}=-0,09406317
```

```math
\Delta s=\frac{b_1c_1-b_0c_2}{\det}=0,10409868.
```

Portanto:

```math
\boxed{r_1=0,00593683}
```

```math
\boxed{s_1=1,00409868}
```

A primeira atualização já aproxima o fator `x²-1`, associado às raízes `±1` do polinômio.
