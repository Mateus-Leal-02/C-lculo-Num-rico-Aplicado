# Exercício resolvido — Algoritmo de Thomas

Resolva:

```math
\begin{aligned}
2x_1-x_2&=1,\\
-x_1+2x_2-x_3&=0,\\
-x_2+2x_3&=1.
\end{aligned}
```

As diagonais são:

```text
a = [0, -1, -1]
b = [2,  2,  2]
c = [-1, -1, 0]
d = [1,  0,  1]
```

## Varredura direta

Primeira linha:

```math
c'_1=\frac{-1}{2}=-0,5,
\qquad
d'_1=\frac{1}{2}=0,5.
```

Segunda linha:

```math
m_2=b_2-a_2c'_1=2-(-1)(-0,5)=1,5
```

```math
c'_2=\frac{-1}{1,5}=-\frac{2}{3}
```

```math
d'_2=\frac{0-(-1)(0,5)}{1,5}=\frac{1}{3}.
```

Terceira linha:

```math
m_3=2-(-1)\left(-\frac{2}{3}\right)=\frac{4}{3}
```

```math
d'_3=\frac{1-(-1)(1/3)}{4/3}=1.
```

## Substituição regressiva

```math
x_3=d'_3=1
```

```math
x_2=d'_2-c'_2x_3=\frac13+\frac23=1
```

```math
x_1=d'_1-c'_1x_2=0,5+0,5=1.
```

Logo:

```math
\boxed{x=[1,1,1]^T}
```

A substituição no sistema original fornece resíduo nulo.
