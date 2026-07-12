# Exercício resolvido — duas iterações do aclive máximo

Para:

```math
f(x,y)=2xy+2x-x^2-2y^2
```

o gradiente é:

```math
\nabla f=(2y+2-2x,\;2x-4y).
```

## Iteração 0

No ponto `P0=(-2,3)`:

```math
\nabla f(P_0)=(12,-16).
```

A direção de aclive é `d0=(12,-16)`. Para a função quadrática, o passo ótimo é:

```math
h_0=-\frac{g_0^Td_0}{d_0^THd_0}=\frac{5}{26}.
```

Logo:

```math
P_1=P_0+h_0d_0
=\left(\frac{4}{13},-\frac{1}{13}\right).
```

## Iteração 1

```math
\nabla f(P_1)=\left(\frac{16}{13},\frac{12}{13}\right).
```

O novo passo ótimo é:

```math
h_1=\frac{5}{4}.
```

Assim:

```math
P_2=P_1+h_1\nabla f(P_1)
=\left(\frac{24}{13},\frac{14}{13}\right).
```

Em decimais:

```text
P1 ≈ (0,307692; -0,076923)
P2 ≈ (1,846154; 1,076923)
```

O caminho alterna direções devido à geometria alongada das curvas de nível. O aclive máximo converge, mas não atinge o ótimo em apenas duas iterações.
