# Exercício resolvido — equações de diferenças

Considere uma malha quadrada com `Delta x = Delta y = Delta`.

## 1. Nó interno

A equação de Laplace discretizada é:

```math
\frac{T_E-2T_P+T_W}{\Delta^2}
+\frac{T_N-2T_P+T_S}{\Delta^2}=0.
```

Portanto:

```math
\boxed{4T_P-T_E-T_W-T_N-T_S=0}.
```

## 2. Nó na extremidade direita convectiva

A condição é:

```math
-k\frac{\partial T}{\partial x}=h(T_P-T_\infty).
```

Com nó fictício a leste e derivada centrada, elimina-se `T_E`. Definindo:

```math
Bi_\Delta=\frac{h\Delta}{k},
```

resulta:

```math
\boxed{(4+2Bi_\Delta)T_P-2T_W-T_N-T_S
=2Bi_\Delta T_\infty}.
```

## 3. Canto superior direito

Há convecção simultânea nas direções `x` e `y`. A eliminação dos dois nós fictícios fornece:

```math
\boxed{(4+4Bi_\Delta)T_P-2T_W-2T_S
=4Bi_\Delta T_\infty}.
```

## Interpretação

Os coeficientes duplicados dos vizinhos internos aparecem porque cada condição de Robin substitui um nó externo por uma expressão que envolve o vizinho oposto e a temperatura ambiente.
