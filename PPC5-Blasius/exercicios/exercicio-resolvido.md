# Exercício resolvido — primeiro passo do RK4

Considere:

```text
s = f''(0) = 0,332057337205
Delta eta = 0,1
f(0) = 0
f'(0) = 0
f''(0) = s
```

O sistema é:

```math
y_1'=y_2,\qquad y_2'=y_3,\qquad y_3'=-\frac12y_1y_3.
```

Usando coeficientes já multiplicados pelo passo:

```text
k1 = [0;
      0,0332057337;
      0]
```

```text
k2 ≈ [0,0016602867;
      0,0332057337;
      0]
```

```text
k3 ≈ [0,0016602867;
      0,0332057337;
     -0,0000137828]
```

```text
k4 ≈ [0,0033205734;
      0,0332043554;
     -0,0000275644]
```

A atualização:

```math
\mathbf y_1=\mathbf y_0+\frac16
(\mathbf k_1+2\mathbf k_2+2\mathbf k_3+\mathbf k_4)
```

fornece aproximadamente:

```math
\boxed{f(0,1)\approx0,00166029}
```

```math
\boxed{f'(0,1)\approx0,03320550}
```

```math
\boxed{f''(0,1)\approx0,33204815}
```

A pequena redução de `f''` é coerente com a equação `f'''=-ff''/2`, pois `f` se torna positivo após o início da integração.
