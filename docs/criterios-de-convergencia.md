# Critérios de convergência

Um método numérico deve possuir um critério de parada explícito e verificável.

## Incremento entre iterações

```math
\|x^{(k+1)}-x^{(k)}\| < \varepsilon
```

É simples, mas um incremento pequeno não garante sozinho que a equação foi satisfeita.

## Resíduo

Para um sistema `A x = b`:

```math
r^{(k)} = b-Ax^{(k)}
```

Um critério robusto é:

```math
\|r^{(k)}\|_\infty < \varepsilon_r
```

## Erro da condição de contorno

No Método do Tiro do PPC5:

```math
E(s)=f'(\eta_{max};s)-1
```

A convergência é aceita quando `|E(s)| < tolerância`.

## Norma do gradiente

Em otimização sem restrições:

```math
\|\nabla f(x_k)\|_2 < \varepsilon_g
```

A condição indica proximidade de um ponto estacionário. A Hessiana ou a estrutura do problema determina se ele é máximo, mínimo ou sela.

## Critério em métodos de malha

No Liebmann/SOR, pode-se usar:

```math
\max_{i,j}|T_{i,j}^{novo}-T_{i,j}^{velho}| < \varepsilon
```

Também é recomendável calcular o resíduo discreto da equação de Laplace.

## Máximo de iterações

Todo método iterativo deve possuir um limite de segurança. Ao atingir o limite, o programa deve informar que não convergiu, em vez de apresentar o último valor como solução garantida.

## Convergência de malha

Compare sucessivos refinamentos e uma grandeza de interesse `Q`:

```math
\frac{|Q_{fina}-Q_{grossa}|}{|Q_{fina}|} < \varepsilon_m
```

A tolerância iterativa deve ser suficientemente menor que o erro de discretização; caso contrário, o estudo de malha mistura duas fontes de erro.
