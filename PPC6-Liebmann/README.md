# PPC6 — Condução 2D em Aleta Retangular

## Problema físico

A aleta possui temperatura prescrita na base e troca calor por convecção nas superfícies superior, inferior e na extremidade livre.

```math
\frac{\partial^2T}{\partial x^2}
+\frac{\partial^2T}{\partial y^2}=0.
```

Condições:

```math
T(0,y)=T_b
```

```math
-k\frac{\partial T}{\partial n}=h(T-T_\infty).
```

## Métodos implementados

1. diferenças finitas em malha retangular;
2. condição de Robin por nós fictícios;
3. eliminação de Gauss com pivoteamento parcial;
4. Liebmann, equivalente a Gauss–Seidel na malha;
5. SOR com fator `omega`;
6. estudo de relaxação;
7. refinamento de malha;
8. comparação com a solução analítica unidimensional.

## Arquivos

| Arquivo | Descrição |
|---|---|
| `Liebmann.py` | código principal considerado nesta documentação |
| `docs/APC6.pdf` | derivação das equações discretas |
| `docs/PPC6.pdf` | relatório computacional |
| `exercicios/exercicio-resolvido.md` | equações de nós interno e convectivos |
| `desafio.md` | investigação do fator ótimo de relaxação |

> Caso ainda exista `Laplace.py`, compare os dois códigos e mantenha somente a versão completa para evitar ambiguidade.

## Dependências

```bash
pip install numpy matplotlib
```

## Como executar

```bash
python Liebmann.py
```

## Parâmetros padrão documentados

| Parâmetro | Valor |
|---|---:|
| `L` | `0,12 m` |
| `H` | `0,02 m` |
| `k` | `205 W/(m·K)` |
| `h` | `25 W/(m²·K)` |
| `T_b` | `100 °C` |
| `T_inf` | `25 °C` |
| `nx` | `21` |
| `ny` | `9` |
| `omega` | `1,60` |

## Saídas esperadas

- campos de temperatura por Gauss, Liebmann e SOR;
- histórico ou resumo de convergência;
- linha central comparada com solução 1D;
- estudo de `omega`;
- estudo de malha;
- mapa térmico e isotermas.

## Validação

- equivalência entre os três resolvedores;
- resíduo das equações discretas;
- comparação com a solução 1D;
- comportamento com refinamento de malha;
- temperaturas limitadas entre ambiente e base para o caso sem geração.

## Material complementar

- [Exercício resolvido](exercicios/exercicio-resolvido.md)
- [Desafio](desafio.md)
- [APC6](docs/APC6.pdf)
- [PPC6](docs/PPC6.pdf)
