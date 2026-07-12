# PPC5 — Equação de Blasius pelo Método do Tiro

## Problema

```math
f'''(\eta)+\frac12 f(\eta)f''(\eta)=0
```

com:

```math
f(0)=0,\qquad f'(0)=0,\qquad f'(\infty)=1.
```

Fisicamente:

```math
f'(\eta)=\frac{u}{U_\infty}.
```

## Transformação em sistema

```math
y_1=f,\qquad y_2=f',\qquad y_3=f''
```

```math
\frac{dy_1}{d\eta}=y_2,
\qquad
\frac{dy_2}{d\eta}=y_3,
\qquad
\frac{dy_3}{d\eta}=-\frac12y_1y_3.
```

O parâmetro desconhecido é:

```math
s=f''(0).
```

## Métodos

- Método do Tiro;
- RK4 implementado explicitamente para três EDOs;
- Newton–Raphson para corrigir `s`;
- derivada da função erro por diferença finita;
- interpolação para determinar `eta_99`.

## Arquivos

| Arquivo | Descrição |
|---|---|
| `Blasius.py` | código principal existente no GitHub |
| `docs/APC5.pdf` | formulação e coeficientes do RK4 |
| `docs/PPC5.pdf` | relatório da implementação |
| `exercicios/exercicio-resolvido.md` | primeiro passo do RK4 |
| `desafio.md` | comparação entre Newton e secante |

## Dependências

A solução numérica usa a biblioteca padrão. Para gráficos:

```bash
pip install matplotlib
```

## Como executar

```bash
python Blasius.py
```

## Resultado de referência

```math
f''(0)\approx0,332057337205.
```

```math
\eta_{99}\approx4,909989486.
```

```math
C_f=\frac{2f''(0)}{\sqrt{Re_x}}.
```

## Validação

- condição `f'(eta_max) ≈ 1`;
- valor clássico de `f''(0)`;
- monotonicidade de `f'`;
- aproximação assintótica de `f''` para zero;
- estudo com `Delta eta` e `eta_max`.

## Fontes de erro

- truncamento do domínio infinito;
- passo do RK4;
- tolerância do Método do Tiro;
- derivada numérica de `E(s)`;
- interpolação de `eta_99`.

## Material complementar

- [Exercício resolvido](exercicios/exercicio-resolvido.md)
- [Desafio](desafio.md)
- [APC5](docs/APC5.pdf)
- [PPC5](docs/PPC5.pdf)
