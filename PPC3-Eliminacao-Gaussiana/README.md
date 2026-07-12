# PPC3 — Difusão de Calor Transiente em Pastilha de UO₂

> A pasta mantém o nome histórico `PPC3-Eliminacao-Gaussiana`, mas o método principal implementado é o **Algoritmo de Thomas (TDMA)** para sistemas tridiagonais.

## Problema físico

O programa modela a condução de calor unidimensional transiente em uma placa plana simétrica, representando a meia-espessura de uma pastilha combustível de dióxido de urânio com geração volumétrica e convecção na superfície.

```math
\rho C_p\frac{\partial T}{\partial t}
=k\frac{\partial^2T}{\partial x^2}+\dot q.
```

Condições:

```math
\left.\frac{\partial T}{\partial x}\right|_{x=0}=0
```

```math
-k\left.\frac{\partial T}{\partial x}\right|_{x=L}=h(T_L-T_\infty).
```

## Método numérico

- diferenças finitas no espaço;
- esquema totalmente implícito no tempo;
- tratamento específico dos nós de simetria e convecção;
- solução do sistema tridiagonal pelo TDMA.

## Arquivos

| Arquivo | Descrição |
|---|---|
| `Reator.py` | código principal existente |
| `docs/APC3.pdf` | exercícios analíticos de sistemas lineares |
| `docs/PPC3.pdf` | relatório do modelo transiente |
| `exercicios/exercicio-resolvido.md` | solução manual de um sistema pelo TDMA |
| `desafio.md` | implementação de Crank–Nicolson |

## Dependências

```bash
pip install numpy matplotlib
```

## Como executar

```bash
python Reator.py
```

Exemplo com parâmetros numéricos:

```bash
python Reator.py --N 81 --dt 0.005 --t-final 20
```

Sem figuras:

```bash
python Reator.py --sem-plots
```

## Parâmetros físicos usados no caso padrão

| Parâmetro | Valor | Unidade |
|---|---:|---|
| `L` | `5,0e-3` | m |
| `k` | `4,0` | W/(m·K) |
| `rho` | `10500` | kg/m³ |
| `Cp` | `300` | J/(kg·K) |
| `h` | `30000` | W/(m²·K) |
| `T_inf` | `300` | °C |
| `q_dot` | `3,0e8` | W/m³ |

## Validação

No caso sem geração interna, o resultado numérico é comparado com a solução em série para parede plana com simetria e convecção. No caso com geração, o perfil deve tender à solução parabólica permanente.

## Fontes de erro

- discretização espacial;
- discretização temporal;
- truncamento da série analítica;
- representação unidimensional;
- propriedades constantes.

## Material complementar

- [Exercício resolvido](exercicios/exercicio-resolvido.md)
- [Desafio](desafio.md)
- [APC3](docs/APC3.pdf)
- [PPC3](docs/PPC3.pdf)
