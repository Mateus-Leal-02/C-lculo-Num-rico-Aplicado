# PPC3 — Difusão de Calor Transiente em Pastilha de UO₂

Este diretório contém a implementação computacional do **PPC3 de Cálculo Numérico Aplicado**.

O problema resolvido é a condução de calor unidimensional transiente em uma placa plana simétrica, representando uma idealização da meia-espessura de uma pastilha combustível de dióxido de urânio (`UO₂`) submetida a geração volumétrica interna e resfriamento convectivo na superfície.

## Modelo físico

A equação governante é:

```math
\rho C_p \frac{\partial T}{\partial t}
=
k \frac{\partial^2 T}{\partial x^2}
+
\dot{q}
```

com as condições de contorno:

```math
\left.\frac{\partial T}{\partial x}\right|_{x=0}=0
```

```math
-k\left.\frac{\partial T}{\partial x}\right|_{x=L}
=
h(T_L-T_\infty)
```

## Método numérico

A solução usa:

- diferenças finitas no espaço;
- esquema implícito no tempo;
- montagem de um sistema tridiagonal;
- solução do sistema pelo Algoritmo de Thomas, também chamado de TDMA.

O esquema implícito foi usado porque é estável para passos de tempo maiores que os normalmente permitidos em formulações explícitas.

## Arquivos

| Arquivo | Descrição |
|---|---|
| `ppc3_reator_implicito.py` | Código principal do PPC3 |
| `ppc3_perfis_transientes.dat` | Perfis de temperatura em diferentes tempos |
| `ppc3_validacao_sem_geracao.dat` | Comparação entre solução numérica e solução analítica sem geração |
| `ppc3_perfis_transientes.png` | Figura dos perfis transientes, gerada se `matplotlib` estiver instalado |
| `ppc3_validacao.png` | Figura de validação, gerada se `matplotlib` estiver instalado |

## Como executar

Na pasta deste PPC:

```bash
python ppc3_reator_implicito.py
```

Para alterar a malha, o passo de tempo e o tempo final:

```bash
python ppc3_reator_implicito.py --N 81 --dt 0.005 --t-final 20
```

Para executar sem gerar figuras:

```bash
python ppc3_reator_implicito.py --sem-plots
```

Para salvar os arquivos de saída em outra pasta:

```bash
python ppc3_reator_implicito.py --outdir resultados
```

## Dependências

Obrigatória:

```bash
pip install numpy
```

Opcional para geração de gráficos:

```bash
pip install matplotlib
```

## Parâmetros adotados

| Parâmetro | Valor | Unidade |
|---|---:|---|
| `L` | `5.0e-3` | m |
| `k` | `4.0` | W/(m·K) |
| `rho` | `10500.0` | kg/m³ |
| `Cp` | `300.0` | J/(kg·K) |
| `h` | `30000.0` | W/(m²·K) |
| `T_inf` | `300.0` | °C |
| `q_dot` | `3.0e8` | W/m³ |

## Resultado esperado

O programa calcula o aquecimento transiente da pastilha e compara o resultado numérico, no caso sem geração interna, com a solução analítica por série. No caso com geração, o perfil tende ao regime permanente parabólico esperado para condução unidimensional com geração volumétrica.
