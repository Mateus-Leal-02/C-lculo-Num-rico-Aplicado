# PPC5 — Equação de Blasius pelo Método do Tiro

Este diretório contém a implementação computacional do **PPC5 de Cálculo Numérico Aplicado**.

O programa resolve numericamente a equação de Blasius:

```math
f'''(\eta)+\frac{1}{2}f(\eta)f''(\eta)=0
```

sujeita às condições de contorno:

```math
f(0)=0,\qquad f'(0)=0,\qquad f'(\infty)=1.
```

Como a condição sobre `f'` é especificada longe da parede, o problema é tratado como um **Problema de Valor de Contorno (PVC)**. O valor desconhecido

```math
s=f''(0)
```

é ajustado até que:

```math
\left|f'(\eta_{\max})-1\right|<\varepsilon.
```

## Métodos implementados

1. **Método do Tiro**
   - transforma o PVC em uma sequência de Problemas de Valor Inicial;
   - utiliza `s = f''(0)` como parâmetro de ajuste;
   - define a função erro `E(s) = f'(η_max; s) - 1`.

2. **Runge–Kutta de quarta ordem**
   - integra simultaneamente `f`, `f'` e `f''`;
   - os 12 coeficientes do RK4 são calculados diretamente no código;
   - não são utilizadas bibliotecas de integração numérica pronta.

3. **Newton–Raphson**
   - atualiza os sucessivos chutes de `s`;
   - a derivada de `E(s)` é aproximada por diferença finita progressiva.

## Sistema de primeira ordem

São definidas as variáveis:

```math
y_1=f,\qquad y_2=f',\qquad y_3=f''.
```

O sistema integrado pelo programa é:

```math
\frac{dy_1}{d\eta}=y_2,
```

```math
\frac{dy_2}{d\eta}=y_3,
```

```math
\frac{dy_3}{d\eta}=-\frac{1}{2}y_1y_3.
```

As condições iniciais do PVI equivalente são:

```math
y_1(0)=0,\qquad y_2(0)=0,\qquad y_3(0)=s.
```

## Arquivos

| Arquivo | Descrição |
|---|---|
| `ppc5_blasius.py` | Código principal do PPC5 |
| `blasius_solution.dat` | Valores de `η`, `f`, `f'` e `f''` |
| `shooting_log.dat` | Histórico das iterações do Método do Tiro |
| `perfil_f.png` | Gráfico de `f(η)` |
| `perfil_fp.png` | Gráfico de `f'(η) = u/U∞` |
| `perfil_fpp.png` | Gráfico de `f''(η)` |

Os arquivos `.dat` e `.png` são criados durante a execução.

## Como executar

Na pasta do projeto:

```bash
python ppc5_blasius.py
```

A execução padrão utiliza:

```text
s0 = 0.3
Delta eta = 0.01
eta_max = 10
tolerância = 1e-10
máximo de iterações = 50
Re_x = 1e5
```

### Informar parâmetros pela linha de comando

```bash
python ppc5_blasius.py \
  --s0 0.3 \
  --deta 0.01 \
  --eta-max 10 \
  --tol 1e-10 \
  --max-iter 50 \
  --re-x 1e5
```

### Executar no modo interativo

```bash
python ppc5_blasius.py --interativo
```

Nesse modo, o programa solicita os parâmetros pelo terminal. Pressionar `Enter` mantém o valor padrão indicado.

### Salvar os resultados em outra pasta

```bash
python ppc5_blasius.py --outdir resultados
```

### Executar sem gerar gráficos

```bash
python ppc5_blasius.py --sem-plots
```

### Consultar todas as opções

```bash
python ppc5_blasius.py --help
```

## Dependências

A solução numérica e os arquivos `.dat` utilizam apenas a biblioteca padrão do Python.

Recomenda-se **Python 3.9 ou superior**.

Para gerar os gráficos, instale o `matplotlib`:

```bash
pip install matplotlib
```

Caso o `matplotlib` não esteja instalado, o cálculo ainda será executado normalmente, mas os arquivos PNG não serão gerados.

## Formato dos arquivos `.dat`

### `blasius_solution.dat`

```text
eta f fp fpp
```

onde:

- `eta` é a variável de similaridade;
- `f` é a função de Blasius;
- `fp` representa `f'(η) = u/U∞`;
- `fpp` representa `f''(η)`.

### `shooting_log.dat`

```text
iter s fp_eta_max erro
```

onde:

- `iter` é o número da iteração do Método do Tiro;
- `s` é o chute corrente de `f''(0)`;
- `fp_eta_max` é o valor de `f'(η_max)`;
- `erro` é `f'(η_max) - 1`.

## Grandezas calculadas

Além da solução da equação de Blasius, o programa calcula:

### Coeficiente local de atrito

```math
C_f=\frac{2f''(0)}{\sqrt{Re_x}}.
```

### Posição adimensional de 99% da velocidade externa

```math
f'(\eta_{99})=0.99.
```

O valor de `η_99` é determinado por interpolação linear entre os pontos da malha.

### Coeficiente de espessura da camada limite

```math
\frac{\delta}{x}=\frac{C_\delta}{\sqrt{Re_x}},
\qquad
C_\delta=\eta_{99}.
```

O resultado numérico é comparado com a correlação clássica:

```math
\frac{\delta}{x}=\frac{4.92}{\sqrt{Re_x}}.
```

## Resultado esperado

Com os parâmetros padrão, o valor convergido deve ser aproximadamente:

```math
f''(0)=0.332057337205.
```

Também são esperados:

```math
\eta_{99}\approx 4.909989486136,
```

```math
C_f\approx 2.100114998676\times 10^{-3}
\quad \text{para}\quad Re_x=10^5.
```

O Método do Tiro converge em aproximadamente quatro iterações para o chute inicial `s0 = 0.3`.

## Exemplo de saída

```text
Resultados finais
-----------------
Convergência             = sim
f''(0) convergido        = 0.332057337205
Iterações do tiro        = 4
Erro final               = 4.742872761199e-13
f'(eta_max)              = 1.000000000000
eta_99                   = 4.909989486136
C_delta = eta_99         = 4.909989486136
C_f para Re_x informado  = 2.100114998676e-03
```

Pequenas diferenças nos últimos algarismos podem ocorrer conforme a versão do Python e os parâmetros numéricos adotados.
