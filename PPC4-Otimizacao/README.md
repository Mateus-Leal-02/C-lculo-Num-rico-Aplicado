# PPC4 — Otimização Bidimensional sem Restrições

Este diretório contém a implementação computacional do **PPC4 de Cálculo Numérico Aplicado**.

O problema consiste em maximizar a função bidimensional:

```math
f(x,y)=2xy+2x-x^2-2y^2
```

a partir de um ponto inicial, usando dois métodos de otimização baseados em gradiente.

## Métodos implementados

1. **Aclive máximo**
   - usa a direção do gradiente como direção de busca;
   - em cada iteração, calcula o melhor passo na direção escolhida.

2. **Gradientes conjugados de Fletcher-Reeves**
   - combina o gradiente atual com a direção anterior;
   - tende a convergir em menos iterações para funções quadráticas.

## Arquivos

| Arquivo | Descrição |
|---|---|
| `Otimizacao.py` | Código principal do PPC4 |
| `output1.dat` | Histórico iterativo do método do aclive máximo |
| `output2.dat` | Histórico iterativo do método de Fletcher-Reeves |
| `function.dat` | Malha de valores `x`, `y`, `f(x,y)` para curvas de nível |
| `ppc4_caminhos.png` | Figura dos caminhos iterativos, gerada se `matplotlib` estiver instalado |

## Como executar

Na pasta deste PPC:

```bash
python Otimizacao.py
```

Por padrão, o código usa:

```text
x0 = -2
y0 = 3
tol = 1e-8
```

Para alterar o ponto inicial:

```bash
python Otimizacao.py --x0 -2 --y0 3
```

Para alterar a tolerância e o número máximo de iterações:

```bash
python Otimizacao.py --tol 1e-10 --max-iter 200
```

Para executar sem gerar figura:

```bash
python Otimizacao.py --sem-plots
```

Para salvar os arquivos em outra pasta:

```bash
python Otimizacao.py --outdir resultados
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

## Formato dos arquivos `.dat`

Os arquivos `output1.dat` e `output2.dat` seguem o formato:

```text
iter erro h x y dfdx dfdy f
```

onde:

- `iter` é o número da iteração;
- `erro` é a norma do gradiente;
- `h` é o passo ótimo na direção de busca;
- `x` e `y` são as coordenadas do ponto atual;
- `dfdx` e `dfdy` são as componentes do gradiente;
- `f` é o valor da função objetivo.

O arquivo `function.dat` segue o formato:

```text
x y f
```

## Resultado esperado

O ponto ótimo analítico é:

```math
(x^*,y^*)=(2,1)
```

com:

```math
f(x^*,y^*)=2
```

O método de Fletcher-Reeves converge rapidamente para esse ponto por se tratar de uma função quadrática côncava.

