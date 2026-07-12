# PPC6 — Condução 2D em Aleta Retangular por Diferenças Finitas

Este repositório contém a implementação computacional do **PPC6 de Cálculo Numérico Aplicado**.

O programa resolve numericamente o problema bidimensional de condução de calor em regime permanente em uma aleta retangular de seção transversal constante. A base da aleta possui temperatura prescrita e as demais superfícies trocam calor por convecção com o ambiente.

A equação governante é:

```math
\frac{\partial^2 T}{\partial x^2}+\frac{\partial^2 T}{\partial y^2}=0
```

com as condições de contorno:

```math
T(0,y)=T_b
```

na base da aleta, e

```math
-k\frac{\partial T}{\partial n}=h(T-T_\infty)
```

nas superfícies superior, inferior e na extremidade livre.

## Métodos implementados

1. **Diferenças finitas**
   - discretiza a equação de Laplace em uma malha retangular uniforme;
   - usa diferenças centradas para os nós internos;
   - usa nós fictícios para aplicar a condição convectiva de Robin.

2. **Eliminação de Gauss**
   - monta o sistema linear completo `[A]{T}={b}`;
   - resolve o sistema por eliminação direta;
   - usa pivoteamento parcial para melhorar a estabilidade numérica.

3. **Método de Liebmann**
   - resolve iterativamente as equações de diferenças;
   - corresponde ao método de Gauss-Seidel aplicado à malha;
   - usa imediatamente os valores recém-atualizados.

4. **Método de Liebmann com relaxação**
   - implementa a forma SOR;
   - atualiza a temperatura por:

```math
T^{novo}=(1-\omega)T^{velho}+\omega T^{calculado}
```

5. **Comparação com solução analítica 1D**
   - extrai a temperatura ao longo da linha central da aleta;
   - compara com a solução analítica clássica de aleta com ponta convectiva;
   - calcula o erro percentual médio.

## Arquivos

| Arquivo | Descrição |
|---|---|
| `ppc6_aleta.py` | Código principal do PPC6 |
| `temperature_gauss.dat` | Campo de temperatura obtido por eliminação de Gauss |
| `temperature_liebmann.dat` | Campo de temperatura obtido por Liebmann sem relaxação |
| `temperature_sor.dat` | Campo de temperatura obtido por Liebmann com relaxação |
| `linha_central.dat` | Comparação entre solução numérica 2D e solução analítica 1D |
| `estudo_relaxacao.dat` | Efeito do fator de relaxação `omega` na convergência |
| `estudo_malha.dat` | Estudo de refinamento de malha |
| `resumo_resultados.txt` | Resumo dos resultados principais |
| `mapa_temperatura.png` | Mapa de temperatura bidimensional |
| `contornos_isotermicos.png` | Curvas de nível da temperatura |
| `linha_central.png` | Temperatura na linha central da aleta |

Os arquivos `.dat`, `.txt` e `.png` são gerados automaticamente durante a execução.

## Como executar

Na pasta do projeto:

```bash
python ppc6_aleta.py
```

A execução padrão usa:

```text
L = 0.12 m
H = 0.02 m
k = 205 W/(m.K)
h = 25 W/(m².K)
T_b = 100 °C
T_inf = 25 °C
nx = 21
ny = 9
tolerância = 1e-8
omega = 1.60
```

## Informar parâmetros pela linha de comando

```bash
python ppc6_aleta.py \
  --L 0.12 \
  --H 0.02 \
  --k 205 \
  --h 25 \
  --Tb 100 \
  --Tinf 25 \
  --nx 21 \
  --ny 9 \
  --tol 1e-8 \
  --omega 1.6
```

## Executar no modo interativo

```bash
python ppc6_aleta.py --interativo
```

Nesse modo, o programa solicita os dados pelo terminal. Pressionar `Enter` mantém o valor padrão indicado.

## Salvar os resultados em outra pasta

```bash
python ppc6_aleta.py --outdir resultados
```

## Executar sem gerar gráficos

```bash
python ppc6_aleta.py --sem-plots
```

## Pular os estudos de relaxação e malha

```bash
python ppc6_aleta.py --pular-estudos
```

Essa opção é útil para testes rápidos.

## Consultar todas as opções

```bash
python ppc6_aleta.py --help
```

## Dependências

A parte numérica principal usa apenas a biblioteca padrão do Python.

Recomenda-se **Python 3.9 ou superior**.

Para gerar gráficos, instale o `matplotlib`:

```bash
pip install matplotlib
```

Caso o `matplotlib` não esteja instalado, o cálculo será executado normalmente, mas os arquivos PNG não serão gerados.

## Formato dos arquivos `.dat`

### `temperature_gauss.dat`, `temperature_liebmann.dat` e `temperature_sor.dat`

```text
x y T
```

onde:

- `x` é a coordenada horizontal do nó;
- `y` é a coordenada vertical do nó;
- `T` é a temperatura calculada.

### `linha_central.dat`

```text
x T_numerico T_analitico erro_percentual
```

onde:

- `T_numerico` é a temperatura obtida pela solução bidimensional;
- `T_analitico` é a temperatura da solução unidimensional clássica;
- `erro_percentual` é o erro relativo local.

### `estudo_relaxacao.dat`

```text
omega iteracoes erro_final tempo_s convergiu
```

### `estudo_malha.dat`

```text
nx ny iteracoes erro_final tempo_s erro_medio_percentual convergiu
```

## Solução analítica 1D usada na comparação

A distribuição adimensional clássica para uma aleta de seção constante com extremidade convectiva é:

```math
\theta(x)=\frac{T(x)-T_\infty}{T_b-T_\infty}
```

com:

```math
m=\sqrt{\frac{hP}{kA_c}}
```

A solução usada é:

```math
\theta(x)=
\frac{
\cosh[m(L-x)] + \frac{h}{mk}\sinh[m(L-x)]
}{
\cosh(mL)+\frac{h}{mk}\sinh(mL)
}.
```

No código, considera-se largura unitária fora do plano, de modo que a comparação analítica representa uma aproximação unidimensional compatível com a troca convectiva pelas faces superior e inferior.

## Resultado esperado com os parâmetros padrão

A execução padrão produz resultados próximos de:

```text
Gauss: tempo aproximado de 0.10 s
Liebmann: aproximadamente 11012 iterações
SOR com omega = 1.60: aproximadamente 3097 iterações
Erro percentual médio na linha central: aproximadamente 0.0218%
```

Pequenas diferenças podem ocorrer por causa da máquina utilizada e da versão do Python.

## Observações numéricas

O método direto de Gauss fornece uma solução de referência para a malha montada. O método de Liebmann converge para a mesma solução, mas exige muitas iterações. A sobre-relaxação reduz significativamente o número de iterações, desde que o valor de `omega` seja escolhido dentro de uma faixa estável.

Para problemas maiores, a matriz de coeficientes se torna esparsa. O código monta a matriz densa apenas para fins didáticos e para mostrar claramente a relação entre diferenças finitas e sistema linear. Em aplicações reais de maior porte, seria mais eficiente usar armazenamento esparso.


