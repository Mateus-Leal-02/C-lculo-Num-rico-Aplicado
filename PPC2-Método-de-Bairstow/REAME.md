# 🔢 Raízes de Polinômios — Método de Bairstow

> **Programa para Casa #2 · Cálculo Numérico Aplicado**
> Prof. Dr. Rafael Gabler Gontijo · UnB · abril de 2026

---

## 📌 Sobre o projeto

Este repositório contém a implementação do **método de Bairstow** para determinação de todas as raízes (reais e complexas) de polinômios de grau arbitrário. O método é aplicado à análise do sistema massa-mola-amortecedor com 2 graus de liberdade (GDL) estudado na APC2, cujo polinômio característico é:

$$P(\lambda) = 2\lambda^4 + 5\lambda^3 + 12\lambda^2 + 8\lambda + 8$$

Adicionalmente, é gerado o **fractal de Bairstow**: um mapa de convergência no plano $(r_0, s_0)$ que revela a estrutura geométrica das bacias de atração do método.

**Nenhuma biblioteca de fatoração de polinômios é utilizada** — o método é implementado do zero a partir dos conceitos de divisão sintética e Newton-Raphson.

---

## ⚙️ O método de Bairstow

O método de Bairstow divide iterativamente o polinômio $f_n(x)$ por um fator quadrático real $D(x) = x^2 - rx - s$, ajustando $(r, s)$ para zerar o resto da divisão.

### Divisão sintética (coeficientes $b_i$)

Para $f_n(x) = a_n x^n + a_{n-1}x^{n-1} + \cdots + a_0$, a divisão por $D(x)$ produz:

$$f_n(x) = \big(b_2 + b_3 x + \cdots + b_n x^{n-2}\big)(x^2 - rx - s) + \big[b_0 + b_1(x - r)\big]$$

com a recorrência:

$$b_n = a_n, \quad b_{n-1} = a_{n-1} + r\,b_n, \quad b_k = a_k + r\,b_{k+1} + s\,b_{k+2}, \; k = n{-}2, \ldots, 0$$

### Derivadas parciais (coeficientes $c_i$)

Uma segunda divisão sintética sobre os $b_i$ fornece os $c_i$, que representam as derivadas parciais de $b_0$ e $b_1$ em relação a $r$ e $s$:

$$c_n = b_n, \quad c_{n-1} = b_{n-1} + r\,c_n, \quad c_k = b_k + r\,c_{k+1} + s\,c_{k+2}, \; k = n{-}2, \ldots, 1$$

com o mapeamento:
$$c_1 = \frac{\partial b_0}{\partial r}, \qquad c_2 = \frac{\partial b_0}{\partial s} = \frac{\partial b_1}{\partial r}, \qquad c_3 = \frac{\partial b_1}{\partial s}$$

### Passo de Newton-Raphson

O sistema $2 \times 2$ para os incrementos $(\Delta r, \Delta s)$ é:

$$\begin{bmatrix} c_2 & c_3 \\ c_1 & c_2 \end{bmatrix} \begin{bmatrix} \Delta r \\ \Delta s \end{bmatrix} = \begin{bmatrix} -b_1 \\ -b_0 \end{bmatrix}$$

Solução pela regra de Cramer com $\det = c_2^2 - c_1 c_3$:

$$\Delta r = \frac{b_0 c_3 - b_1 c_2}{\det}, \qquad \Delta s = \frac{b_1 c_1 - b_0 c_2}{\det}$$

### Extração das raízes

Quando $|b_0|$ e $|b_1|$ são suficientemente pequenos, as raízes do fator convergido $x^2 - rx - s = 0$ são:

$$x_{r} = \frac{r \pm \sqrt{r^2 + 4s}}{2}$$

O polinômio é então **deflacionado** (dividido pelo fator convergido) e o processo é repetido até todas as raízes serem encontradas.

---

## 📁 Estrutura do repositório

```
.
├── bairstow.py          # Código principal (método + análises)
├── resultados/          # Gerado automaticamente ao rodar o código
│   ├── analise_1_validacao.png   # Validação grau 7
│   ├── analise_2_apc2.png        # Autovalores do sistema 2-GDL
│   ├── analise_4_fractal.png     # Fractal de Bairstow
│   └── raizes_apc2.csv           # Raízes exportadas
└── README.md
```

---

## 🚀 Como executar

### Pré-requisitos

- Python 3.8 ou superior
- Bibliotecas: `numpy` e `matplotlib`

```bash
pip install numpy matplotlib
```

### Rodando o programa completo

```bash
python bairstow.py
```

O programa executa as quatro análises automaticamente e salva todos os resultados em `resultados/`. A geração do fractal pode levar alguns minutos dependendo da resolução.

---

## 📊 Análises realizadas

### Análise 1 — Validação com polinômio de grau 7

Constrói um polinômio mônico com raízes pré-definidas (3 reais + 2 pares complexos conjugados) e verifica se o método recupera as raízes originais dentro de tolerância numérica.

```
Raízes de entrada:  1, -2, 3,  (1 ± 2i),  (-1 ± i)
Tolerância: 1e-10 | Erro máximo obtido: ~10⁻¹⁵ (precisão de máquina)
```

---

### Análise 2 — Sistema dinâmico 2-GDL (APC2)

Determina os autovalores do polinômio característico da APC2:

$$P(\lambda) = 2\lambda^4 + 5\lambda^3 + 12\lambda^2 + 8\lambda + 8$$

Obtido a partir do sistema com parâmetros:

| Parâmetro | Valor |
|---|---|
| $m_1, m_2$ | 2 kg, 1 kg |
| $k_1, k_2$ | 4 N/m, 2 N/m |
| $c_1, c_2$ | 2 N·s/m, 1 N·s/m |

Os autovalores encontrados são dois pares complexos conjugados, todos com parte real negativa, confirmando a **estabilidade do sistema**.

---

### Análise 3 — Sensibilidade ao chute inicial

Varre uma grade de valores $(r_0, s_0)$ e registra o número de iterações para convergência, identificando regiões problemáticas antes da geração do fractal completo.

---

### Análise 4 — Fractal de Bairstow

Gera um mapa de cor no plano $(r_0, s_0) \in [-3, 3]^2$ onde cada pixel é colorido pelo número de iterações para convergência do primeiro fator quadrático. Regiões pretas correspondem a pontos que não convergiram dentro do limite de iterações.

```
Resolução padrão: 350 × 350 pontos
Máximo de iterações: 80
Colormap: inferno (claro = poucos iters, escuro = muitos iters, preto = divergiu)
```

---

## 🔧 Usando as funções individualmente

```python
from bairstow import bairstow

# Polinômio de grau arbitrário: P(x) = x³ - 6x² + 11x - 6 = (x-1)(x-2)(x-3)
coeffs = [1, -6, 11, -6]
raizes, iters = bairstow(coeffs, r0=0.5, s0=0.5, tol=1e-10)

print(raizes)    # array([3.+0.j, 1.+0.j, 2.+0.j])
print(iters)     # [N_iters_deflacao1, 0]  (0 = fórmula direta para grau 1)
```

```python
# Gerar apenas o fractal para outro polinômio
from bairstow import analise_fractal

# f(x) = x^5 - 1
analise_fractal(
    coeffs=[1, 0, 0, 0, 0, -1],
    r_range=(-2, 2),
    s_range=(-2, 2),
    resolucao=300,
    max_iter=60
)
```

```python
# Exportar raízes para CSV
from bairstow import exportar_raizes

exportar_raizes([2, 5, 12, 8, 8], r0=0.5, s0=0.5, nome="resultados/raizes.csv")
```

---

## 📐 Tabela de funções principais

| Função | Descrição |
|---|---|
| `_calc_b(coeffs, r, s)` | Divisão sintética → coeficientes $b_i$ |
| `_calc_c(b, r, s)` | Segunda divisão sintética → coeficientes $c_i$ |
| `_passo_bairstow(coeffs, r, s)` | Um passo Newton-Raphson → $(\Delta r, \Delta s)$ |
| `bairstow(coeffs, r0, s0, tol, max_iter)` | Método completo → todas as raízes |
| `analise_validacao()` | Validação com polinômio grau 7 |
| `analise_apc2()` | Autovalores do sistema 2-GDL |
| `analise_convergencia()` | Grade de chutes iniciais |
| `analise_fractal(...)` | Mapa fractal no plano $(r_0, s_0)$ |
| `exportar_raizes(...)` | Salva raízes em `.csv` |

---

## 🧩 Convergência e limitações

- O método converge **quadraticamente** (como o Newton-Raphson) próximo à solução.
- A convergência depende do chute inicial $(r_0, s_0)$ — o fractal de Bairstow mapeia visualmente essa sensibilidade.
- Para polinômios com raízes de multiplicidade alta, a convergência pode degradar para linear.
- Se o determinante do sistema $2 \times 2$ for próximo de zero, o código aplica uma pequena perturbação em $(r, s)$ e reinicia a iteração local.
- Validação cruzada automática com `numpy.roots` é exibida na Análise 2.

---

## 📚 Referências

1. Chapra, S. C., Canale, R. P. — *Métodos Numéricos para Engenharia*, McGraw-Hill, 5ª ed. (2008) — Capítulo 9.
2. Bairstow, L. — *Applied Aerodynamics*, Longmans (1920) — Apêndice.
3. Gontijo, R. G. — *Notas de Aula: Cálculo Numérico Aplicado*, UnB (2026).

---

## 👤 Autor

**Mateus Leal Silva** — Matrícula: 221028134
Engenharia Mecânica · Universidade de Brasília
Disciplina: Cálculo Numérico Aplicado
