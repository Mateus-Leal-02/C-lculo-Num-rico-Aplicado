# 🌊 Sedimentação de Esfera — Runge-Kutta de 4ª Ordem

> **Programa para Casa #1 · Cálculo Numérico Aplicado**  
> Prof. Dr. Rafael Gabler Gontijo · UnB · março de 2026

---

## 📌 Sobre o projeto

Este repositório contém a implementação do **método de Runge-Kutta de 4ª ordem (RK4)** aplicado à equação de movimento de uma esfera sedimentando em um fluido viscoso em regime de baixo Reynolds.

O problema físico consiste em uma esfera de raio $a$ e densidade $\rho_s$ se deslocando em um fluido de viscosidade $\eta$ sob ação da gravidade. Dependendo do regime de escoamento, diferentes modelos de arrasto são considerados:

| Regime | Força de Arrasto | Solução |
|---|---|---|
| $Re \to 0$ (Stokes) | $f_d = -6\pi\eta a v$ | Analítica: $v^* = 1 - e^{-t^*/St}$ |
| $Re_s \lesssim 1$ (Oseen) | Stokes + termo quadrático | Analítica (Ricatti) + RK4 |
| $Re \gg 1$ (turbulento) | Quadrático não-linear | Somente numérica |

---

## 🧮 Equação adimensional do movimento

Após adimensionalização com a velocidade de Stokes $U_s = 2a^2 \Delta\rho\, g / (9\eta)$ e o tempo convectivo $a/U_s$, a equação do movimento assume a forma:

$$
\frac{dv^*}{dt^*}
= 1 - v^* - \frac{3}{8}\,\mathrm{Re}_s\,(v^*)^2,
\qquad
v^*(0)=0
$$

onde:
- $St = \dfrac{2\rho_s a U_s}{9\eta}$ — **número de Stokes**: inércia da partícula vs. arrasto viscoso
- $Re_s = \dfrac{\rho_f U_s a}{\eta}$ — **Reynolds de partícula**: importância do termo quadrático de Oseen

Para $Re_s = 0$, recupera-se a equação linear com solução exata $v^* = 1 - e^{-t^*/St}$.

---

## ⚙️ Método numérico — RK4

O avanço temporal é feito pelo clássico Runge-Kutta de 4ª ordem:

$$k_1 = f(t_i,\; v_i)$$
$$k_2 = f\!\left(t_i + \tfrac{h}{2},\; v_i + \tfrac{h}{2}k_1\right)$$
$$k_3 = f\!\left(t_i + \tfrac{h}{2},\; v_i + \tfrac{h}{2}k_2\right)$$
$$k_4 = f\!\left(t_i + h,\; v_i + h\,k_3\right)$$
$$v_{i+1} = v_i + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

O erro de truncamento local é $\mathcal{O}(h^5)$, e o erro global é $\mathcal{O}(h^4)$.  
**Nenhuma biblioteca de integração numérica externa é utilizada** — o RK4 é implementado do zero.

---

## 📁 Estrutura do repositório

```
.
├── rk4_sedimentation.py   # Código principal
├── resultados/            # Gerado automaticamente ao rodar o código
│   ├── analise_1_stokes.png
│   ├── analise_2_convergencia.png
│   ├── analise_3_reynolds.png
│   └── dados_St1.0_Re0.0_h0.05.csv
└── README.md
```

---

## 🚀 Como executar

### Pré-requisitos

- Python 3.8 ou superior
- Bibliotecas: `numpy` e `matplotlib`

Instale as dependências com:

```bash
pip install numpy matplotlib
```

### Rodando o programa

```bash
python rk4_sedimentation.py
```

O programa executa automaticamente as três análises e salva os resultados na pasta `resultados/`.

---

## 📊 Análises realizadas

### Análise 1 — Validação com a solução analítica ($Re_s = 0$)

Compara a solução RK4 com a solução exata $v^* = 1 - e^{-t^*/St}$ para diferentes valores do número de Stokes.

```
Parâmetros: Re_s = 0, h = 0.05, t_max = 15
St testados: [0.5, 1.0, 2.0, 4.0]
```

> Quanto maior o $St$, mais lenta é a relaxação para a velocidade terminal.

---

### Análise 2 — Convergência em $h$

Reduz progressivamente o passo de tempo e mede o erro máximo em relação à solução analítica. Espera-se decaimento proporcional a $h^4$.

```
Parâmetros: St = 1.0, Re_s = 0, t_max = 10
h testados: [1.0, 0.5, 0.1, 0.05, 0.01]
```

> O gráfico log-log mostra a linha de referência $\mathcal{O}(h^4)$ confirmando a ordem do método.

---

### Análise 3 — Efeito do número de Reynolds

Resolve o PVI com correção de Oseen para diferentes $Re_s$ e compara com o limite assintótico de Stokes.

```
Parâmetros: St = 1.0, h = 0.05, t_max = 15
Re_s testados: [0.0, 0.2, 0.5, 1.0]
```

> O aumento de $Re_s$ reduz a velocidade terminal (o arrasto quadrático é maior), gerando desvio crescente em relação ao regime de Stokes.

---

## 🔧 Usando funções individualmente

As funções podem ser importadas e chamadas separadamente:

```python
from rk4_sedimentation import rk4, f, solucao_analitica_stokes
import numpy as np

# Resolver para St = 2, Re_s = 0.5, h = 0.01
t, v = rk4(f, v0=0.0, t0=0.0, t_max=20.0, h=0.01, St=2.0, Re_s=0.5)

# Solução analítica de Stokes para comparação
v_ex = solucao_analitica_stokes(t, St=2.0)

print(f"Velocidade terminal numérica:  {v[-1]:.6f}")
print(f"Velocidade terminal analítica: {v_ex[-1]:.6f}")
```

---

## 📐 Definições dos parâmetros adimensionais

| Símbolo | Expressão | Significado físico |
|---|---|---|
| $U_s$ | $2a^2 \Delta\rho\, g / (9\eta)$ | Velocidade terminal de Stokes |
| $v^*$ | $v_z / U_s$ | Velocidade adimensional |
| $t^*$ | $t \cdot U_s / a$ | Tempo adimensional |
| $St$ | $2\rho_s a U_s / (9\eta)$ | Número de Stokes |
| $Re_s$ | $\rho_f U_s a / \eta$ | Reynolds de partícula |

---

## 📚 Referências

1. Sobral, Y. D., Oliveira, T. F., Cunha, F. R. — *On the unsteady forces during the motion of a sedimenting particle*, Powder Technology, 178 (2007), 129–141.
2. Chapra, S. C., Canale, R. P. — *Métodos Numéricos para Engenharia*, McGraw-Hill, 5ª ed. (2008).

---

## 👤 Autor

**Mateus Leal** — Matrícula: 221028134  
Engenharia Mecânica · Universidade de Brasília  
Disciplina: Cálculo Numérico Aplicado

