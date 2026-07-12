# 🔢 Cálculo Numérico Aplicado

> Universidade de Brasília — Faculdade de Tecnologia  
> Disciplina: Cálculo Numérico Aplicado — 2026/1  
> Professor: Prof. Dr. Rafael Gabler Gontijo  
> Autor: Mateus Leal Silva

Repositório didático com implementações computacionais desenvolvidas ao longo da disciplina de **Cálculo Numérico Aplicado**. Cada Programa para Casa (PPC) combina modelagem matemática, implementação do algoritmo, validação numérica, documentação, exercício resolvido e um desafio de extensão.

O objetivo não é apenas apresentar resultados finais, mas manter os métodos como **caixas transparentes**: as etapas numéricas são implementadas explicitamente e comentadas para permitir auditoria, estudo e reprodução.

---

## 🧭 Projetos disponíveis

| Projeto | Aplicação | Métodos principais | Acesso |
|---|---|---|---|
| PPC1 | Sedimentação de uma partícula esférica | Runge–Kutta de 4ª ordem | [Abrir PPC1](PPC1-Runge-Kutta/) |
| PPC2 | Raízes de polinômios e autovalores | Método de Bairstow | [Abrir PPC2](PPC2-Método-de-Bairstow/) |
| PPC3 | Difusão de calor transiente em pastilha de UO₂ | Diferenças finitas, esquema implícito e TDMA | [Abrir PPC3](PPC3-Eliminacao-Gaussiana/) |
| PPC4 | Maximização bidimensional sem restrições | Aclive máximo e Fletcher–Reeves | [Abrir PPC4](PPC4-Otimizacao/) |
| PPC5 | Camada limite laminar de Blasius | Método do Tiro, Newton–Raphson e RK4 | [Abrir PPC5](PPC5-Blasius/) |
| PPC6 | Condução 2D em aleta retangular | Diferenças finitas, Gauss, Liebmann e SOR | [Abrir PPC6](PPC6-Liebmann/) |

---

## 📁 Organização

```text
C-lculo-Num-rico-Aplicado/
├── README.md
├── requirements.txt
├── docs/
│   ├── roteiro-de-estudos.md
│   ├── erros-numericos.md
│   ├── criterios-de-convergencia.md
│   ├── glossario.md
│   └── privacidade-e-direitos-autorais.md
├── PPC1-Runge-Kutta/
│   ├── README.md
│   ├── rk4_sedimentation.py
│   ├── docs/
│   ├── exercicios/
│   ├── inputs/
│   └── resultados/
├── PPC2-Método-de-Bairstow/
├── PPC3-Eliminacao-Gaussiana/
├── PPC4-Otimizacao/
├── PPC5-Blasius/
└── PPC6-Liebmann/
```

Em cada PPC:

- `README.md`: formulação, método, entradas, saídas, execução e validação;
- arquivo `.py`: implementação principal;
- `docs/`: relatórios APC/PPC e orientações de documentação;
- `exercicios/`: desenvolvimento manual de um exemplo autoral;
- `desafio.md`: extensão proposta ao leitor;
- `inputs/`: parâmetros ou dados de entrada, quando aplicável;
- `resultados/`: arquivos e gráficos gerados pela execução.

---

## 🚀 Preparação do ambiente

### 1. Clonar o repositório

```bash
git clone https://github.com/Mateus-Leal-02/C-lculo-Num-rico-Aplicado.git
cd C-lculo-Num-rico-Aplicado
```

### 2. Criar um ambiente virtual

No Windows:

```powershell
python -m venv .venv
.venv\Scripts\activate
```

No Linux ou WSL:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Instalar as dependências

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

---

## ▶️ Execução rápida

Entre na pasta do projeto desejado e execute o código principal. Exemplo:

```bash
cd PPC5-Blasius
python Blasius.py
```

Para consultar os argumentos disponíveis em códigos com interface de linha de comando:

```bash
python NomeDoCodigo.py --help
```

As instruções específicas, parâmetros padrão e arquivos de saída estão documentados no `README.md` de cada PPC.

---

## ✅ Estratégia de validação

As implementações são verificadas por uma ou mais das seguintes abordagens:

1. comparação com solução analítica conhecida;
2. comparação com ponto ótimo, raízes ou coeficientes de referência;
3. estudo de refinamento de passo ou malha;
4. análise do resíduo e do critério de convergência;
5. comparação entre métodos independentes;
6. verificação de comportamento físico esperado.

Consulte [`docs/criterios-de-convergencia.md`](docs/criterios-de-convergencia.md) e [`docs/erros-numericos.md`](docs/erros-numericos.md).

---

## 🧩 Exercícios e desafios

Cada pasta contém:

- um **exercício resolvido**, destinado a reproduzir manualmente uma etapa central do método;
- um **desafio para o leitor**, destinado a ampliar o código, comparar algoritmos ou avaliar sensibilidade numérica.

Esses materiais são autorais e complementares aos relatórios da disciplina.

---
