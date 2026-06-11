# GRÁFICOS EM PYTHON — COLA PARA A PROVA

Este arquivo serve para quando o enunciado pedir algo como:

```text
- plote o gráfico da função
- gere a curva de convergência
- mostre o caminho das iterações
- faça gráfico de x por y
- gere curvas de nível
- salve a figura
```

A regra de segurança para prova em Linux/WSL é:

```text
Use plt.savefig("nome.png") em vez de depender de plt.show().
```

Assim o gráfico é gerado como arquivo de imagem.

---

## 1. Testar se Matplotlib está instalado

```bash
python3 -c "import matplotlib; print(matplotlib.__version__)"
```

Se der erro, tente instalar:

```bash
sudo apt update
sudo apt install python3-matplotlib -y
```

Ou, se for permitido:

```bash
python3 -m pip install matplotlib
```

---

## 2. Estrutura básica de qualquer gráfico

```python
import numpy as np
import matplotlib.pyplot as plt

# 1. Criar os dados
x = np.linspace(0, 10, 100)
y = x**2

# 2. Criar o gráfico
plt.figure()
plt.plot(x, y)

# 3. Dar nomes
plt.xlabel("x")
plt.ylabel("y")
plt.title("Gráfico de y = x²")
plt.grid(True)

# 4. Salvar
plt.savefig("grafico.png", dpi=300, bbox_inches="tight")

# 5. Opcional
plt.show()
```

Rodar:

```bash
python3 arquivo.py
```

Ver se a figura apareceu:

```bash
ls
```

---

## 3. Exemplo completo: gráfico de uma função

Crie:

```bash
nano grafico_funcao.py
```

Cole:

```python
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x**3 - x - 2

x = np.linspace(-3, 3, 400)
y = f(x)

plt.figure()
plt.plot(x, y, label="f(x) = x³ - x - 2")
plt.axhline(0, linestyle="--", linewidth=1)
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Gráfico da função")
plt.grid(True)
plt.legend()
plt.savefig("grafico_funcao.png", dpi=300, bbox_inches="tight")
plt.show()
```

Salve:

```text
Ctrl + O
Enter
Ctrl + X
```

Rode:

```bash
python3 grafico_funcao.py
```

Resultado esperado:

```text
grafico_funcao.png
```

---

## 4. Gráfico de pontos

Use quando tiver dados tabulados.

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.array([0, 1, 2, 3, 4], dtype=float)
y = np.array([1, 2, 4, 8, 16], dtype=float)

plt.figure()
plt.plot(x, y, "o", label="dados")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Gráfico de pontos")
plt.grid(True)
plt.legend()
plt.savefig("grafico_pontos.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## 5. Gráfico de linha com pontos

Bom para mostrar convergência, temperatura ao longo do tempo, erro por iteração etc.

```python
import numpy as np
import matplotlib.pyplot as plt

iteracao = np.array([0, 1, 2, 3, 4, 5])
erro = np.array([100, 40, 15, 5, 1, 0.1])

plt.figure()
plt.plot(iteracao, erro, marker="o")
plt.xlabel("Iteração")
plt.ylabel("Erro")
plt.title("Convergência do método")
plt.grid(True)
plt.savefig("convergencia.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## 6. Gráfico com escala logarítmica

Útil quando o erro cai muito rápido.

```python
import numpy as np
import matplotlib.pyplot as plt

iteracao = np.array([0, 1, 2, 3, 4, 5])
erro = np.array([100, 40, 15, 5, 1, 0.1])

plt.figure()
plt.semilogy(iteracao, erro, marker="o")
plt.xlabel("Iteração")
plt.ylabel("Erro")
plt.title("Convergência em escala log")
plt.grid(True)
plt.savefig("convergencia_log.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## 7. Gráfico com duas curvas

Use quando comparar método numérico com solução analítica, Euler com RK4, etc.

```python
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, 5, 100)
y1 = np.exp(-t)
y2 = 1 / (1 + t)

plt.figure()
plt.plot(t, y1, label="exp(-t)")
plt.plot(t, y2, label="1/(1+t)")
plt.xlabel("t")
plt.ylabel("y")
plt.title("Comparação de curvas")
plt.grid(True)
plt.legend()
plt.savefig("comparacao.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## 8. Gráfico de barras

Use para comparar valores discretos.

```python
import matplotlib.pyplot as plt

metodos = ["Bisseção", "Newton", "Secante"]
iteracoes = [25, 6, 8]

plt.figure()
plt.bar(metodos, iteracoes)
plt.xlabel("Método")
plt.ylabel("Número de iterações")
plt.title("Comparação entre métodos")
plt.grid(axis="y")
plt.savefig("barras.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## 9. Curvas de nível de função f(x,y)

Muito útil para otimização multidimensional.

Exemplo:

```text
f(x,y) = 2xy + 2x - x² - 2y²
```

```python
import numpy as np
import matplotlib.pyplot as plt

def f(x, y):
    return 2*x*y + 2*x - x**2 - 2*y**2

x = np.linspace(-4, 4, 200)
y = np.linspace(-2, 5, 200)

X, Y = np.meshgrid(x, y)
Z = f(X, Y)

plt.figure()
cont = plt.contour(X, Y, Z, levels=20)
plt.clabel(cont, inline=True, fontsize=8)

plt.xlabel("x")
plt.ylabel("y")
plt.title("Curvas de nível de f(x,y)")
plt.grid(True)
plt.savefig("curvas_nivel.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## 10. Curvas de nível com caminho das iterações

Use quando o método gera pontos:

```text
z0, z1, z2, z3, ...
```

Exemplo para otimização:

```python
import numpy as np
import matplotlib.pyplot as plt

def f(x, y):
    return 2*x*y + 2*x - x**2 - 2*y**2

# Exemplo de caminho. Troque pelos pontos reais do seu método.
pontos = np.array([
    [-2.0, 3.0],
    [1.0, 0.0],
    [1.7, 0.8],
    [2.0, 1.0]
], dtype=float)

x = np.linspace(-4, 4, 200)
y = np.linspace(-2, 5, 200)
X, Y = np.meshgrid(x, y)
Z = f(X, Y)

plt.figure()
cont = plt.contour(X, Y, Z, levels=20)
plt.clabel(cont, inline=True, fontsize=8)

plt.plot(pontos[:, 0], pontos[:, 1], marker="o", label="iterações")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Caminho das iterações")
plt.grid(True)
plt.legend()
plt.savefig("caminho_iteracoes.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## 11. Gráfico de superfície 3D

Use somente se o professor pedir superfície.

```python
import numpy as np
import matplotlib.pyplot as plt

def f(x, y):
    return 2*x*y + 2*x - x**2 - 2*y**2

x = np.linspace(-4, 4, 100)
y = np.linspace(-2, 5, 100)

X, Y = np.meshgrid(x, y)
Z = f(X, Y)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(X, Y, Z)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("f(x,y)")
ax.set_title("Superfície da função")

plt.savefig("superficie_3d.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## 12. Gráfico para EDO/RK4

Quando você tiver vetores `ts` e `ys`.

```python
import numpy as np
import matplotlib.pyplot as plt

# Exemplo fictício
ts = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ys = np.array([0.5, 0.829, 1.214, 1.649, 2.127, 2.641])

plt.figure()
plt.plot(ts, ys, marker="o")
plt.xlabel("t")
plt.ylabel("y(t)")
plt.title("Solução numérica da EDO")
plt.grid(True)
plt.savefig("edo_rk4.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## 13. Gráfico para sistema de EDOs

Exemplo: estado `ys` tem duas colunas, posição e velocidade.

```python
import numpy as np
import matplotlib.pyplot as plt

ts = np.linspace(0, 10, 100)
x = np.cos(ts)
v = -np.sin(ts)

plt.figure()
plt.plot(ts, x, label="x(t)")
plt.plot(ts, v, label="v(t)")
plt.xlabel("t")
plt.ylabel("estado")
plt.title("Sistema de EDOs")
plt.grid(True)
plt.legend()
plt.savefig("sistema_edo.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## 14. Salvar dados em CSV e gráfico em PNG

Útil se o enunciado pedir arquivo de saída.

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 10, 50)
y = np.sin(x)

dados = np.column_stack((x, y))
np.savetxt("dados.csv", dados, delimiter=",", header="x,y", comments="")

plt.figure()
plt.plot(x, y)
plt.xlabel("x")
plt.ylabel("sen(x)")
plt.title("Seno")
plt.grid(True)
plt.savefig("seno.png", dpi=300, bbox_inches="tight")
plt.show()
```

Arquivos gerados:

```text
dados.csv
seno.png
```

---

## 15. Ver arquivos gerados

```bash
ls
```

Procurar imagens:

```bash
find . -name "*.png"
```

---

## 16. Abrir imagem pelo terminal

Em alguns ambientes funciona:

```bash
xdg-open grafico.png
```

Se não abrir, pelo menos confirme que o arquivo foi criado:

```bash
ls -lh grafico.png
```

---

## 17. Erros comuns

### Erro: No module named matplotlib

Instale:

```bash
sudo apt update
sudo apt install python3-matplotlib -y
```

### O programa roda, mas não aparece janela

Normal no WSL/terminal. Por isso use:

```python
plt.savefig("grafico.png", dpi=300, bbox_inches="tight")
```

### O arquivo PNG não apareceu

Confira se `plt.savefig(...)` veio antes de `plt.show()`.

Melhor ordem:

```python
plt.savefig("grafico.png", dpi=300, bbox_inches="tight")
plt.show()
```

### O gráfico ficou cortado

Use:

```python
plt.savefig("grafico.png", dpi=300, bbox_inches="tight")
```

### O gráfico anterior misturou com o novo

Use sempre antes de cada gráfico:

```python
plt.figure()
```

---

## 18. Modelo pronto para colar na prova

```python
import numpy as np
import matplotlib.pyplot as plt

# Defina os dados
x = np.linspace(0, 10, 100)
y = x**2

# Faça o gráfico
plt.figure()
plt.plot(x, y, label="resultado")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Resultado numérico")
plt.grid(True)
plt.legend()

# Salve
plt.savefig("resultado.png", dpi=300, bbox_inches="tight")
plt.show()
```

Rodar:

```bash
python3 prova.py
```

Ver figura:

```bash
ls
```
