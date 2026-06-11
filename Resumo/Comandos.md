# COMANDOS LINUX — GITHUB, ARQUIVOS E PYTHON

Cola rápida para usar no terminal Linux/WSL durante a prova.

Foco:
- Entrar no repositório do GitHub
- Navegar por pastas
- Criar, abrir, editar, renomear e excluir arquivos
- Rodar códigos Python
- Salvar alterações no GitHub

---

## 1. Começar no terminal

### Ir para sua pasta principal

```bash
cd ~
```

### Ver onde você está

```bash
pwd
```

### Listar arquivos e pastas

```bash
ls
```

### Listar com detalhes

```bash
ls -la
```

### Limpar a tela

```bash
clear
```

---

## 2. Clonar o repositório do GitHub

Use o link HTTPS do repositório.

Formato:

```bash
git clone https://github.com/USUARIO/NOME_DO_REPOSITORIO.git
```

Exemplo:

```bash
git clone https://github.com/Mateus-Leal-02/C-lculo-Num-rico-Aplicado.git
```

Depois entre na pasta:

```bash
cd C-lculo-Num-rico-Aplicado
```

Confira os arquivos:

```bash
ls
```

---

## 3. Entrar em um repositório já clonado

```bash
cd ~
ls
cd NOME_DO_REPOSITORIO
```

Atualizar com o GitHub:

```bash
git pull
```

---

## 4. Navegar entre pastas

### Entrar em uma pasta

```bash
cd nome_da_pasta
```

### Voltar uma pasta

```bash
cd ..
```

### Voltar para a pasta principal

```bash
cd ~
```

### Ver caminho atual

```bash
pwd
```

### Ver arquivos da pasta atual

```bash
ls
```

### Autocompletar

Digite o começo do nome e aperte `TAB`.

```bash
cd Cal<TAB>
```

---

## 5. Criar arquivos e pastas

### Criar pasta

```bash
mkdir nome_da_pasta
```

Exemplo:

```bash
mkdir prova
```

### Criar várias pastas

```bash
mkdir PPC3 PPC4 Prova
```

### Criar pasta dentro de pasta

```bash
mkdir -p prova/testes
```

### Criar arquivo vazio

```bash
touch prova.py
```

### Criar e abrir arquivo Python

```bash
nano prova.py
```

---

## 6. Abrir e editar com nano

### Abrir arquivo

```bash
nano prova.py
```

### Criar arquivo novo

```bash
nano novo_arquivo.py
```

### Salvar e sair

Dentro do nano:

```text
Ctrl + O
Enter
Ctrl + X
```

Significado:

```text
Ctrl + O  -> salvar
Enter     -> confirmar nome
Ctrl + X  -> sair
```

### Procurar palavra dentro do nano

```text
Ctrl + W
```

---

## 7. Processo interativo para programar

Ciclo principal:

```bash
nano prova.py
python3 prova.py
nano prova.py
python3 prova.py
```

Ou seja:

```text
editar -> salvar -> rodar -> ver erro/resultado -> editar de novo -> rodar de novo
```

---

## 8. Rodar código Python

### Rodar arquivo

```bash
python3 prova.py
```

### Rodar exemplos do kit

```bash
python3 exemplos_kit_prova.py
```

### Ver versão do Python

```bash
python3 --version
```

### Testar NumPy

```bash
python3 -c "import numpy as np; print(np.__version__)"
```

---

## 9. Criar um código Python simples

```bash
nano teste.py
```

Digite:

```python
print("Funcionou!")

a = 2
b = 3
print("soma =", a + b)
```

Salve:

```text
Ctrl + O
Enter
Ctrl + X
```

Rode:

```bash
python3 teste.py
```

---

## 10. Criar arquivo rapidamente sem nano

```bash
cat > teste.py << 'PY'
print("Funcionou!")
PY
```

Rodar:

```bash
python3 teste.py
```

Exemplo com NumPy:

```bash
cat > teste_numpy.py << 'PY'
import numpy as np

A = np.array([[2, 1], [1, 3]], dtype=float)
b = np.array([1, 2], dtype=float)

x = np.linalg.solve(A, b)

print("x =", x)
print("A @ x =", A @ x)
PY
```

Rodar:

```bash
python3 teste_numpy.py
```

---

## 11. Ver conteúdo de arquivos

### Mostrar arquivo inteiro

```bash
cat prova.py
```

### Ler página por página

```bash
less prova.py
```

Dentro do `less`:

```text
q       -> sair
/termo  -> procurar
n       -> próxima ocorrência
```

### Ver início

```bash
head prova.py
```

### Ver final

```bash
tail prova.py
```

### Ver primeiras 40 linhas

```bash
head -n 40 prova.py
```

### Ver últimas 40 linhas

```bash
tail -n 40 prova.py
```

---

## 12. Copiar arquivos

### Copiar arquivo

```bash
cp origem.py copia.py
```

Exemplo:

```bash
cp exemplos_kit_prova.py prova.py
```

### Copiar o kit para editar sem alterar o original

```bash
cp kit_prova_calculo_numerico.py prova.py
```

### Copiar arquivo para pasta

```bash
cp prova.py Prova/
```

### Copiar pasta inteira

```bash
cp -r pasta_original pasta_copia
```

---

## 13. Renomear e mover arquivos

### Renomear arquivo

```bash
mv nome_antigo.py nome_novo.py
```

Exemplo:

```bash
mv teste.py prova.py
```

### Mover arquivo para pasta

```bash
mv prova.py Prova/
```

### Mover e renomear ao mesmo tempo

```bash
mv prova.py Prova/prova_final.py
```

---

## 14. Excluir arquivos e pastas

Cuidado com esses comandos.

### Excluir arquivo

```bash
rm arquivo.py
```

### Excluir vários arquivos

```bash
rm arquivo1.py arquivo2.py
```

### Excluir pasta vazia

```bash
rmdir nome_da_pasta
```

### Excluir pasta com arquivos dentro

```bash
rm -r nome_da_pasta
```

Evite `rm -r` se não tiver certeza.

---

## 15. Procurar arquivos

### Procurar todos os arquivos Python

```bash
find . -name "*.py"
```

### Procurar arquivos com "prova" no nome

```bash
find . -name "*prova*"
```

### Procurar o kit

```bash
find . -name "kit_prova_calculo_numerico.py"
```

### Procurar a partir da pasta principal

```bash
find ~ -name "kit_prova_calculo_numerico.py"
```

---

## 16. Procurar funções dentro dos arquivos

### Procurar função específica

```bash
grep -R "def bissecao" -n .
```

### Procurar palavra

```bash
grep -R "thomas" -n .
```

### Procurar ignorando maiúsculas/minúsculas

```bash
grep -Ri "gauss" .
```

### Principais buscas do kit

```bash
grep -R "def bissecao" -n .
grep -R "def falsa_posicao" -n .
grep -R "def newton_raphson" -n .
grep -R "def secante" -n .
grep -R "def newton_sistema" -n .
grep -R "def bairstow" -n .
grep -R "def eliminacao" -n .
grep -R "def decomposicao_lu" -n .
grep -R "def thomas" -n .
grep -R "def cholesky" -n .
grep -R "def gauss_seidel" -n .
grep -R "def razao_aurea" -n .
grep -R "def interpolacao_quadratica" -n .
grep -R "def aclive_maximo" -n .
grep -R "def fletcher" -n .
grep -R "def euler" -n .
grep -R "def rk4" -n .
```

---

## 17. Salvar saída do programa

### Salvar resultado em arquivo

```bash
python3 prova.py > resultado.txt
```

### Ver resultado salvo

```bash
cat resultado.txt
```

### Salvar resultado e erros

```bash
python3 prova.py > resultado.txt 2>&1
```

### Ver resultado página por página

```bash
less resultado.txt
```

---

## 18. Git: ver, salvar e enviar alterações

### Ver status

```bash
git status
```

### Adicionar todos os arquivos modificados

```bash
git add .
```

### Criar commit

```bash
git commit -m "Adiciona arquivos de prova"
```

### Enviar para o GitHub

```bash
git push
```

### Fluxo completo

```bash
git status
git add .
git commit -m "Atualiza arquivos de calculo numerico"
git push
```

---

## 19. Baixar atualizações do GitHub

```bash
git pull
```

Ou:

```bash
cd ~/NOME_DO_REPOSITORIO
git pull
```

---

## 20. Criar solução usando o kit

Dentro da pasta do repositório:

```bash
nano prova.py
```

Coloque:

```python
import numpy as np
from kit_prova_calculo_numerico import *

print("Arquivo da prova rodando")
```

Salve:

```text
Ctrl + O
Enter
Ctrl + X
```

Rode:

```bash
python3 prova.py
```

---

## 21. Exemplo usando função do kit

```bash
nano prova.py
```

Cole:

```python
from kit_prova_calculo_numerico import bissecao

def f(x):
    return x**3 - x - 2

raiz = bissecao(f, 1, 2)

print("raiz =", raiz)
print("f(raiz) =", f(raiz))
```

Salve e rode:

```bash
python3 prova.py
```

---

## 22. Corrigir erros comuns

### Arquivo não encontrado

Mensagem:

```text
No such file or directory
```

Faça:

```bash
pwd
ls
```

Você provavelmente está na pasta errada.

### Python não encontrado

Se isso der erro:

```bash
python prova.py
```

Use:

```bash
python3 prova.py
```

### Módulo do kit não encontrado

Mensagem:

```text
ModuleNotFoundError: No module named 'kit_prova_calculo_numerico'
```

Confira se o kit está na mesma pasta:

```bash
ls
```

Se não estiver, procure:

```bash
find ~ -name "kit_prova_calculo_numerico.py"
```

### Erro de indentação

Mensagem:

```text
IndentationError
```

Exemplo certo:

```python
for i in range(5):
    print(i)
```

Exemplo errado:

```python
for i in range(5):
print(i)
```

### Programa travou

Interrompa:

```text
Ctrl + C
```

---

## 23. Atalhos úteis

```text
TAB              -> autocompletar
seta para cima   -> comando anterior
seta para baixo  -> próximo comando
Ctrl + C         -> parar programa travado
Ctrl + L         -> limpar tela
Ctrl + A         -> início da linha
Ctrl + E         -> fim da linha
Ctrl + U         -> apagar antes do cursor
Ctrl + K         -> apagar depois do cursor
```

---

## 24. Checklist da prova

```text
1. Entre no repositório.
2. Rode ls para ver os arquivos.
3. Copie um exemplo parecido para prova.py.
4. Edite com nano.
5. Rode com python3 prova.py.
6. Corrija erros.
7. Salve saída se precisar.
```

Comandos:

```bash
cd ~
cd NOME_DO_REPOSITORIO
git pull
ls
cp exemplos_kit_prova.py prova.py
nano prova.py
python3 prova.py
```

---

## 25. Sequência mais importante para decorar

```bash
cd ~
cd NOME_DO_REPOSITORIO
ls
cp exemplos_kit_prova.py prova.py
nano prova.py
python3 prova.py
nano prova.py
python3 prova.py
```

Salvar no nano:

```text
Ctrl + O
Enter
Ctrl + X
```
