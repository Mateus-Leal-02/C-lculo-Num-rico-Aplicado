# Cálculo Numérico Aplicado

## Como usar no Linux

```bash
git clone LINK_DO_SEU_REPOSITORIO
cd NOME_DO_REPOSITORIO
ls
python3 Resumo.py
python3 -i Resumo.py
```

Para achar uma função rapidamente:

```bash
grep -R "def thomas" -n .
grep -R "def bissecao" -n .
grep -R "def gauss_seidel" -n .
grep -R "def fletcher" -n .
```

Para copiar e editar:

```bash
cp Resumo.py prova.py
nano prova.py
python3 prova.py
```

## Como decidir o método

### Raiz de função
- Tem intervalo `[a,b]` com troca de sinal: **bisseção** ou **falsa posição**.
- Tem derivada ou é fácil derivar: **Newton-Raphson**.
- Não tem derivada, mas tem dois chutes: **secante**.
- É polinômio de grau alto: **Bairstow** ou `np.roots`, se permitido.

### Sistema linear
- Matriz geral pequena/média: **eliminação de Gauss com pivotamento**.
- Mesmo A para vários b: **LU**.
- Sistema tridiagonal: **Thomas/TDMA**.
- Matriz simétrica definida positiva: **Cholesky**.
- Matriz grande/esparsa e diagonal dominante: **Gauss-Seidel**.

### Otimização
- Função 1D sem restrição: **razão áurea** ou **interpolação quadrática**.
- Função 2D/multidimensional com gradiente: **aclive máximo**.
- Função quadrática com gradiente: **Fletcher-Reeves** costuma convergir rápido.

## Estrutura mental para qualquer código da prova

1. Definir funções matemáticas: `f`, `df`, `grad`, `A`, `b`, etc.
2. Definir entrada/chute inicial/tolerância/número máximo de iterações.
3. Rodar método iterativo.
4. Imprimir tabela/log: iteração, aproximação, erro, valor da função/resíduo.
5. Validar resultado: substituir na equação, calcular resíduo ou comparar com resposta esperada.

## Critérios de parada que quase sempre funcionam

- Raiz: `abs(f(x)) < tol` ou erro relativo aproximado pequeno.
- Sistema linear: `norma(A@x - b) < tol` ou erro relativo entre iterações.
- Otimização: `norma(grad(x)) < tol`.
- EDO: número de passos `n = int((tf-t0)/h)`.

## Fórmulas-chave

### Bisseção
`xr = (a+b)/2`

### Falsa posição
`xr = b - f(b)*(a-b)/(f(a)-f(b))`

### Newton-Raphson
`x_{k+1} = x_k - f(x_k)/f'(x_k)`

### Secante
`x_{k+1} = x_k - f(x_k)*(x_{k-1}-x_k)/(f(x_{k-1})-f(x_k))`

### Gauss-Seidel
`x_i = (b_i - soma(a_ij*x_j, j!=i))/a_ii`

Com relaxação:
`x_i_novo = lambda*x_i_GS + (1-lambda)*x_i_antigo`

### Thomas/TDMA
Sistema: `a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i`.
Usar vetores `a`, `b`, `c`, `d`, com `a[0]=0` e `c[-1]=0`.

### Otimização por gradiente
Atualização geral:
`z_{k+1} = z_k + h_k p_k`

Aclive máximo:
`p_k = grad(f)(z_k)`

Fletcher-Reeves:
`beta = ||g_{k+1}||^2 / ||g_k||^2`
`p_{k+1} = g_{k+1} + beta*p_k`

### RK4
`k1=f(t,y)`
`k2=f(t+h/2, y+h*k1/2)`
`k3=f(t+h/2, y+h*k2/2)`
`k4=f(t+h, y+h*k3)`
`y_{n+1}=y_n + h*(k1+2*k2+2*k3+k4)/6`
