"""
EXEMPLOS - Cálculo Numérico Aplicado

Como usar:
1) Deixe este arquivo na mesma pasta do Exemplos.py
2) Rode:
       python3 Exemplos.py

A ideia é copiar o BLOCO do método que encaixa no enunciado
e adaptar f, df, grad, A, b, intervalo, chute inicial etc.
"""

import numpy as np

from kit_prova_calculo_numerico import (
    norma,
    erro_relativo_atual,
    bissecao,
    falsa_posicao,
    newton_raphson,
    secante,
    newton_sistema,
    bairstow,
    eliminacao_gauss_pivot,
    decomposicao_lu,
    substituicao_progressiva,
    substituicao_regressiva,
    resolver_lu,
    thomas_algorithm,
    cholesky_solve,
    gauss_seidel,
    razao_aurea,
    interpolacao_quadratica_1d,
    aclive_maximo,
    fletcher_reeves_max,
    euler,
    rk4,
)


def titulo(txt):
    print("\n" + "=" * 70)
    print(txt)
    print("=" * 70)


# ============================================================
# 0) FUNÇÕES ÚTEIS
# ============================================================

titulo("0.1) norma(v)")

v = np.array([3, 4])
print("v =", v)
print("norma(v) =", norma(v))
# Esperado: 5


titulo("0.2) erro_relativo_atual(novo, antigo)")

novo = np.array([1.10, 2.20, 3.30])
antigo = np.array([1.00, 2.00, 3.00])
print("novo =", novo)
print("antigo =", antigo)
print("erro relativo aproximado percentual =", erro_relativo_atual(novo, antigo), "%")


# ============================================================
# 1) RAÍZES DE FUNÇÕES
# Problema exemplo comum:
# Encontrar a raiz de f(x) = x^3 - x - 2
# A raiz real está entre 1 e 2.
# ============================================================

titulo("1.1) bissecao(f, a, b)")

f = lambda x: x**3 - x - 2
raiz = bissecao(f, a=1, b=2, tol=1e-6, verbose=False)
print("raiz por bisseção =", raiz)
print("f(raiz) =", f(raiz))


titulo("1.2) falsa_posicao(f, a, b)")

f = lambda x: x**3 - x - 2
raiz = falsa_posicao(f, a=1, b=2, tol=1e-6, verbose=False)
print("raiz por falsa posição =", raiz)
print("f(raiz) =", f(raiz))


titulo("1.3) newton_raphson(f, df, x0)")

f = lambda x: x**3 - x - 2
df = lambda x: 3*x**2 - 1
raiz = newton_raphson(f, df, x0=1.5, tol=1e-6, verbose=False)
print("raiz por Newton-Raphson =", raiz)
print("f(raiz) =", f(raiz))


titulo("1.4) secante(f, x0, x1)")

f = lambda x: x**3 - x - 2
raiz = secante(f, x0=1, x1=2, tol=1e-6, verbose=False)
print("raiz por secante =", raiz)
print("f(raiz) =", f(raiz))


titulo("1.5) newton_sistema(F, J, x0)")

# Sistema não linear:
# F1(x,y) = x^2 + y^2 - 4 = 0       -> círculo de raio 2
# F2(x,y) = x - y = 0               -> reta x = y
# Solução positiva esperada: x = y = sqrt(2)

def F(z):
    x, y = z
    return np.array([
        x**2 + y**2 - 4,
        x - y
    ], dtype=float)

def J(z):
    x, y = z
    return np.array([
        [2*x, 2*y],
        [1, -1]
    ], dtype=float)

sol = newton_sistema(F, J, x0=[1.0, 1.5], tol=1e-8, verbose=False)
print("solução do sistema =", sol)
print("F(sol) =", F(sol))


# ============================================================
# 2) RAÍZES DE POLINÔMIOS: BAIRSTOW
# ============================================================

titulo("2.1) bairstow(coef, r, s)")

# Polinômio: p(x) = x^3 - 6x^2 + 11x - 6
# Raízes exatas: 1, 2, 3
coef = [1, -6, 11, -6]
raizes = bairstow(coef, r=0.5, s=-0.5, tol=1e-6, verbose=False)
print("raízes por Bairstow =", raizes)


# ============================================================
# 3) SISTEMAS LINEARES
# Problema exemplo:
# 3x - 0.1y - 0.2z = 7.85
# 0.1x + 7y - 0.3z = -19.3
# 0.3x - 0.2y + 10z = 71.4
# Solução esperada: [3, -2.5, 7]
# ============================================================

A = np.array([
    [3, -0.1, -0.2],
    [0.1, 7, -0.3],
    [0.3, -0.2, 10]
], dtype=float)

b = np.array([7.85, -19.3, 71.4], dtype=float)


titulo("3.1) eliminacao_gauss_pivot(A, b)")

x = eliminacao_gauss_pivot(A, b)
print("x =", x)
print("A @ x =", A @ x)


titulo("3.2) decomposicao_lu(A)")

L, U = decomposicao_lu(A)
print("L =")
print(L)
print("U =")
print(U)
print("verificação L @ U =")
print(L @ U)


titulo("3.3) substituicao_progressiva(L, b)")

# Exemplo triangular inferior:
L_ex = np.array([
    [1, 0, 0],
    [2, 1, 0],
    [3, -1, 1]
], dtype=float)
b_ex = np.array([1, 4, 2], dtype=float)

y = substituicao_progressiva(L_ex, b_ex)
print("L_ex =")
print(L_ex)
print("b_ex =", b_ex)
print("y =", y)
print("verificação L_ex @ y =", L_ex @ y)


titulo("3.4) substituicao_regressiva(U, y)")

# Exemplo triangular superior:
U_ex = np.array([
    [2, -1, 1],
    [0, 3, -2],
    [0, 0, 4]
], dtype=float)
y_ex = np.array([2, 5, 8], dtype=float)

x = substituicao_regressiva(U_ex, y_ex)
print("U_ex =")
print(U_ex)
print("y_ex =", y_ex)
print("x =", x)
print("verificação U_ex @ x =", U_ex @ x)


titulo("3.5) resolver_lu(A, b)")

x = resolver_lu(A, b)
print("x por LU =", x)
print("A @ x =", A @ x)


titulo("3.6) thomas_algorithm(a, b, c, d)")

# Sistema tridiagonal:
# 2x1 - x2      = 1
# -x1 + 2x2 -x3 = 0
#      -x2 +2x3 = 1
#
# Solução: [1, 1, 1]
a_tri = [0, -1, -1]   # subdiagonal, a[0]=0
b_tri = [2, 2, 2]     # diagonal principal
c_tri = [-1, -1, 0]   # superdiagonal, c[-1]=0
d_tri = [1, 0, 1]     # lado direito

x = thomas_algorithm(a_tri, b_tri, c_tri, d_tri)
print("x =", x)


titulo("3.7) cholesky_solve(A, b)")

# Cholesky exige matriz simétrica definida positiva.
A_spd = np.array([
    [4, 1],
    [1, 3]
], dtype=float)
b_spd = np.array([1, 2], dtype=float)

x = cholesky_solve(A_spd, b_spd)
print("x =", x)
print("A_spd @ x =", A_spd @ x)


titulo("3.8) gauss_seidel(A, b, x0, lamb)")

# Matriz diagonalmente dominante, boa para Gauss-Seidel.
A_gs = np.array([
    [10, 2, 1],
    [1, 5, 1],
    [2, 3, 10]
], dtype=float)
b_gs = np.array([7, -8, 6], dtype=float)

x = gauss_seidel(A_gs, b_gs, x0=[0, 0, 0], tol=1e-6, lamb=1.0, verbose=False)
print("x sem relaxação =", x)
print("A_gs @ x =", A_gs @ x)

x_relax = gauss_seidel(A_gs, b_gs, x0=[0, 0, 0], tol=1e-6, lamb=1.2, verbose=False)
print("x com relaxação lambda=1.2 =", x_relax)


# ============================================================
# 4) OTIMIZAÇÃO 1D E 2D
# ============================================================

titulo("4.1) razao_aurea(f, xl, xu, modo='max')")

# Maximizar f(x) = -(x - 2)^2 + 4
# Máximo em x = 2, f = 4
f1d = lambda x: -(x - 2)**2 + 4
xopt, fopt = razao_aurea(f1d, xl=0, xu=5, tol=1e-6, modo="max", verbose=False)
print("x ótimo =", xopt)
print("f ótimo =", fopt)


titulo("4.2) interpolacao_quadratica_1d(f, x0, x1, x2, modo='max')")

# Mesma função. Usamos três pontos que cercam o máximo.
f1d = lambda x: -(x - 2)**2 + 4
xopt, fopt = interpolacao_quadratica_1d(f1d, x0=0, x1=1, x2=4, tol=1e-8, modo="max", verbose=False)
print("x ótimo por interpolação quadrática =", xopt)
print("f ótimo =", fopt)


titulo("4.3) aclive_maximo(f, grad, z0)")

# Exemplo do PPC4:
# Maximizar f(x,y) = 2xy + 2x - x^2 - 2y^2
# Ótimo analítico: (2, 1), f = 2

def f2d(z):
    x, y = z
    return 2*x*y + 2*x - x**2 - 2*y**2

def grad2d(z):
    x, y = z
    return np.array([
        2*y + 2 - 2*x,
        2*x - 4*y
    ], dtype=float)

# Neste kit, o aclive_maximo usa uma interpolação quadrática simples.
# Para este exemplo quadrático, limitamos a 5 iterações para evitar degeneração numérica
# quando os três pontos da interpolação ficam praticamente coincidentes perto do ótimo.
sol, log = aclive_maximo(f2d, grad2d, z0=[-2, 3], tol=1e-8, max_iter=5, verbose=False)
print("solução aproximada por aclive máximo =", sol)
print("f(sol) =", f2d(sol))
print("iterações registradas =", len(log))


titulo("4.4) fletcher_reeves_max(f, grad, z0)")

sol, log = fletcher_reeves_max(f2d, grad2d, z0=[-2, 3], tol=1e-8, verbose=False)
print("solução por Fletcher-Reeves =", sol)
print("f(sol) =", f2d(sol))
print("iterações registradas =", len(log))


# ============================================================
# 5) EDO/PVI: EULER E RK4
# ============================================================

titulo("5.1) euler(f, t0, y0, h, n)")

# PVI: dy/dt = y - t^2 + 1, y(0)=0.5
# Exemplo clássico de EDO de 1ª ordem.
def edo(t, y):
    return y - t**2 + 1

ts, ys = euler(edo, t0=0, y0=[0.5], h=0.2, n=5)
print("t =", ts)
print("y por Euler =", ys[:, 0])


titulo("5.2) rk4(f, t0, y0, h, n)")

ts, ys = rk4(edo, t0=0, y0=[0.5], h=0.2, n=5)
print("t =", ts)
print("y por RK4 =", ys[:, 0])


titulo("5.3) rk4 para sistema de EDOs")

# Sistema massa-mola simples:
# x' = v
# v' = -x
# Estado y = [x, v]
# Condição inicial: x(0)=1, v(0)=0

def oscilador(t, y):
    x, v = y
    return np.array([v, -x], dtype=float)

ts, ys = rk4(oscilador, t0=0, y0=[1, 0], h=0.1, n=10)
print("último t =", ts[-1])
print("último estado [x, v] =", ys[-1])


print("\nTodos os exemplos rodaram.")
