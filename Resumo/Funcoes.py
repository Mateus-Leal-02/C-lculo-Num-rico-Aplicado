"""
Cálculo Numérico Aplicado
Autor: Mateus Leal Silva - versão de revisão rápida

Ideia: deixar funções prontas no GitHub e, na prova, COPIAR/ADAPTAR.
Dependência: numpy.
Como rodar:
    python3 Resumo.py
Como usar interativamente:
    python3 -i Resumo.py

"""

import numpy as np

# ============================================================
# 0) FUNÇÕES ÚTEIS
# ============================================================

def norma(v):
    """Norma euclidiana."""
    v = np.array(v, dtype=float)
    return float(np.sqrt(np.dot(v, v)))


def erro_relativo_atual(novo, antigo):
    """Erro relativo percentual aproximado. Evita divisão por zero."""
    novo = np.array(novo, dtype=float)
    antigo = np.array(antigo, dtype=float)
    den = np.maximum(np.abs(novo), 1e-15)
    return np.max(np.abs((novo - antigo) / den)) * 100.0


# ============================================================
# 1) RAÍZES DE FUNÇÕES: bisseção, falsa posição, Newton, secante
# ============================================================

def bissecao(f, a, b, tol=1e-8, max_iter=100, verbose=True):
    """
    Resolve f(x)=0 em [a,b]. Exige f(a)*f(b)<0.
    Vantagem: robusto. Desvantagem: mais lento.
    """
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError("Bisseção precisa de intervalo com troca de sinal: f(a)*f(b)<0")

    xr_old = None
    for k in range(1, max_iter + 1):
        xr = (a + b) / 2.0
        fr = f(xr)
        ea = None if xr_old is None else abs((xr - xr_old) / max(abs(xr), 1e-15)) * 100.0

        if verbose:
            print(f"{k:3d}  a={a:.8g}  b={b:.8g}  xr={xr:.12g}  f(xr)={fr:.3e}  ea={ea}")

        if abs(fr) < tol or (ea is not None and ea < tol):
            return xr

        if fa * fr < 0:
            b, fb = xr, fr
        else:
            a, fa = xr, fr
        xr_old = xr
    return xr


def falsa_posicao(f, a, b, tol=1e-8, max_iter=100, verbose=True):
    """Método da falsa posição/regula falsi. Exige f(a)*f(b)<0."""
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError("Falsa posição precisa de intervalo com troca de sinal")

    xr_old = None
    for k in range(1, max_iter + 1):
        xr = b - fb * (a - b) / (fa - fb)
        fr = f(xr)
        ea = None if xr_old is None else abs((xr - xr_old) / max(abs(xr), 1e-15)) * 100.0

        if verbose:
            print(f"{k:3d}  xr={xr:.12g}  f(xr)={fr:.3e}  ea={ea}")

        if abs(fr) < tol or (ea is not None and ea < tol):
            return xr

        if fa * fr < 0:
            b, fb = xr, fr
        else:
            a, fa = xr, fr
        xr_old = xr
    return xr


def newton_raphson(f, df, x0, tol=1e-8, max_iter=50, verbose=True):
    """Newton-Raphson: x_{k+1}=x_k-f(x_k)/f'(x_k). Rápido, mas depende de x0."""
    x = float(x0)
    for k in range(1, max_iter + 1):
        dfx = df(x)
        if abs(dfx) < 1e-15:
            raise ZeroDivisionError("Derivada quase zero no Newton-Raphson")
        x_new = x - f(x) / dfx
        ea = abs((x_new - x) / max(abs(x_new), 1e-15)) * 100.0

        if verbose:
            print(f"{k:3d}  x={x_new:.12g}  f(x)={f(x_new):.3e}  ea={ea:.3e}%")

        if abs(f(x_new)) < tol or ea < tol:
            return x_new
        x = x_new
    return x


def secante(f, x0, x1, tol=1e-8, max_iter=50, verbose=True):
    """Método da secante: não precisa da derivada, mas precisa de dois chutes."""
    x_ant, x = float(x0), float(x1)
    for k in range(1, max_iter + 1):
        f_ant, fx = f(x_ant), f(x)
        if abs(f_ant - fx) < 1e-15:
            raise ZeroDivisionError("Denominador quase zero na secante")
        x_new = x - fx * (x_ant - x) / (f_ant - fx)
        ea = abs((x_new - x) / max(abs(x_new), 1e-15)) * 100.0

        if verbose:
            print(f"{k:3d}  x={x_new:.12g}  f(x)={f(x_new):.3e}  ea={ea:.3e}%")

        if abs(f(x_new)) < tol or ea < tol:
            return x_new
        x_ant, x = x, x_new
    return x


def newton_sistema(F, J, x0, tol=1e-8, max_iter=30, verbose=True):
    """
    Newton para sistemas não-lineares:
        J(x_k) * delta = -F(x_k)
        x_{k+1} = x_k + delta
    """
    x = np.array(x0, dtype=float)
    for k in range(1, max_iter + 1):
        delta = np.linalg.solve(J(x), -F(x))
        x_new = x + delta
        erro = norma(delta)

        if verbose:
            print(f"{k:3d}  x={x_new}  ||delta||={erro:.3e}  ||F||={norma(F(x_new)):.3e}")

        if erro < tol or norma(F(x_new)) < tol:
            return x_new
        x = x_new
    return x


# ============================================================
# 2) RAÍZES DE POLINÔMIOS: Bairstow
# ============================================================

def bairstow(coef, r=0.0, s=0.0, tol=1e-8, max_iter=100, verbose=True):
    """
    Método de Bairstow para raízes de polinômios reais.
    coef deve estar em ordem decrescente: [a_n, a_{n-1}, ..., a_0].
    Fatora termos quadráticos x^2 - r*x - s.
    Retorna uma lista de raízes, possivelmente complexas.
    """
    a = np.array(coef, dtype=float)
    roots = []

    while len(a) > 3:
        n = len(a) - 1
        rr, ss = float(r), float(s)
        for it in range(1, max_iter + 1):
            b = np.zeros(n + 1)
            c = np.zeros(n + 1)

            b[0] = a[0]
            b[1] = a[1] + rr * b[0]
            for i in range(2, n + 1):
                b[i] = a[i] + rr * b[i - 1] + ss * b[i - 2]

            c[0] = b[0]
            c[1] = b[1] + rr * c[0]
            for i in range(2, n + 1):
                c[i] = b[i] + rr * c[i - 1] + ss * c[i - 2]

            # Sistema para correções dr e ds
            # [c[n-1] c[n-2]; c[n] c[n-1]] [dr; ds] = [-b[n]; -b[n-1]]
            M = np.array([[c[n - 1], c[n - 2]], [c[n], c[n - 1]]], dtype=float)
            rhs = np.array([-b[n], -b[n - 1]], dtype=float)
            try:
                dr, ds = np.linalg.solve(M, rhs)
            except np.linalg.LinAlgError:
                rr += 0.1
                ss += 0.1
                continue

            rr += dr
            ss += ds
            erro = max(abs(dr / max(abs(rr), 1e-15)), abs(ds / max(abs(ss), 1e-15))) * 100.0

            if verbose:
                print(f"Bairstow it={it:3d}  r={rr:.8g}  s={ss:.8g}  erro={erro:.3e}%")

            if erro < tol:
                break

        # Raízes de x^2 - r*x - s = 0
        disc = rr**2 + 4.0 * ss
        roots.extend(np.roots([1.0, -rr, -ss]))

        # Deflação: b[0:n-1] são coeficientes do quociente
        a = b[:n - 1]

    if len(a) == 3:
        roots.extend(np.roots(a))
    elif len(a) == 2:
        roots.append(-a[1] / a[0])

    return np.array(roots)


# ============================================================
# 3) SISTEMAS LINEARES
# ============================================================

def eliminacao_gauss_pivot(A, b):
    """Eliminação de Gauss com pivotamento parcial."""
    A = np.array(A, dtype=float).copy()
    b = np.array(b, dtype=float).copy()
    n = len(b)

    # Eliminação progressiva
    for k in range(n - 1):
        pivo = k + np.argmax(np.abs(A[k:, k]))
        if abs(A[pivo, k]) < 1e-15:
            raise ZeroDivisionError("Matriz singular ou pivô nulo")
        if pivo != k:
            A[[k, pivo]] = A[[pivo, k]]
            b[[k, pivo]] = b[[pivo, k]]
        for i in range(k + 1, n):
            fator = A[i, k] / A[k, k]
            A[i, k:] -= fator * A[k, k:]
            b[i] -= fator * b[k]

    # Substituição regressiva
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        soma = np.dot(A[i, i + 1:], x[i + 1:])
        x[i] = (b[i] - soma) / A[i, i]
    return x


def decomposicao_lu(A):
    """LU sem pivotamento: A=L@U. Use quando não houver pivô zero."""
    A = np.array(A, dtype=float)
    n = A.shape[0]
    L = np.eye(n)
    U = np.zeros_like(A)

    for k in range(n):
        for j in range(k, n):
            U[k, j] = A[k, j] - np.dot(L[k, :k], U[:k, j])
        if abs(U[k, k]) < 1e-15:
            raise ZeroDivisionError("Pivô zero na LU sem pivotamento")
        for i in range(k + 1, n):
            L[i, k] = (A[i, k] - np.dot(L[i, :k], U[:k, k])) / U[k, k]
    return L, U


def substituicao_progressiva(L, b):
    L = np.array(L, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    y = np.zeros(n)
    for i in range(n):
        y[i] = (b[i] - np.dot(L[i, :i], y[:i])) / L[i, i]
    return y


def substituicao_regressiva(U, y):
    U = np.array(U, dtype=float)
    y = np.array(y, dtype=float)
    n = len(y)
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i + 1:], x[i + 1:])) / U[i, i]
    return x


def resolver_lu(A, b):
    L, U = decomposicao_lu(A)
    y = substituicao_progressiva(L, b)
    return substituicao_regressiva(U, y)


def thomas_algorithm(a, b, c, d):
    """
    TDMA/Thomas para sistema tridiagonal:
        a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
    Convenção: a[0]=0 e c[-1]=0.
    """
    a = np.array(a, dtype=float).copy()
    b = np.array(b, dtype=float).copy()
    c = np.array(c, dtype=float).copy()
    d = np.array(d, dtype=float).copy()
    n = len(d)

    c_star = np.zeros(n)
    d_star = np.zeros(n)
    x = np.zeros(n)

    c_star[0] = c[0] / b[0]
    d_star[0] = d[0] / b[0]

    for i in range(1, n):
        m = b[i] - a[i] * c_star[i - 1]
        if abs(m) < 1e-15:
            raise ZeroDivisionError("Pivô zero no algoritmo de Thomas")
        if i < n - 1:
            c_star[i] = c[i] / m
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m

    x[-1] = d_star[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d_star[i] - c_star[i] * x[i + 1]
    return x


def cholesky_solve(A, b):
    """Resolve A x=b usando Cholesky. A deve ser simétrica definida positiva."""
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    L = np.linalg.cholesky(A)
    y = substituicao_progressiva(L, b)
    x = substituicao_regressiva(L.T, y)
    return x


def gauss_seidel(A, b, x0=None, tol=1e-8, max_iter=500, lamb=1.0, verbose=True):
    """
    Gauss-Seidel com relaxação.
    lamb=1 sem relaxação; lamb>1 sobre-relaxação; lamb<1 sub-relaxação.
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)

    for k in range(1, max_iter + 1):
        x_old = x.copy()
        for i in range(n):
            soma1 = np.dot(A[i, :i], x[:i])          # valores novos
            soma2 = np.dot(A[i, i + 1:], x_old[i + 1:])  # valores antigos
            x_gs = (b[i] - soma1 - soma2) / A[i, i]
            x[i] = lamb * x_gs + (1.0 - lamb) * x_old[i]

        erro = erro_relativo_atual(x, x_old)
        if verbose:
            print(f"{k:3d}  x={x}  ea={erro:.3e}%")
        if erro < tol:
            return x
    return x


# ============================================================
# 4) OTIMIZAÇÃO 1D E 2D
# ============================================================

def razao_aurea(f, xl, xu, tol=1e-8, max_iter=100, modo="min", verbose=True):
    """
    Busca da razão áurea para mínimo ou máximo em [xl,xu].
    modo='min' ou modo='max'.
    """
    R = (np.sqrt(5.0) - 1.0) / 2.0
    d = R * (xu - xl)
    x1 = xl + d
    x2 = xu - d
    f1, f2 = f(x1), f(x2)

    for k in range(1, max_iter + 1):
        if (modo == "min" and f1 < f2) or (modo == "max" and f1 > f2):
            xl = x2
            x2 = x1
            f2 = f1
            d = R * (xu - xl)
            x1 = xl + d
            f1 = f(x1)
        else:
            xu = x1
            x1 = x2
            f1 = f2
            d = R * (xu - xl)
            x2 = xu - d
            f2 = f(x2)

        xopt = x1 if ((modo == "min" and f1 < f2) or (modo == "max" and f1 > f2)) else x2
        erro = (1.0 - R) * abs((xu - xl) / max(abs(xopt), 1e-15)) * 100.0
        if verbose:
            print(f"{k:3d}  xl={xl:.8g}  xu={xu:.8g}  xopt={xopt:.12g}  f={f(xopt):.8g}  ea={erro:.3e}%")
        if erro < tol:
            return xopt, f(xopt)
    return xopt, f(xopt)


def interpolacao_quadratica_1d(f, x0, x1, x2, tol=1e-8, max_iter=50, modo="max", verbose=True):
    """
    Interpolação quadrática com três pontos. Retorna vértice da parábola.
    Útil para linha de busca em otimização.
    """
    for k in range(1, max_iter + 1):
        f0, f1, f2 = f(x0), f(x1), f(x2)
        num = f0 * (x1**2 - x2**2) + f1 * (x2**2 - x0**2) + f2 * (x0**2 - x1**2)
        den = 2*f0*(x1 - x2) + 2*f1*(x2 - x0) + 2*f2*(x0 - x1)
        if abs(den) < 1e-15:
            raise ZeroDivisionError("Denominador nulo na interpolação quadrática")
        x3 = num / den

        if verbose:
            print(f"{k:3d}  x3={x3:.12g}  f(x3)={f(x3):.8g}")

        if abs(x3 - x1) < tol:
            return x3, f(x3)

        # Mantém três pontos em torno do melhor valor conhecido
        pts = sorted([(x0, f0), (x1, f1), (x2, f2), (x3, f(x3))], key=lambda p: p[0])
        melhor_idx = np.argmax([p[1] for p in pts]) if modo == "max" else np.argmin([p[1] for p in pts])
        i0 = max(0, melhor_idx - 1)
        i2 = min(len(pts) - 1, melhor_idx + 1)
        if i2 - i0 < 2:
            if i0 == 0:
                i2 = 2
            else:
                i0 = len(pts) - 3
        x0, x1, x2 = pts[i0][0], pts[melhor_idx][0], pts[i2][0]
    return x3, f(x3)


def aclive_maximo(f, grad, z0, tol=1e-8, max_iter=100, verbose=True):
    """Steepest ascent / aclive máximo para maximização sem restrições."""
    z = np.array(z0, dtype=float)
    log = []
    for k in range(max_iter):
        g = np.array(grad(z), dtype=float)
        erro = norma(g)
        if erro < tol:
            break
        p = g

        def phi(h):
            return f(z + h * p)

        # três pontos para linha de busca; ajuste se o passo explodir
        h, _ = interpolacao_quadratica_1d(phi, 0.0, 1.0, 2.0, tol=1e-12, max_iter=3, modo="max", verbose=False)
        log.append([k, erro, h, *z, *g])
        if verbose:
            print(f"{k:3d}  erro={erro:.3e}  h={h:.8g}  z={z}  grad={g}")
        z = z + h * p
    return z, np.array(log, dtype=float)


def fletcher_reeves_max(f, grad, z0, tol=1e-8, max_iter=100, verbose=True):
    """Gradientes conjugados de Fletcher-Reeves para maximização."""
    z = np.array(z0, dtype=float)
    g = np.array(grad(z), dtype=float)
    p = g.copy()
    log = []

    for k in range(max_iter):
        erro = norma(g)
        if erro < tol:
            break

        def phi(h):
            return f(z + h * p)

        h, _ = interpolacao_quadratica_1d(phi, 0.0, 1.0, 2.0, tol=1e-12, max_iter=3, modo="max", verbose=False)
        log.append([k, erro, h, *z, *g])
        if verbose:
            print(f"{k:3d}  erro={erro:.3e}  h={h:.8g}  z={z}  grad={g}")

        z_new = z + h * p
        g_new = np.array(grad(z_new), dtype=float)
        beta = np.dot(g_new, g_new) / max(np.dot(g, g), 1e-15)
        p = g_new + beta * p
        z, g = z_new, g_new
    return z, np.array(log, dtype=float)


# ============================================================
# 5) EDO/PVI: Euler e Runge-Kutta de 4ª ordem
# ============================================================

def euler(f, t0, y0, h, n):
    """Euler explícito: y_{i+1}=y_i+h*f(t_i,y_i)."""
    t = float(t0)
    y = np.array(y0, dtype=float)
    ts = [t]
    ys = [y.copy()]
    for _ in range(n):
        y = y + h * np.array(f(t, y), dtype=float)
        t = t + h
        ts.append(t)
        ys.append(y.copy())
    return np.array(ts), np.array(ys)


def rk4(f, t0, y0, h, n):
    """Runge-Kutta clássico de 4ª ordem para escalar ou sistema."""
    t = float(t0)
    y = np.array(y0, dtype=float)
    ts = [t]
    ys = [y.copy()]
    for _ in range(n):
        k1 = np.array(f(t, y), dtype=float)
        k2 = np.array(f(t + h/2.0, y + h*k1/2.0), dtype=float)
        k3 = np.array(f(t + h/2.0, y + h*k2/2.0), dtype=float)
        k4 = np.array(f(t + h, y + h*k3), dtype=float)
        y = y + (h/6.0) * (k1 + 2*k2 + 2*k3 + k4)
        t = t + h
        ts.append(t)
        ys.append(y.copy())
    return np.array(ts), np.array(ys)


# ============================================================
# 6) EXEMPLOS RÁPIDOS PARA TESTAR SE O ARQUIVO ESTÁ FUNCIONANDO
# ============================================================

if __name__ == "__main__":
    print("\n=== Teste 1: raiz por bisseção de x^3 - x - 2 ===")
    f_raiz = lambda x: x**3 - x - 2
    raiz = bissecao(f_raiz, 1, 2, tol=1e-6, verbose=False)
    print("raiz ≈", raiz, "f(raiz)=", f_raiz(raiz))

    print("\n=== Teste 2: sistema tridiagonal por Thomas ===")
    a = [0, -1, -1]
    b = [2, 2, 2]
    c = [-1, -1, 0]
    d = [1, 0, 1]
    print("x =", thomas_algorithm(a, b, c, d))

    print("\n=== Teste 3: otimização PPC4 ===")
    def f2(z):
        x, y = z
        return 2*x*y + 2*x - x**2 - 2*y**2

    def grad2(z):
        x, y = z
        return np.array([2*y + 2 - 2*x, 2*x - 4*y], dtype=float)

    sol, log = fletcher_reeves_max(f2, grad2, [-2, 3], verbose=False)
    print("solução ≈", sol, "f=", f2(sol))

    print("\nArquivo OK. Na prova, copie a função do método pedido e adapte f, grad, A, b, etc.")
