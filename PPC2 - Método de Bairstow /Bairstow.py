"""
=============================================================================
  Determinação de Raízes de Polinômios — Método de Bairstow
=============================================================================

Descrição:
    Este programa implementa o método de Bairstow para determinação de todas
    as raízes (reais e complexas conjugadas) de um polinômio de grau arbitrário.
    O método divide iterativamente o polinômio por fatores quadráticos reais
    x² - rx - s, ajustando (r, s) pelo método de Newton-Raphson até que o
    resto da divisão sintética seja nulo.

    A implementação é aplicada a dois problemas:
      (a) Validação com polinômio de grau 7 de raízes conhecidas.
      (b) Polinômio característico do sistema 2-GDL da APC2:
              P(λ) = 2λ⁴ + 5λ³ + 12λ² + 8λ + 8

    Adicionalmente, gera-se o fractal de Bairstow: um mapa de convergência
    no plano (r₀, s₀) que ilustra a sensibilidade do método ao chute inicial.

Referências:
    - Chapra & Canale, Métodos Numéricos para Engenharia, 5ª ed., McGraw-Hill (2008)
    - Bairstow, L., Applied Aerodynamics, Longmans, 1920 (Apêndice)

Autor: Mateus Leal Silva | Matrícula: 221028134
Disciplina: Cálculo Numérico Aplicado — Prof. Dr. Rafael Gabler Gontijo
Data: abril de 2026
=============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import warnings

# ---------------------------------------------------------------------------
# 1. DIVISÃO SINTÉTICA — COEFICIENTES b
# ---------------------------------------------------------------------------

def _calc_b(coeffs, r, s):
    """
    Divide o polinômio f(x) pelo fator quadrático (x² - rx - s) e retorna
    o vetor de coeficientes b via divisão sintética.

    Convenção de armazenamento (potências decrescentes):
        coeffs = [a_n, a_{n-1}, ..., a_1, a_0]  →  f(x) = a_n x^n + ...

    Relações de recorrência (equação 13 do enunciado):
        b[0] = a_n
        b[1] = a_{n-1} + r · b[0]
        b[j] = a_{n-j} + r · b[j-1] + s · b[j-2],  j = 2, ..., n

    O resto da divisão é: b[n-1]·(x − r) + b[n],
    ou seja, b_1 = b[n-1]  e  b_0 = b[n]  na notação do enunciado.

    Parâmetros
    ----------
    coeffs : array-like — coeficientes do polinômio (potências decrescentes)
    r, s   : float      — parâmetros do divisor quadrático x² - rx - s

    Retorna
    -------
    b : np.ndarray — vetor de coeficientes b (tamanho n+1)
    """
    n = len(coeffs) - 1
    b = np.zeros(n + 1)
    b[0] = coeffs[0]
    if n >= 1:
        b[1] = coeffs[1] + r * b[0]
    for j in range(2, n + 1):
        b[j] = coeffs[j] + r * b[j - 1] + s * b[j - 2]
    return b


# ---------------------------------------------------------------------------
# 2. DIVISÃO SINTÉTICA — COEFICIENTES c (derivadas parciais de b)
# ---------------------------------------------------------------------------

def _calc_c(b, r, s):
    """
    Aplica uma segunda divisão sintética sobre os coeficientes b para obter
    os coeficientes c, que representam as derivadas parciais de b₀ e b₁
    em relação a r e s (equação 16 do enunciado):

        c[0] = b[0]
        c[1] = b[1] + r · c[0]
        c[j] = b[j] + r · c[j-1] + s · c[j-2],  j = 2, ..., n-1

    Mapeamento para as derivadas (equação 17 do enunciado):
        ∂b₀/∂r = c[n-1]  (≡ c₁ na notação do enunciado)
        ∂b₀/∂s = ∂b₁/∂r = c[n-2]  (≡ c₂)
        ∂b₁/∂s = c[n-3]  (≡ c₃)

    Parâmetros
    ----------
    b    : np.ndarray — coeficientes b (saída de _calc_b)
    r, s : float      — parâmetros do divisor quadrático

    Retorna
    -------
    c : np.ndarray — vetor de coeficientes c (tamanho n, índices 0 a n-1)
    """
    n = len(b) - 1
    c = np.zeros(n)
    c[0] = b[0]
    if n >= 2:
        c[1] = b[1] + r * c[0]
    for j in range(2, n):
        c[j] = b[j] + r * c[j - 1] + s * c[j - 2]
    return c


# ---------------------------------------------------------------------------
# 3. UMA ITERAÇÃO DO MÉTODO DE BAIRSTOW (passo de Newton-Raphson)
# ---------------------------------------------------------------------------

def _passo_bairstow(coeffs, r, s):
    """
    Executa um passo do método de Newton-Raphson sobre os parâmetros (r, s)
    para reduzir o resto da divisão sintética (b₀ e b₁) a zero.

    Sistema 2×2 resolvido pela regra de Cramer (equação 15 do enunciado):

        ⎡ c₂  c₃ ⎤ ⎡ Δr ⎤   ⎡ -b₁ ⎤
        ⎣ c₁  c₂ ⎦ ⎣ Δs ⎦ = ⎣ -b₀ ⎦

    Solução:
        det = c₂² - c₁·c₃
        Δr  = (b₀·c₃ - b₁·c₂) / det
        Δs  = (b₁·c₁ - b₀·c₂) / det

    Parâmetros
    ----------
    coeffs : np.ndarray — coeficientes do polinômio (potências decrescentes)
    r, s   : float      — estimativa atual dos parâmetros

    Retorna
    -------
    (delta_r, delta_s) : floats  — passos de atualização, ou (None, None) se
                                   o sistema linear for degenerado
    b0, b1             : floats  — restos da divisão sintética atual
    """
    n = len(coeffs) - 1
    b = _calc_b(coeffs, r, s)
    c = _calc_c(b, r, s)

    b0 = b[n]       # termo independente do resto (≡ b₀ no enunciado)
    b1 = b[n - 1]   # coeficiente linear do resto  (≡ b₁ no enunciado)

    # c₁, c₂, c₃ na notação do enunciado (mapeados de índices Python)
    c1 = c[n - 1]   # ∂b₀/∂r
    c2 = c[n - 2]   # ∂b₀/∂s = ∂b₁/∂r
    c3 = c[n - 3]   # ∂b₁/∂s

    det = c2 ** 2 - c1 * c3

    if abs(det) < 1e-30:
        return None, None, b1, b0

    delta_r = (b0 * c3 - b1 * c2) / det
    delta_s = (b1 * c1 - b0 * c2) / det

    return delta_r, delta_s, b1, b0


# ---------------------------------------------------------------------------
# 4. MÉTODO DE BAIRSTOW — IMPLEMENTAÇÃO COMPLETA
# ---------------------------------------------------------------------------

def bairstow(coeffs, r0=0.5, s0=0.5, tol=1e-10, max_iter=1000):
    """
    Método de Bairstow para determinação de todas as raízes de um polinômio
    de grau arbitrário, incluindo raízes complexas conjugadas.

    Algoritmo (resumo do fluxograma do enunciado):
      1. Supor valores iniciais (r₀, s₀).
      2. Calcular coeficientes b e c por divisão sintética.
      3. Resolver sistema 2×2 para (Δr, Δs).
      4. Atualizar r ← r + Δr, s ← s + Δs.
      5. Repetir até |Δr/r| e |Δs/s| < tol.
      6. Extrair par de raízes via fórmula quadrática em x² - rx - s.
      7. Deflacionar o polinômio e repetir para o quociente.

    Parâmetros
    ----------
    coeffs   : array-like — coeficientes [a_n, ..., a_0] (potências decrescentes)
    r0, s0   : float      — chutes iniciais para os parâmetros do divisor
    tol      : float      — tolerância relativa de convergência
    max_iter : int        — máximo de iterações por deflação

    Retorna
    -------
    roots  : np.ndarray (complex) — todas as raízes do polinômio
    iters  : list of int          — nº de iterações por passo de deflação
    """
    coeffs = np.array(coeffs, dtype=float)
    coeffs = coeffs / coeffs[0]           # normalizar pelo coeficiente líder

    roots = []
    iter_counts = []
    poly = coeffs.copy()

    while len(poly) > 3:
        n_local = len(poly) - 1
        r, s = r0, s0
        converged = False

        for it in range(max_iter):
            dr, ds, b1, b0 = _passo_bairstow(poly, r, s)

            if dr is None:
                # Sistema degenerado: pequena perturbação e continua
                r += 1e-6
                s += 1e-6
                continue

            r += dr
            s += ds

            # Critério de parada: convergência relativa
            crit_r = abs(dr / r) if abs(r) > 1e-14 else abs(dr)
            crit_s = abs(ds / s) if abs(s) > 1e-14 else abs(ds)

            if crit_r < tol and crit_s < tol:
                converged = True
                break

        if not converged:
            warnings.warn(
                f"Bairstow não convergiu para o fator (r={r:.4f}, s={s:.4f}) "
                f"após {max_iter} iterações. Tente outros valores de r0, s0."
            )

        iter_counts.append(it + 1)

        # --- Extrair par de raízes via fórmula quadrática: x² - rx - s = 0 ---
        discriminante = r ** 2 + 4.0 * s
        if discriminante >= 0:
            raiz1 = (r + np.sqrt(discriminante)) / 2.0
            raiz2 = (r - np.sqrt(discriminante)) / 2.0
        else:
            parte_real = r / 2.0
            parte_imag = np.sqrt(-discriminante) / 2.0
            raiz1 = complex(parte_real,  parte_imag)
            raiz2 = complex(parte_real, -parte_imag)

        roots.extend([raiz1, raiz2])

        # --- Deflação: substituir poly pelo quociente de grau n_local - 2 ---
        b = _calc_b(poly, r, s)
        poly = b[:n_local - 1]            # coeficientes b[0] ... b[n_local-2]

    # --- Tratar polinômio residual ---
    if len(poly) == 3:
        # Grau 2: fórmula de Bhaskara
        a2, a1, a0 = poly
        disc = a1 ** 2 - 4.0 * a2 * a0
        if disc >= 0:
            roots.append((-a1 + np.sqrt(disc)) / (2.0 * a2))
            roots.append((-a1 - np.sqrt(disc)) / (2.0 * a2))
        else:
            roots.append((-a1 + 1j * np.sqrt(-disc)) / (2.0 * a2))
            roots.append((-a1 - 1j * np.sqrt(-disc)) / (2.0 * a2))
        iter_counts.append(0)

    elif len(poly) == 2:
        # Grau 1: raiz direta  xr = -a0 / a1
        roots.append(-poly[1] / poly[0])
        iter_counts.append(0)

    return np.array(roots, dtype=complex), iter_counts


# ---------------------------------------------------------------------------
# 5. ANÁLISE 1 — Validação com polinômio de grau 7 (raízes conhecidas)
# ---------------------------------------------------------------------------

def analise_validacao(salvar=True):
    """
    Constrói um polinômio mônico de grau 7 a partir de raízes pré-definidas
    (3 reais + 2 pares complexos conjugados), aplica o método de Bairstow e
    compara as raízes recuperadas com as originais.

    Raízes usadas:  1, -2, 3,  (1±2i),  (-1±i)
    """
    raizes_ref = np.array([
        1.0, -2.0, 3.0,
        complex( 1.0,  2.0), complex( 1.0, -2.0),
        complex(-1.0,  1.0), complex(-1.0, -1.0)
    ])

    # np.poly constrói o polinômio mônico a partir das raízes
    coeffs = np.poly(raizes_ref)  # coeficientes em potências decrescentes

    raizes_calc, iters = bairstow(coeffs, r0=0.5, s0=-0.5, tol=1e-10)

    # Ordenar ambos os conjuntos para comparação
    chave = lambda z: (round(z.real, 6), round(z.imag, 6))
    ref_ord  = sorted(raizes_ref,  key=chave)
    calc_ord = sorted(raizes_calc, key=chave)

    print("\n" + "=" * 70)
    print("  ANÁLISE 1 — Validação com polinômio de grau 7")
    print("=" * 70)
    print(f"  Coeficientes: {np.round(coeffs, 4)}")
    print(f"  Iterações por deflação: {iters}")
    print()
    print(f"{'Raiz conhecida':>30}  {'Raiz calculada':>30}  {'Erro abs':>12}")
    print("-" * 76)
    for zk, zc in zip(ref_ord, calc_ord):
        print(f"{str(np.round(zk, 6)):>30}  {str(np.round(zc, 6)):>30}  {abs(zk-zc):12.3e}")

    # --- Gráfico no plano complexo ---
    fig, ax = plt.subplots(figsize=(7, 6))

    rk = np.array(raizes_ref)
    rc = np.array(raizes_calc)

    ax.scatter(rk.real, rk.imag, s=120, marker="o", color="#d7191c",
               zorder=5, label="Raízes conhecidas", edgecolors="k", linewidths=0.8)
    ax.scatter(rc.real, rc.imag, s=70, marker="x", color="#2c7bb6",
               zorder=6, label="Raízes Bairstow", linewidths=2.5)

    ax.axhline(0, color="gray", lw=0.8, ls="--", alpha=0.6)
    ax.axvline(0, color="gray", lw=0.8, ls="--", alpha=0.6)
    ax.set_xlabel(r"Re($x$)", fontsize=13)
    ax.set_ylabel(r"Im($x$)", fontsize=13)
    ax.set_title("Validação do método de Bairstow — polinômio grau 7", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if salvar:
        os.makedirs("resultados", exist_ok=True)
        fig.savefig("resultados/analise_1_validacao.png", dpi=150)
        print("\n[✓] resultados/analise_1_validacao.png salvo.")
    plt.show()


# ---------------------------------------------------------------------------
# 6. ANÁLISE 2 — Polinômio característico do sistema dinâmico (APC2)
# ---------------------------------------------------------------------------

def analise_apc2(salvar=True):
    """
    Aplica o método de Bairstow ao polinômio característico do sistema
    massa-mola-amortecedor com N=2 GDL obtido na APC2:

        P(λ) = 2λ⁴ + 5λ³ + 12λ² + 8λ + 8

    Parâmetros do sistema (APC2):
        m₁ = 2 kg,  m₂ = 1 kg
        k₁ = 4 N/m, k₂ = 2 N/m
        c₁ = 2 N·s/m, c₂ = 1 N·s/m

    Interpreta os autovalores λ = σ ± iωd em termos de:
      - σ < 0 → modo amortecido (estável)
      - ωd    → frequência de oscilação amortecida [rad/s adimensional]
      - ζ     → coeficiente de amortecimento modal
    """
    # Coeficientes em potências decrescentes: [a4, a3, a2, a1, a0]
    coeffs_apc2 = [2.0, 5.0, 12.0, 8.0, 8.0]

    raizes, iters = bairstow(coeffs_apc2, r0=0.5, s0=0.5, tol=1e-10)

    # Validação cruzada com numpy
    raizes_np = np.roots(coeffs_apc2)

    print("\n" + "=" * 70)
    print("  ANÁLISE 2 — Sistema dinâmico 2 GDL (APC2)")
    print("  P(λ) = 2λ⁴ + 5λ³ + 12λ² + 8λ + 8")
    print("=" * 70)
    print(f"\n  Iterações por deflação: {iters}")
    print()
    print(f"{'k':>4}  {'λₖ':>28}  {'σₖ = Re(λ)':>12}  "
          f"{'ωd = |Im(λ)|':>14}  {'|λₖ|':>8}  {'Status':>10}")
    print("-" * 82)

    for k, lam in enumerate(raizes):
        sigma  = lam.real
        omega_d = abs(lam.imag)
        modulo = abs(lam)
        if sigma < -1e-8:
            status = "Estável"
        elif abs(sigma) < 1e-8:
            status = "Marginal"
        else:
            status = "INSTÁVEL"
        print(f"{k+1:>4}  {str(np.round(lam, 5)):>28}  {sigma:>12.5f}  "
              f"{omega_d:>14.5f}  {modulo:>8.4f}  {status:>10}")

    print(f"\n  Verificação numpy.roots: {np.round(np.sort_complex(raizes_np), 5)}")

    # --- Gráfico no plano complexo ---
    fig, ax = plt.subplots(figsize=(7, 6))

    cores = ["#d7191c", "#d7191c", "#2c7bb6", "#2c7bb6"]
    for k, lam in enumerate(raizes):
        lbl = f"$\\lambda_{k+1}$ = {np.round(lam, 4)}"
        ax.scatter(lam.real, lam.imag, s=140, color=cores[k],
                   edgecolors="k", linewidths=0.8, zorder=5, label=lbl)

    # Semiplano estável sombreado
    ylims = (-3.5, 3.5)
    ax.fill_betweenx(ylims, -6, 0, alpha=0.06, color="green", label="Semiplano estável")
    ax.set_ylim(ylims)
    ax.axvline(0, color="k", lw=1.2, ls="--", label="Eixo imaginário")
    ax.set_xlabel(r"Re($\lambda$)", fontsize=13)
    ax.set_ylabel(r"Im($\lambda$)", fontsize=13)
    ax.set_title("Autovalores do sistema 2-GDL (APC2)\n"
                 r"$P(\lambda) = 2\lambda^4 + 5\lambda^3 + 12\lambda^2 + 8\lambda + 8$",
                 fontsize=11)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if salvar:
        os.makedirs("resultados", exist_ok=True)
        fig.savefig("resultados/analise_2_apc2.png", dpi=150)
        print("\n[✓] resultados/analise_2_apc2.png salvo.")
    plt.show()

    return raizes


# ---------------------------------------------------------------------------
# 7. ANÁLISE 3 — Convergência para diferentes chutes iniciais (r₀, s₀)
# ---------------------------------------------------------------------------

def analise_convergencia(coeffs=None, salvar=True):
    """
    Testa o método de Bairstow com uma grade de chutes iniciais (r₀, s₀) e
    registra o número total de iterações e se houve convergência.

    Permite identificar regiões sensíveis no espaço de parâmetros iniciais
    antes da geração completa do fractal.
    """
    if coeffs is None:
        coeffs = [2.0, 5.0, 12.0, 8.0, 8.0]

    r0_vals = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
    s0_vals = [-2.0, -1.0,  0.0, 0.5, 1.0, 2.0]

    print("\n" + "=" * 60)
    print("  ANÁLISE 3 — Sensibilidade ao chute inicial (r₀, s₀)")
    print("=" * 60)
    print(f"\n{'r₀':>7}  {'s₀':>7}  {'Iters':>7}  {'Convergiu?':>11}  {'|λ₁| Bairstow':>15}")
    print("-" * 54)

    for r0 in r0_vals:
        for s0 in s0_vals:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                raizes, iters = bairstow(coeffs, r0=r0, s0=s0,
                                         tol=1e-10, max_iter=300)
                conv = "Sim" if len(w) == 0 else "Não"
            total = sum(iters)
            lam1 = abs(raizes[0])
            print(f"{r0:>7.2f}  {s0:>7.2f}  {total:>7}  {conv:>11}  {lam1:>15.6f}")


# ---------------------------------------------------------------------------
# 8. ANÁLISE 4 — Fractal de Bairstow
# ---------------------------------------------------------------------------

def _iters_fractal(coeffs, r0, s0, tol=1e-8, max_iter=100):
    """
    Conta as iterações para convergência do PRIMEIRO fator quadrático
    a partir de (r0, s0). Retorna max_iter se não converge ou diverge.
    Função auxiliar para geração do fractal.
    """
    poly = np.array(coeffs, dtype=float)
    poly = poly / poly[0]
    r, s = r0, s0

    for it in range(max_iter):
        # Verificar divergência
        if abs(r) > 1e8 or abs(s) > 1e8:
            return max_iter

        try:
            dr, ds, _, _ = _passo_bairstow(poly, r, s)
        except Exception:
            return max_iter

        if dr is None:
            return max_iter

        r += dr
        s += ds

        crit_r = abs(dr / r) if abs(r) > 1e-14 else abs(dr)
        crit_s = abs(ds / s) if abs(s) > 1e-14 else abs(ds)

        if crit_r < tol and crit_s < tol:
            return it + 1

    return max_iter


def analise_fractal(coeffs=None, r_range=(-3.0, 3.0), s_range=(-3.0, 3.0),
                    resolucao=350, max_iter=80, salvar=True):
    """
    Gera o fractal de Bairstow: um mapa de cor no plano (r₀, s₀) em que
    cada pixel é colorido pelo número de iterações necessárias para que o
    primeiro fator quadrático convirja.

    Pontos que não convergem dentro de max_iter são pintados de preto,
    revelando a geometria fractal das bacias de convergência do método.

    Parâmetros
    ----------
    coeffs     : array-like   — coeficientes do polinômio
    r_range    : (float,float) — intervalo de r₀ no eixo x
    s_range    : (float,float) — intervalo de s₀ no eixo y
    resolucao  : int           — nº de pontos por eixo (resolucao × resolucao)
    max_iter   : int           — máximo de iterações (define "divergência")
    salvar     : bool          — se True, salva a figura em resultados/
    """
    if coeffs is None:
        coeffs = [2.0, 5.0, 12.0, 8.0, 8.0]

    r_vals = np.linspace(r_range[0], r_range[1], resolucao)
    s_vals = np.linspace(s_range[0], s_range[1], resolucao)
    mapa   = np.zeros((resolucao, resolucao), dtype=int)

    print(f"\n[Fractal] Calculando grade {resolucao}×{resolucao} = "
          f"{resolucao**2:,} pontos...")

    for i, s0 in enumerate(s_vals):
        if (i + 1) % (resolucao // 10) == 0:
            print(f"  {100*(i+1)//resolucao:3d}% concluído...", flush=True)
        for j, r0 in enumerate(r_vals):
            mapa[i, j] = _iters_fractal(coeffs, r0, s0,
                                        tol=1e-8, max_iter=max_iter)

    print("  100% — pronto.")

    # --- Plotar ---
    fig, ax = plt.subplots(figsize=(8, 7))

    mapa_plot = mapa.astype(float)
    mapa_plot[mapa == max_iter] = np.nan   # divergência → NaN → preto

    cmap = plt.get_cmap("inferno").copy()
    cmap.set_bad("black")

    im = ax.imshow(
        mapa_plot,
        extent=[r_range[0], r_range[1], s_range[0], s_range[1]],
        origin="lower", cmap=cmap, aspect="equal",
        vmin=1, vmax=max_iter - 1
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Nº de iterações para convergência", fontsize=11)

    ax.axhline(0, color="white", lw=0.6, ls="--", alpha=0.4)
    ax.axvline(0, color="white", lw=0.6, ls="--", alpha=0.4)
    ax.set_xlabel(r"$r_0$", fontsize=14)
    ax.set_ylabel(r"$s_0$", fontsize=14)
    ax.set_title(
        "Fractal de Bairstow\n"
        r"$P(\lambda)=2\lambda^4+5\lambda^3+12\lambda^2+8\lambda+8$"
        f"\n(preto = não convergiu em {max_iter} iterações)",
        fontsize=11
    )
    fig.tight_layout()

    if salvar:
        os.makedirs("resultados", exist_ok=True)
        fig.savefig("resultados/analise_4_fractal.png", dpi=150)
        print("[✓] resultados/analise_4_fractal.png salvo.")
    plt.show()

    return mapa


# ---------------------------------------------------------------------------
# 9. EXPORTAR RAÍZES PARA CSV
# ---------------------------------------------------------------------------

def exportar_raizes(coeffs, r0=0.5, s0=0.5, tol=1e-10,
                    nome="resultados/raizes_apc2.csv"):
    """
    Calcula e salva as raízes do polinômio em arquivo CSV com colunas:
        indice, Re(λ), Im(λ), |λ|, ∠λ [rad], estabilidade
    """
    raizes, _ = bairstow(coeffs, r0=r0, s0=s0, tol=tol)
    os.makedirs("resultados", exist_ok=True)

    with open(nome, "w") as fout:
        fout.write("indice,Re_lambda,Im_lambda,modulo,fase_rad,estabilidade\n")
        for k, lam in enumerate(raizes):
            sigma = lam.real
            if sigma < -1e-8:
                status = "Estavel"
            elif abs(sigma) < 1e-8:
                status = "Marginal"
            else:
                status = "Instavel"
            fout.write(
                f"{k+1},{lam.real:.10f},{lam.imag:.10f},"
                f"{abs(lam):.10f},{np.angle(lam):.10f},{status}\n"
            )
    print(f"[✓] Raízes exportadas: {nome}")


# ---------------------------------------------------------------------------
# 10. PONTO DE ENTRADA
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("  PPC2 — Método de Bairstow — Raízes de Polinômios")
    print("  Cálculo Numérico Aplicado — Prof. Dr. Rafael Gabler Gontijo")
    print("=" * 60)

    print("\n[1/4] Validação com polinômio de grau 7 (raízes conhecidas)...")
    analise_validacao()

    print("\n[2/4] Polinômio característico do sistema 2-GDL (APC2)...")
    analise_apc2()

    print("\n[3/4] Análise de convergência para diferentes (r₀, s₀)...")
    analise_convergencia()

    print("\n[4/4] Gerando fractal de Bairstow (pode levar alguns minutos)...")
    analise_fractal(resolucao=350, max_iter=80)

    print("\n[+] Exportando raízes do APC2 em CSV...")
    exportar_raizes([2.0, 5.0, 12.0, 8.0, 8.0])

    print("\nConcluído. Todos os resultados foram salvos em ./resultados/")
