"""
=============================================================================
  Sedimentação de Esfera em Baixo Reynolds — Runge-Kutta de 4ª Ordem (RK4)
=============================================================================

Descrição:
    Este programa resolve numericamente a equação adimensional do movimento
    de uma esfera sedimentando em um fluido viscoso (regime de baixo Reynolds),
    com e sem a correção de arrasto de Oseen.

    Equação adimensional (com correção de Oseen):
        St * dv*/dt* = 1 - v* - (3/8) * Re_s * v*²

    Para Re_s → 0 (arrasto de Stokes puro):
        St * dv*/dt* = 1 - v*
        Solução analítica: v*(t*) = 1 - exp(-t*/St)

Referências:
    - Chapra & Canale, Métodos Numéricos para Engenharia, 5ª ed., McGraw-Hill (2008)
    - Sobral et al., Powder Technology, 178 (2007), 129–141

Autor: [Seu Nome] | Matrícula: XXXXXXXX
Disciplina: Cálculo Numérico Aplicado — Prof. Dr. Rafael Gabler Gontijo
Data: março de 2026
=============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ---------------------------------------------------------------------------
# 1. FUNÇÃO QUE DEFINE O LADO DIREITO DO PVI
# ---------------------------------------------------------------------------

def f(t, v, St, Re_s):
    """
    Lado direito da equação adimensional do movimento:
        dv*/dt* = (1 - v* - (3/8)*Re_s*v*²) / St

    Parâmetros
    ----------
    t    : float — instante de tempo adimensional (não aparece explicitamente,
                   incluído para manter a assinatura genérica do RK4)
    v    : float — velocidade adimensional v* no instante t*
    St   : float — número de Stokes
    Re_s : float — número de Reynolds de partícula

    Retorna
    -------
    float — derivada dv*/dt* no ponto (t*, v*)
    """
    return (1.0 - v - (3.0 / 8.0) * Re_s * v**2) / St


# ---------------------------------------------------------------------------
# 2. INTEGRADOR RK4 GENÉRICO
# ---------------------------------------------------------------------------

def rk4(f, v0, t0, t_max, h, St, Re_s):
    """
    Método de Runge-Kutta de 4ª ordem para um PVI escalar.

    Relações de recorrência:
        k1 = f(t_i,          v_i)
        k2 = f(t_i + h/2,   v_i + h/2 * k1)
        k3 = f(t_i + h/2,   v_i + h/2 * k2)
        k4 = f(t_i + h,     v_i + h   * k3)
        v_{i+1} = v_i + (h/6) * (k1 + 2k2 + 2k3 + k4)

    Parâmetros
    ----------
    f      : callable — função f(t, v, St, Re_s)
    v0     : float    — condição inicial v*(0)
    t0     : float    — tempo inicial
    t_max  : float    — tempo final
    h      : float    — passo de tempo adimensional
    St     : float    — número de Stokes
    Re_s   : float    — número de Reynolds de partícula

    Retorna
    -------
    t_arr : np.ndarray — vetor de instantes de tempo
    v_arr : np.ndarray — vetor de velocidades numéricas
    """
    N = int((t_max - t0) / h)          # número de passos
    t_arr = np.zeros(N + 1)
    v_arr = np.zeros(N + 1)

    t_arr[0] = t0
    v_arr[0] = v0

    for i in range(N):
        ti = t_arr[i]
        vi = v_arr[i]

        k1 = f(ti,           vi,                   St, Re_s)
        k2 = f(ti + h / 2,   vi + h / 2 * k1,     St, Re_s)
        k3 = f(ti + h / 2,   vi + h / 2 * k2,     St, Re_s)
        k4 = f(ti + h,       vi + h * k3,          St, Re_s)

        v_arr[i + 1] = vi + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        t_arr[i + 1] = ti + h

    return t_arr, v_arr


# ---------------------------------------------------------------------------
# 3. SOLUÇÃO ANALÍTICA (LIMITE Re_s → 0)
# ---------------------------------------------------------------------------

def solucao_analitica_stokes(t, St):
    """
    Solução exata para Re_s = 0 (arrasto de Stokes puro):
        v*(t*) = 1 - exp(-t*/St)

    Parâmetros
    ----------
    t  : np.ndarray — vetor de tempos adimensionais
    St : float      — número de Stokes

    Retorna
    -------
    np.ndarray — velocidade analítica em cada instante
    """
    return 1.0 - np.exp(-t / St)


# ---------------------------------------------------------------------------
# 4. ANÁLISE 1 — Comparação numérico vs. analítico para vários St (Re_s = 0)
# ---------------------------------------------------------------------------

def analise_stokes(t_max=15.0, h=0.05, salvar=True):
    """
    Plota v*(t*) numérico (RK4) vs. solução analítica de Stokes
    para diferentes valores do número de Stokes (Re_s = 0).
    """
    valores_St = [0.5, 1.0, 2.0, 4.0]
    Re_s = 0.0
    cores = ["#2c7bb6", "#d7191c", "#1a9641", "#f07d02"]

    fig, ax = plt.subplots(figsize=(8, 5))

    for St, cor in zip(valores_St, cores):
        t_num, v_num = rk4(f, v0=0.0, t0=0.0, t_max=t_max, h=h, St=St, Re_s=Re_s)
        v_ex = solucao_analitica_stokes(t_num, St)

        ax.plot(t_num, v_num, "-",  color=cor, lw=2,   label=f"RK4  St={St}")
        ax.plot(t_num, v_ex,  "--", color=cor, lw=1.5, label=f"Exata St={St}")

    ax.set_xlabel(r"$t^*$", fontsize=13)
    ax.set_ylabel(r"$v^*(t^*)$", fontsize=13)
    ax.set_title(r"Solução numérica vs. analítica — $Re_s = 0$", fontsize=13)
    ax.legend(fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if salvar:
        os.makedirs("resultados", exist_ok=True)
        fig.savefig("resultados/analise_1_stokes.png", dpi=150)
        print("[✓] resultados/analise_1_stokes.png salvo.")
    plt.show()


# ---------------------------------------------------------------------------
# 5. ANÁLISE 2 — Convergência: efeito do passo de tempo h
# ---------------------------------------------------------------------------

def analise_convergencia(St=1.0, Re_s=0.0, t_max=10.0, salvar=True):
    """
    Varia o passo de tempo h e analisa o erro máximo em relação
    à solução analítica. Verifica a convergência de 4ª ordem do RK4.
    """
    passos = [1.0, 0.5, 0.1, 0.05, 0.01]
    erros  = []

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    for h in passos:
        t_num, v_num = rk4(f, 0.0, 0.0, t_max, h, St, Re_s)
        v_ex  = solucao_analitica_stokes(t_num, St)
        erro  = np.abs(v_num - v_ex)
        erros.append(erro.max())
        axes[0].plot(t_num, erro, label=f"h = {h}")

    axes[0].set_xlabel(r"$t^*$", fontsize=12)
    axes[0].set_ylabel(r"Erro absoluto $|v^*_{RK4} - v^*_{ex}|$", fontsize=11)
    axes[0].set_title("Erro ao longo do tempo", fontsize=12)
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    # Gráfico log-log: erro máximo vs. h
    axes[1].loglog(passos, erros, "o-", color="#2c7bb6", lw=2, ms=7)
    # Linha de referência de ordem 4
    h_ref = np.array([passos[0], passos[-1]])
    axes[1].loglog(h_ref, erros[0] * (h_ref / passos[0])**4,
                   "k--", label=r"$\mathcal{O}(h^4)$")
    axes[1].set_xlabel(r"Passo de tempo $h$", fontsize=12)
    axes[1].set_ylabel("Erro máximo", fontsize=12)
    axes[1].set_title("Convergência do RK4", fontsize=12)
    axes[1].legend(fontsize=10)
    axes[1].grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    if salvar:
        os.makedirs("resultados", exist_ok=True)
        fig.savefig("resultados/analise_2_convergencia.png", dpi=150)
        print("[✓] resultados/analise_2_convergencia.png salvo.")
    plt.show()


# ---------------------------------------------------------------------------
# 6. ANÁLISE 3 — Efeito do número de Reynolds (Re_s ≠ 0)
# ---------------------------------------------------------------------------

def analise_reynolds(St=1.0, t_max=15.0, h=0.05, salvar=True):
    """
    Resolve o PVI com diferentes valores de Re_s e compara com
    o limite assintótico Re_s → 0 (Stokes puro).
    """
    valores_Re = [0.0, 0.2, 0.5, 1.0]
    cores = ["#1a1a2e", "#e94560", "#0f3460", "#533483"]

    fig, ax = plt.subplots(figsize=(8, 5))

    for Re_s, cor in zip(valores_Re, cores):
        t_num, v_num = rk4(f, 0.0, 0.0, t_max, h, St, Re_s)
        label = r"$Re_s = $" + f"{Re_s}" if Re_s > 0 else r"$Re_s \to 0$ (Stokes)"
        ax.plot(t_num, v_num, lw=2, color=cor, label=label)

    # Velocidade terminal analítica de Stokes (linha de referência)
    ax.axhline(1.0, color="gray", lw=1, ls=":", label=r"$v^*_\infty = 1$ (Stokes)")

    ax.set_xlabel(r"$t^*$", fontsize=13)
    ax.set_ylabel(r"$v^*(t^*)$", fontsize=13)
    ax.set_title(r"Efeito do $Re_s$ na dinâmica de sedimentação — $St = $" + f"{St}",
                 fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if salvar:
        os.makedirs("resultados", exist_ok=True)
        fig.savefig("resultados/analise_3_reynolds.png", dpi=150)
        print("[✓] resultados/analise_3_reynolds.png salvo.")
    plt.show()


# ---------------------------------------------------------------------------
# 7. EXPORTAR DADOS PARA ARQUIVO
# ---------------------------------------------------------------------------

def exportar_dados(St=1.0, Re_s=0.0, h=0.05, t_max=15.0):
    """
    Salva os resultados numéricos e analíticos em arquivo .csv.
    """
    t_num, v_num = rk4(f, 0.0, 0.0, t_max, h, St, Re_s)
    v_ex = solucao_analitica_stokes(t_num, St) if Re_s == 0.0 else np.full_like(t_num, np.nan)

    os.makedirs("resultados", exist_ok=True)
    nome = f"resultados/dados_St{St}_Re{Re_s}_h{h}.csv"
    header = "t_adim,v_num,v_analitica"
    np.savetxt(nome, np.column_stack([t_num, v_num, v_ex]),
               delimiter=",", header=header, comments="")
    print(f"[✓] Dados exportados: {nome}")


# ---------------------------------------------------------------------------
# 8. PONTO DE ENTRADA
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("  PPC1 — Sedimentação de Esfera — RK4")
    print("=" * 60)

    print("\n[1/3] Comparação com solução analítica (Re_s = 0)...")
    analise_stokes(t_max=15.0, h=0.05)

    print("\n[2/3] Análise de convergência em h...")
    analise_convergencia(St=1.0, Re_s=0.0, t_max=10.0)

    print("\n[3/3] Efeito do número de Reynolds...")
    analise_reynolds(St=1.0, t_max=15.0, h=0.05)

    print("\n[+] Exportando CSV de exemplo...")
    exportar_dados(St=1.0, Re_s=0.0, h=0.05, t_max=15.0)

    print("\nConcluído. Resultados em ./resultados/")
