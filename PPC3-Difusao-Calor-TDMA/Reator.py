"""
PPC3 - Calculo Numerico Aplicado
Simulacao de difusao de calor transiente 1D em pastilha de UO2.

Modelo:
    rho*Cp*dT/dt = k*d2T/dx2 + q_dot

Condicoes de contorno:
    x = 0: simetria termica, dT/dx = 0
    x = L: conveccao, -k*dT/dx = h*(T_s - T_inf)

Metodo numerico:
    - Diferencas finitas no espaco
    - Esquema implicito no tempo
    - Algoritmo de Thomas para sistema tridiagonal

Autor: Mateus Leal Silva
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np


def thomas_algorithm(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """
    Resolve um sistema tridiagonal Ax = d pelo Algoritmo de Thomas.

    Parametros
    ----------
    a : diagonal inferior, com a[0] = 0
    b : diagonal principal
    c : diagonal superior, com c[-1] = 0
    d : vetor de termos independentes

    Retorno
    -------
    x : vetor solucao
    """
    n = len(d)
    cp = np.zeros(n, dtype=float)
    dp = np.zeros(n, dtype=float)
    x = np.zeros(n, dtype=float)

    if abs(b[0]) < 1.0e-15:
        raise ZeroDivisionError("Pivo nulo na primeira linha do sistema tridiagonal.")

    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]

    for i in range(1, n):
        den = b[i] - a[i] * cp[i - 1]
        if abs(den) < 1.0e-15:
            raise ZeroDivisionError(f"Pivo nulo na linha {i} do sistema tridiagonal.")

        if i < n - 1:
            cp[i] = c[i] / den
        dp[i] = (d[i] - a[i] * dp[i - 1]) / den

    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]

    return x


def montar_sistema_implicito(
    n: int,
    fo: float,
    bi_malha: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Monta as diagonais do sistema tridiagonal implicito.

    Para os nos internos:
        -Fo*T_{i-1}^{p+1} + (1+2Fo)*T_i^{p+1} - Fo*T_{i+1}^{p+1} = d_i

    No centro, usa-se simetria:
        (1+2Fo)*T_0^{p+1} - 2Fo*T_1^{p+1} = d_0

    Na superficie, usa-se conveccao:
        -2Fo*T_{N-2}^{p+1} + (1+2Fo+2Fo*Bi_m)*T_{N-1}^{p+1} = d_{N-1}
    """
    a = np.zeros(n, dtype=float)
    b = np.zeros(n, dtype=float)
    c = np.zeros(n, dtype=float)

    b[0] = 1.0 + 2.0 * fo
    c[0] = -2.0 * fo

    for i in range(1, n - 1):
        a[i] = -fo
        b[i] = 1.0 + 2.0 * fo
        c[i] = -fo

    a[-1] = -2.0 * fo
    b[-1] = 1.0 + 2.0 * fo + 2.0 * fo * bi_malha

    return a, b, c


def simular_reator_implicito(
    *,
    L: float,
    k: float,
    rho: float,
    cp: float,
    h: float,
    T_inf: float,
    T_inicial: float,
    q_dot: float,
    n: int,
    dt: float,
    t_final: float,
    tempos_saida: Iterable[float],
) -> tuple[np.ndarray, list[float], list[np.ndarray]]:
    """
    Simula a conducao transiente 1D com geracao volumetrica e conveccao.

    Retorna:
        x: coordenadas espaciais
        tempos_salvos: tempos efetivamente armazenados
        historico: lista dos perfis de temperatura
    """
    if n < 3:
        raise ValueError("Use pelo menos 3 nos espaciais.")
    if dt <= 0.0 or t_final <= 0.0:
        raise ValueError("dt e t_final devem ser positivos.")

    alpha = k / (rho * cp)
    dx = L / (n - 1)
    x = np.linspace(0.0, L, n)

    fo = alpha * dt / dx**2
    bi_malha = h * dx / k
    termo_geracao = q_dot * dt / (rho * cp)

    a, b, c = montar_sistema_implicito(n, fo, bi_malha)

    T = np.full(n, T_inicial, dtype=float)
    tempos_saida = sorted(set(float(t) for t in tempos_saida if 0.0 <= float(t) <= t_final))

    tempos_salvos = [0.0]
    historico = [T.copy()]

    indice_saida = 0
    while indice_saida < len(tempos_saida) and tempos_saida[indice_saida] <= 0.0:
        indice_saida += 1

    tempo = 0.0
    n_passos = int(np.ceil(t_final / dt))

    for _ in range(n_passos):
        passo = min(dt, t_final - tempo)
        if passo <= 0.0:
            break

        # Remonta se o ultimo passo for menor que dt.
        if abs(passo - dt) > 1.0e-15:
            fo = alpha * passo / dx**2
            termo_geracao = q_dot * passo / (rho * cp)
            a, b, c = montar_sistema_implicito(n, fo, bi_malha)

        d = T + termo_geracao
        d[-1] += 2.0 * fo * bi_malha * T_inf

        T = thomas_algorithm(a, b, c, d)
        tempo += passo

        while indice_saida < len(tempos_saida) and tempo >= tempos_saida[indice_saida] - 1.0e-12:
            tempos_salvos.append(tempo)
            historico.append(T.copy())
            indice_saida += 1

    if tempos_salvos[-1] < t_final - 1.0e-12:
        tempos_salvos.append(tempo)
        historico.append(T.copy())

    return x, tempos_salvos, historico


def perfil_permanente_com_geracao(
    x: np.ndarray,
    *,
    L: float,
    k: float,
    h: float,
    T_inf: float,
    q_dot: float,
) -> np.ndarray:
    """
    Solucao analitica em regime permanente para placa plana com simetria em x=0
    e conveccao em x=L.

    T(x) = T_inf + q_dot*L/h + q_dot*(L^2 - x^2)/(2k)
    """
    return T_inf + q_dot * L / h + q_dot * (L**2 - x**2) / (2.0 * k)


def eigenvalues_plane_wall(Bi: float, n_roots: int = 30) -> np.ndarray:
    """
    Calcula raizes positivas de lambda*tan(lambda)=Bi por bissecao.
    Usado para a solucao analitica transiente sem geracao.
    """
    roots = []
    eps = 1.0e-10

    for m in range(n_roots):
        left = m * np.pi + eps
        right = m * np.pi + np.pi / 2.0 - eps

        def f(lam: float) -> float:
            return lam * np.tan(lam) - Bi

        f_left = f(left)
        for _ in range(100):
            mid = 0.5 * (left + right)
            f_mid = f(mid)

            if f_left * f_mid <= 0.0:
                right = mid
            else:
                left = mid
                f_left = f_mid

        roots.append(0.5 * (left + right))

    return np.array(roots)


def solucao_exata_sem_geracao(
    x: np.ndarray,
    t: float,
    *,
    L: float,
    k: float,
    rho: float,
    cp: float,
    h: float,
    T_inf: float,
    T_inicial: float,
    n_terms: int = 30,
) -> np.ndarray:
    """
    Solucao por serie para parede plana com simetria no centro e conveccao
    na superficie, sem geracao interna.
    """
    if t <= 0.0:
        return np.full_like(x, T_inicial, dtype=float)

    alpha = k / (rho * cp)
    Bi = h * L / k
    Fo = alpha * t / L**2
    x_star = x / L

    theta = np.zeros_like(x_star, dtype=float)
    lambdas = eigenvalues_plane_wall(Bi, n_terms)

    for lam in lambdas:
        coef = 4.0 * np.sin(lam) / (2.0 * lam + np.sin(2.0 * lam))
        theta += coef * np.cos(lam * x_star) * np.exp(-(lam**2) * Fo)

    return T_inf + theta * (T_inicial - T_inf)


def salvar_tabela_perfis(
    path: Path,
    x: np.ndarray,
    tempos: list[float],
    historico: list[np.ndarray],
) -> None:
    """Salva perfis transientes em arquivo .dat."""
    matriz = np.column_stack([x] + historico)
    header = "x_m " + " ".join(f"T_t_{t:.6g}_s" for t in tempos)
    np.savetxt(path, matriz, header=header, fmt="%.12e")


def salvar_validacao(
    path: Path,
    x: np.ndarray,
    T_num: np.ndarray,
    T_exata: np.ndarray,
) -> None:
    """Salva comparacao numerica versus analitica sem geracao."""
    erro = T_num - T_exata
    matriz = np.column_stack([x, T_num, T_exata, erro])
    np.savetxt(path, matriz, header="x_m T_numerica_C T_exata_C erro_C", fmt="%.12e")


def tentar_plotar(
    outdir: Path,
    x: np.ndarray,
    tempos: list[float],
    historico: list[np.ndarray],
    T_perm: np.ndarray,
    x_val: np.ndarray,
    T_num_val: np.ndarray,
    T_exata_val: np.ndarray,
) -> None:
    """Gera graficos se matplotlib estiver instalado."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib nao encontrado; os arquivos .dat foram gerados normalmente.")
        return

    plt.figure()
    for tempo, T in zip(tempos, historico):
        plt.plot(x * 1000.0, T, label=f"t = {tempo:.2f} s")
    plt.plot(x * 1000.0, T_perm, "--", label="Regime permanente")
    plt.xlabel("x [mm]")
    plt.ylabel("Temperatura [C]")
    plt.title("PPC3 - Perfis transientes com geracao interna")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "ppc3_perfis_transientes.png", dpi=200)
    plt.close()

    plt.figure()
    plt.plot(x_val * 1000.0, T_num_val, "o", label="Numerico")
    plt.plot(x_val * 1000.0, T_exata_val, "-", label="Exato")
    plt.xlabel("x [mm]")
    plt.ylabel("Temperatura [C]")
    plt.title("PPC3 - Validacao sem geracao interna")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "ppc3_validacao.png", dpi=200)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="PPC3 - Difusao de calor transiente 1D com esquema implicito e TDMA."
    )
    parser.add_argument("--outdir", default=".", help="Pasta de saida dos arquivos gerados.")
    parser.add_argument("--N", type=int, default=51, help="Numero de nos espaciais.")
    parser.add_argument("--dt", type=float, default=0.01, help="Passo de tempo [s].")
    parser.add_argument("--t-final", type=float, default=20.0, help="Tempo final da simulacao [s].")
    parser.add_argument("--sem-plots", action="store_true", help="Nao gera imagens PNG.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parametros fisicos adotados no relatorio PPC3.
    L = 5.0e-3          # m
    k = 4.0             # W/(m K)
    rho = 10500.0       # kg/m^3
    cp = 300.0          # J/(kg K)
    h = 30000.0         # W/(m^2 K)
    T_inf = 300.0       # graus Celsius
    q_dot = 3.0e8       # W/m^3

    # Caso principal: aquecimento com geracao interna.
    tempos_saida = [0.0, 1.0, 2.0, 5.0, 10.0, 12.0, args.t_final]
    x, tempos, historico = simular_reator_implicito(
        L=L,
        k=k,
        rho=rho,
        cp=cp,
        h=h,
        T_inf=T_inf,
        T_inicial=T_inf,
        q_dot=q_dot,
        n=args.N,
        dt=args.dt,
        t_final=args.t_final,
        tempos_saida=tempos_saida,
    )

    T_perm = perfil_permanente_com_geracao(x, L=L, k=k, h=h, T_inf=T_inf, q_dot=q_dot)
    salvar_tabela_perfis(outdir / "ppc3_perfis_transientes.dat", x, tempos, historico)

    # Caso de validacao: sem geracao interna, comparado com a solucao em serie.
    t_validacao = 10.0
    x_val, tempos_val, hist_val = simular_reator_implicito(
        L=L,
        k=k,
        rho=rho,
        cp=cp,
        h=h,
        T_inf=T_inf,
        T_inicial=500.0,
        q_dot=0.0,
        n=args.N,
        dt=args.dt,
        t_final=t_validacao,
        tempos_saida=[t_validacao],
    )
    T_num_val = hist_val[-1]
    T_exata_val = solucao_exata_sem_geracao(
        x_val,
        t_validacao,
        L=L,
        k=k,
        rho=rho,
        cp=cp,
        h=h,
        T_inf=T_inf,
        T_inicial=500.0,
    )
    salvar_validacao(outdir / "ppc3_validacao_sem_geracao.dat", x_val, T_num_val, T_exata_val)

    erro_max = float(np.max(np.abs(T_num_val - T_exata_val)))
    T_final = historico[-1]

    if not args.sem_plots:
        tentar_plotar(outdir, x, tempos, historico, T_perm, x_val, T_num_val, T_exata_val)

    print("PPC3 - simulacao concluida.")
    print(f"Arquivos gerados em: {outdir.resolve()}")
    print(f"Temperatura maxima no tempo final: {np.max(T_final):.6f} C")
    print(f"Temperatura maxima em regime permanente: {np.max(T_perm):.6f} C")
    print(f"Erro maximo na validacao sem geracao: {erro_max:.6e} C")


if __name__ == "__main__":
    main()
