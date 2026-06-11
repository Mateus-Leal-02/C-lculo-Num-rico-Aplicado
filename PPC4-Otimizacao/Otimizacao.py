"""
PPC4 - Calculo Numerico Aplicado
Otimizacao bidimensional sem restricoes.

Funcao objetivo:
    f(x,y) = 2xy + 2x - x^2 - 2y^2

Metodos implementados:
    1) Aclive maximo
    2) Gradientes conjugados de Fletcher-Reeves

Arquivos gerados:
    output1.dat  -> historico do aclive maximo
    output2.dat  -> historico do Fletcher-Reeves
    function.dat -> malha x, y, f(x,y) para curvas de nivel

Autor: Mateus Leal Silva
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


HESSIANA = np.array(
    [
        [-2.0,  2.0],
        [ 2.0, -4.0],
    ],
    dtype=float,
)


def funcao(z: np.ndarray) -> float:
    """Calcula f(x,y) = 2xy + 2x - x^2 - 2y^2."""
    x, y = z
    return float(2.0 * x * y + 2.0 * x - x**2 - 2.0 * y**2)


def gradiente(z: np.ndarray) -> np.ndarray:
    """Calcula o gradiente de f."""
    x, y = z
    return np.array(
        [
            2.0 * y + 2.0 - 2.0 * x,
            2.0 * x - 4.0 * y,
        ],
        dtype=float,
    )


def norma(v: np.ndarray) -> float:
    """Norma Euclidiana."""
    return float(np.sqrt(np.dot(v, v)))


def passo_otimo_quadratico(z: np.ndarray, direcao: np.ndarray, eps: float = 1.0e-15) -> float:
    """
    Calcula o passo que maximiza f(z + h*d) para uma funcao quadratica.

    Como f e quadratica:
        phi'(h) = grad(f(z))^T d + h d^T H d

    No ponto critico:
        h = - grad(f(z))^T d / (d^T H d)
    """
    g = gradiente(z)
    denominador = float(direcao @ HESSIANA @ direcao)

    if abs(denominador) < eps:
        return 0.0

    return -float(g @ direcao) / denominador


def aclive_maximo(
    z0: np.ndarray,
    tolerancia: float,
    max_iter: int,
) -> tuple[np.ndarray, list[dict[str, float]]]:
    """Executa o metodo do aclive maximo."""
    z = np.array(z0, dtype=float)
    historico = []

    for k in range(max_iter + 1):
        g = gradiente(z)
        erro = norma(g)

        if erro < tolerancia:
            historico.append(
                {"iter": k, "erro": erro, "h": 0.0, "x": z[0], "y": z[1], "dfdx": g[0], "dfdy": g[1], "f": funcao(z)}
            )
            break

        direcao = g.copy()
        h = passo_otimo_quadratico(z, direcao)

        historico.append(
            {"iter": k, "erro": erro, "h": h, "x": z[0], "y": z[1], "dfdx": g[0], "dfdy": g[1], "f": funcao(z)}
        )

        z = z + h * direcao

    return z, historico


def fletcher_reeves(
    z0: np.ndarray,
    tolerancia: float,
    max_iter: int,
) -> tuple[np.ndarray, list[dict[str, float]]]:
    """
    Executa o metodo de gradientes conjugados de Fletcher-Reeves para maximizacao.

    Para maximizacao de uma funcao concava, a primeira direcao e o proprio gradiente.
    """
    z = np.array(z0, dtype=float)
    g = gradiente(z)
    direcao = g.copy()
    historico = []

    for k in range(max_iter + 1):
        erro = norma(g)

        if erro < tolerancia:
            historico.append(
                {"iter": k, "erro": erro, "h": 0.0, "x": z[0], "y": z[1], "dfdx": g[0], "dfdy": g[1], "f": funcao(z)}
            )
            break

        h = passo_otimo_quadratico(z, direcao)
        historico.append(
            {"iter": k, "erro": erro, "h": h, "x": z[0], "y": z[1], "dfdx": g[0], "dfdy": g[1], "f": funcao(z)}
        )

        z_novo = z + h * direcao
        g_novo = gradiente(z_novo)

        denominador = float(g @ g)
        beta = float(g_novo @ g_novo / denominador) if denominador > 0.0 else 0.0

        direcao = g_novo + beta * direcao
        z = z_novo
        g = g_novo

    return z, historico


def salvar_historico(path: Path, historico: list[dict[str, float]]) -> None:
    """Salva o historico iterativo em arquivo .dat."""
    with path.open("w", encoding="utf-8") as f:
        f.write("# iter erro h x y dfdx dfdy f\n")
        for item in historico:
            f.write(
                f"{int(item['iter']):4d} "
                f"{item['erro']: .12e} "
                f"{item['h']: .12e} "
                f"{item['x']: .12e} "
                f"{item['y']: .12e} "
                f"{item['dfdx']: .12e} "
                f"{item['dfdy']: .12e} "
                f"{item['f']: .12e}\n"
            )


def gerar_function_dat(
    path: Path,
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
    nx: int,
    ny: int,
) -> None:
    """Gera arquivo x, y, f(x,y) para visualizacao de curvas de nivel."""
    xs = np.linspace(xmin, xmax, nx)
    ys = np.linspace(ymin, ymax, ny)

    with path.open("w", encoding="utf-8") as f:
        f.write("# x y f\n")
        for x in xs:
            for y in ys:
                z = np.array([x, y], dtype=float)
                f.write(f"{x: .12e} {y: .12e} {funcao(z): .12e}\n")
            f.write("\n")


def tentar_plotar(outdir: Path, hist1: list[dict[str, float]], hist2: list[dict[str, float]]) -> None:
    """Gera uma figura simples do caminho iterativo, se matplotlib estiver instalado."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib nao encontrado; os arquivos .dat foram gerados normalmente.")
        return

    x = np.linspace(-4.0, 4.0, 200)
    y = np.linspace(-3.0, 5.0, 200)
    X, Y = np.meshgrid(x, y)
    Z = 2.0 * X * Y + 2.0 * X - X**2 - 2.0 * Y**2

    p1 = np.array([[item["x"], item["y"]] for item in hist1])
    p2 = np.array([[item["x"], item["y"]] for item in hist2])

    plt.figure()
    plt.contour(X, Y, Z, levels=30)
    plt.plot(p1[:, 0], p1[:, 1], "o-", label="Aclive maximo")
    plt.plot(p2[:, 0], p2[:, 1], "s-", label="Fletcher-Reeves")
    plt.plot(2.0, 1.0, "*", markersize=12, label="Otimo analitico")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("PPC4 - Caminho de otimizacao")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "ppc4_caminhos.png", dpi=200)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="PPC4 - Otimizacao por aclive maximo e Fletcher-Reeves.")
    parser.add_argument("--x0", type=float, default=-2.0, help="Coordenada inicial x0.")
    parser.add_argument("--y0", type=float, default=3.0, help="Coordenada inicial y0.")
    parser.add_argument("--tol", type=float, default=1.0e-8, help="Tolerancia para ||grad f||.")
    parser.add_argument("--max-iter", type=int, default=100, help="Numero maximo de iteracoes.")
    parser.add_argument("--outdir", default=".", help="Pasta de saida dos arquivos .dat.")
    parser.add_argument("--sem-plots", action="store_true", help="Nao gera imagem PNG.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    z0 = np.array([args.x0, args.y0], dtype=float)

    sol_aclive, hist_aclive = aclive_maximo(z0, args.tol, args.max_iter)
    sol_fr, hist_fr = fletcher_reeves(z0, args.tol, args.max_iter)

    salvar_historico(outdir / "output1.dat", hist_aclive)
    salvar_historico(outdir / "output2.dat", hist_fr)
    gerar_function_dat(outdir / "function.dat", xmin=-4.0, xmax=4.0, ymin=-3.0, ymax=5.0, nx=101, ny=101)

    if not args.sem_plots:
        tentar_plotar(outdir, hist_aclive, hist_fr)

    print("PPC4 - otimizacao concluida.")
    print(f"Arquivos gerados em: {outdir.resolve()}")

    print("\nResultado por aclive maximo:")
    print(f"  x = {sol_aclive[0]:.12f}")
    print(f"  y = {sol_aclive[1]:.12f}")
    print(f"  f = {funcao(sol_aclive):.12f}")

    print("\nResultado por Fletcher-Reeves:")
    print(f"  x = {sol_fr[0]:.12f}")
    print(f"  y = {sol_fr[1]:.12f}")
    print(f"  f = {funcao(sol_fr):.12f}")

    print("\nOtimo analitico esperado:")
    print("  x* = 2.000000000000")
    print("  y* = 1.000000000000")
    print("  f* = 2.000000000000")


if __name__ == "__main__":
    main()
