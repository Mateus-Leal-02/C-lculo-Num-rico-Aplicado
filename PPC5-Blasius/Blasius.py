#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PPC5 — Solução numérica da equação de Blasius.

O programa resolve o problema de valor de contorno

    f''' + (1/2) f f'' = 0

com

    f(0) = 0,  f'(0) = 0,  f'(infinito) = 1,

por meio de:

1. Método do Tiro para determinar o valor desconhecido s = f''(0);
2. Runge–Kutta clássico de quarta ordem (RK4) para integrar o sistema
   equivalente de três EDOs de primeira ordem;
3. Newton–Raphson com derivada numérica para atualizar o parâmetro de tiro.

A integração RK4 é implementada diretamente neste arquivo. Não são utilizadas
bibliotecas que realizem a integração numérica pronta.

Autor: Mateus Leal Silva
Disciplina: Cálculo Numérico Aplicado
"""

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


# =============================================================================
# Estruturas de dados
# =============================================================================

@dataclass
class BlasiusSolution:
    """Armazena a solução integrada para um chute específico de f''(0)."""

    eta: list[float]
    f: list[float]
    fp: list[float]
    fpp: list[float]


@dataclass
class ShootingIteration:
    """Dados registrados em uma iteração do Método do Tiro."""

    iteration: int
    s: float
    fp_eta_max: float
    error: float


@dataclass
class ShootingResult:
    """Resultado completo do processo iterativo do Método do Tiro."""

    s: float
    error: float
    iterations: int
    converged: bool
    solution: BlasiusSolution
    history: list[ShootingIteration]


# =============================================================================
# Equação de Blasius e integração por RK4
# =============================================================================

def blasius_rhs(f: float, fp: float, fpp: float) -> tuple[float, float, float]:
    """Calcula o lado direito do sistema equivalente de primeira ordem.

    São definidas as variáveis auxiliares:

        y1 = f
        y2 = f'
        y3 = f''

    Assim, a equação de Blasius é reescrita como:

        y1' = y2
        y2' = y3
        y3' = -(1/2) y1 y3

    Parameters
    ----------
    f, fp, fpp:
        Valores atuais de f, f' e f''.

    Returns
    -------
    tuple[float, float, float]
        Derivadas de f, f' e f'' em relação a eta.
    """
    df_deta = fp
    dfp_deta = fpp
    dfpp_deta = -0.5 * f * fpp
    return df_deta, dfp_deta, dfpp_deta


def rk4_step(
    f: float,
    fp: float,
    fpp: float,
    step: float,
) -> tuple[float, float, float]:
    """Avança simultaneamente f, f' e f'' por um passo usando RK4.

    Os 12 coeficientes do método são escritos explicitamente:
    quatro coeficientes para cada uma das três variáveis do sistema.

    Parameters
    ----------
    f, fp, fpp:
        Estado no início do passo.
    step:
        Tamanho do passo Delta eta.

    Returns
    -------
    tuple[float, float, float]
        Estado atualizado ao final do passo.
    """

    # -------------------------------------------------------------------------
    # Estágio 1: avaliação no início do intervalo
    # -------------------------------------------------------------------------
    rhs1_f, rhs1_fp, rhs1_fpp = blasius_rhs(f, fp, fpp)

    k1_f = step * rhs1_f
    k1_fp = step * rhs1_fp
    k1_fpp = step * rhs1_fpp

    # -------------------------------------------------------------------------
    # Estágio 2: avaliação no ponto médio usando os coeficientes k1
    # -------------------------------------------------------------------------
    f_stage2 = f + 0.5 * k1_f
    fp_stage2 = fp + 0.5 * k1_fp
    fpp_stage2 = fpp + 0.5 * k1_fpp

    rhs2_f, rhs2_fp, rhs2_fpp = blasius_rhs(
        f_stage2,
        fp_stage2,
        fpp_stage2,
    )

    k2_f = step * rhs2_f
    k2_fp = step * rhs2_fp
    k2_fpp = step * rhs2_fpp

    # -------------------------------------------------------------------------
    # Estágio 3: nova avaliação no ponto médio usando os coeficientes k2
    # -------------------------------------------------------------------------
    f_stage3 = f + 0.5 * k2_f
    fp_stage3 = fp + 0.5 * k2_fp
    fpp_stage3 = fpp + 0.5 * k2_fpp

    rhs3_f, rhs3_fp, rhs3_fpp = blasius_rhs(
        f_stage3,
        fp_stage3,
        fpp_stage3,
    )

    k3_f = step * rhs3_f
    k3_fp = step * rhs3_fp
    k3_fpp = step * rhs3_fpp

    # -------------------------------------------------------------------------
    # Estágio 4: avaliação no final do intervalo usando os coeficientes k3
    # -------------------------------------------------------------------------
    f_stage4 = f + k3_f
    fp_stage4 = fp + k3_fp
    fpp_stage4 = fpp + k3_fpp

    rhs4_f, rhs4_fp, rhs4_fpp = blasius_rhs(
        f_stage4,
        fp_stage4,
        fpp_stage4,
    )

    k4_f = step * rhs4_f
    k4_fp = step * rhs4_fp
    k4_fpp = step * rhs4_fpp

    # -------------------------------------------------------------------------
    # Combinação ponderada dos quatro estágios do RK4
    # -------------------------------------------------------------------------
    f_next = f + (k1_f + 2.0 * k2_f + 2.0 * k3_f + k4_f) / 6.0
    fp_next = fp + (k1_fp + 2.0 * k2_fp + 2.0 * k3_fp + k4_fp) / 6.0
    fpp_next = fpp + (
        k1_fpp + 2.0 * k2_fpp + 2.0 * k3_fpp + k4_fpp
    ) / 6.0

    return f_next, fp_next, fpp_next


def integrate_blasius(
    s: float,
    delta_eta: float,
    eta_max: float,
) -> BlasiusSolution:
    """Integra o PVI equivalente para um chute s = f''(0).

    Condições iniciais:

        f(0)   = 0
        f'(0)  = 0
        f''(0) = s

    O último passo é reduzido quando necessário para que a integração termine
    exatamente em eta_max.
    """
    eta_values = [0.0]
    f_values = [0.0]
    fp_values = [0.0]
    fpp_values = [float(s)]

    eta_current = 0.0

    while eta_current < eta_max - 1.0e-15:
        # Evita ultrapassar eta_max caso o intervalo não seja múltiplo do passo.
        step = min(delta_eta, eta_max - eta_current)

        f_next, fp_next, fpp_next = rk4_step(
            f_values[-1],
            fp_values[-1],
            fpp_values[-1],
            step,
        )

        eta_current += step

        eta_values.append(eta_current)
        f_values.append(f_next)
        fp_values.append(fp_next)
        fpp_values.append(fpp_next)

    return BlasiusSolution(
        eta=eta_values,
        f=f_values,
        fp=fp_values,
        fpp=fpp_values,
    )


# =============================================================================
# Método do Tiro
# =============================================================================

def shooting_error(
    s: float,
    delta_eta: float,
    eta_max: float,
) -> tuple[float, BlasiusSolution]:
    """Calcula E(s) = f'(eta_max; s) - 1 para um chute informado."""
    solution = integrate_blasius(s, delta_eta, eta_max)
    error = solution.fp[-1] - 1.0
    return error, solution


def solve_shooting_method(
    initial_guess: float,
    delta_eta: float,
    eta_max: float,
    tolerance: float,
    max_iterations: int,
) -> ShootingResult:
    """Determina f''(0) pelo Método do Tiro.

    O parâmetro s é atualizado por Newton–Raphson:

        s_(m+1) = s_m - E(s_m) / E'(s_m)

    A derivada de E é aproximada por diferença progressiva:

        E'(s) ≈ [E(s + Delta s) - E(s)] / Delta s

    com:

        Delta s = max(1e-6, 1e-5 |s|).
    """
    s = float(initial_guess)
    history: list[ShootingIteration] = []

    # Variáveis inicializadas para manter o último resultado disponível caso o
    # limite máximo de iterações seja atingido sem convergência.
    error = math.inf
    solution = integrate_blasius(s, delta_eta, eta_max)

    for iteration in range(1, max_iterations + 1):
        error, solution = shooting_error(s, delta_eta, eta_max)

        history.append(
            ShootingIteration(
                iteration=iteration,
                s=s,
                fp_eta_max=solution.fp[-1],
                error=error,
            )
        )

        # Critério de parada do Método do Tiro.
        if abs(error) < tolerance:
            return ShootingResult(
                s=s,
                error=error,
                iterations=iteration,
                converged=True,
                solution=solution,
                history=history,
            )

        # Perturbação utilizada para estimar numericamente E'(s).
        delta_s = max(1.0e-6, 1.0e-5 * abs(s))

        error_perturbed, _ = shooting_error(
            s + delta_s,
            delta_eta,
            eta_max,
        )

        derivative = (error_perturbed - error) / delta_s

        if not math.isfinite(derivative) or abs(derivative) < 1.0e-14:
            raise RuntimeError(
                "A derivada numérica de E(s) ficou nula ou muito pequena. "
                "Tente outro chute inicial ou outro passo de integração."
            )

        next_s = s - error / derivative

        # Para a solução física de Blasius, f''(0) deve ser positivo.
        # Caso Newton produza um valor inválido, usa-se um recuo simples.
        if not math.isfinite(next_s) or next_s <= 0.0:
            next_s = 0.5 * s

        s = next_s

    return ShootingResult(
        s=s,
        error=error,
        iterations=max_iterations,
        converged=False,
        solution=solution,
        history=history,
    )


# =============================================================================
# Pós-processamento
# =============================================================================

def find_eta_99(
    eta_values: list[float],
    velocity_profile: list[float],
) -> Optional[float]:
    """Determina eta_99 por interpolação linear.

    eta_99 é definido pela condição:

        f'(eta_99) = 0.99.
    """
    target = 0.99

    for index in range(1, len(eta_values)):
        fp_left = velocity_profile[index - 1]
        fp_right = velocity_profile[index]

        if fp_left <= target <= fp_right:
            eta_left = eta_values[index - 1]
            eta_right = eta_values[index]

            # Interpolação linear entre os dois pontos vizinhos.
            return eta_left + (
                (target - fp_left)
                * (eta_right - eta_left)
                / (fp_right - fp_left)
            )

    return None


def save_solution_file(solution: BlasiusSolution, output_dir: Path) -> Path:
    """Salva eta, f, f' e f'' em blasius_solution.dat."""
    file_path = output_dir / "blasius_solution.dat"

    with file_path.open("w", encoding="utf-8") as file:
        file.write("# eta f fp fpp\n")

        for eta, f, fp, fpp in zip(
            solution.eta,
            solution.f,
            solution.fp,
            solution.fpp,
        ):
            file.write(
                f"{eta:.10f} "
                f"{f:.12e} "
                f"{fp:.12e} "
                f"{fpp:.12e}\n"
            )

    return file_path


def save_shooting_log(
    history: list[ShootingIteration],
    output_dir: Path,
) -> Path:
    """Salva o histórico iterativo em shooting_log.dat."""
    file_path = output_dir / "shooting_log.dat"

    with file_path.open("w", encoding="utf-8") as file:
        file.write("# iter s fp_eta_max erro\n")

        for item in history:
            file.write(
                f"{item.iteration:d} "
                f"{item.s:.12e} "
                f"{item.fp_eta_max:.12e} "
                f"{item.error:.12e}\n"
            )

    return file_path


def generate_plots(
    solution: BlasiusSolution,
    output_dir: Path,
) -> list[Path]:
    """Gera um arquivo PNG separado para cada perfil de similaridade.

    Matplotlib é uma dependência opcional. A ausência da biblioteca não impede
    a execução da solução numérica nem a geração dos arquivos .dat.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "\nAviso: matplotlib não está instalado; "
            "os gráficos não foram gerados."
        )
        return []

    profiles = [
        (solution.f, r"$f(\eta)$", "perfil_f.png"),
        (solution.fp, r"$f'(\eta)=u/U_\infty$", "perfil_fp.png"),
        (solution.fpp, r"$f''(\eta)$", "perfil_fpp.png"),
    ]

    generated_files: list[Path] = []

    for values, y_label, filename in profiles:
        file_path = output_dir / filename

        plt.figure()
        plt.plot(solution.eta, values)
        plt.xlabel(r"$\eta$")
        plt.ylabel(y_label)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(file_path, dpi=200)
        plt.close()

        generated_files.append(file_path)

    return generated_files


# =============================================================================
# Entrada de dados e interface de linha de comando
# =============================================================================

def read_float(prompt: str, default: float) -> float:
    """Lê um número real no modo interativo, aceitando vírgula decimal."""
    text = input(f"{prompt} [{default}]: ").strip()

    if not text:
        return float(default)

    return float(text.replace(",", "."))


def read_int(prompt: str, default: int) -> int:
    """Lê um número inteiro no modo interativo."""
    text = input(f"{prompt} [{default}]: ").strip()

    if not text:
        return int(default)

    return int(text)


def validate_parameters(
    initial_guess: float,
    delta_eta: float,
    eta_max: float,
    tolerance: float,
    max_iterations: int,
    reynolds_x: float,
) -> None:
    """Valida os parâmetros numéricos antes de iniciar os cálculos."""
    if not math.isfinite(initial_guess) or initial_guess <= 0.0:
        raise ValueError("O chute inicial s deve ser positivo e finito.")

    if not math.isfinite(delta_eta) or delta_eta <= 0.0:
        raise ValueError("Delta eta deve ser positivo e finito.")

    if not math.isfinite(eta_max) or eta_max <= 0.0:
        raise ValueError("eta_max deve ser positivo e finito.")

    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("A tolerância deve ser positiva e finita.")

    if max_iterations <= 0:
        raise ValueError("O número máximo de iterações deve ser positivo.")

    if not math.isfinite(reynolds_x) or reynolds_x <= 0.0:
        raise ValueError("Re_x deve ser positivo e finito.")


def build_argument_parser() -> argparse.ArgumentParser:
    """Cria a interface de linha de comando do programa."""
    parser = argparse.ArgumentParser(
        description=(
            "Resolve a equação de Blasius pelo Método do Tiro "
            "com integração RK4."
        )
    )

    parser.add_argument(
        "--s0",
        type=float,
        default=0.3,
        help="chute inicial s = f''(0) (padrão: 0.3)",
    )
    parser.add_argument(
        "--deta",
        type=float,
        default=0.01,
        help="passo de integração Delta eta (padrão: 0.01)",
    )
    parser.add_argument(
        "--eta-max",
        type=float,
        default=10.0,
        help="limite superior do domínio de eta (padrão: 10)",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1.0e-10,
        help="tolerância do Método do Tiro (padrão: 1e-10)",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=50,
        help="número máximo de iterações do tiro (padrão: 50)",
    )
    parser.add_argument(
        "--re-x",
        type=float,
        default=1.0e5,
        help="número de Reynolds local Re_x (padrão: 1e5)",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("."),
        help="pasta para salvar os resultados (padrão: pasta atual)",
    )
    parser.add_argument(
        "--sem-plots",
        action="store_true",
        help="não gerar os arquivos PNG",
    )
    parser.add_argument(
        "--interativo",
        action="store_true",
        help="solicitar os parâmetros pelo terminal",
    )

    return parser


def main() -> int:
    """Executa o fluxo completo do PPC5."""
    parser = build_argument_parser()
    args = parser.parse_args()

    # O modo interativo sobrescreve os valores recebidos pela linha de comando.
    if args.interativo:
        print("PPC5 — Equação de Blasius por Método do Tiro + RK4")
        print("Pressione Enter para manter o valor indicado entre colchetes.\n")

        args.s0 = read_float("Chute inicial s = f''(0)", args.s0)
        args.deta = read_float("Passo de integração Delta eta", args.deta)
        args.eta_max = read_float("Valor máximo de eta", args.eta_max)
        args.tol = read_float("Tolerância de convergência", args.tol)
        args.max_iter = read_int(
            "Número máximo de iterações do Método do Tiro",
            args.max_iter,
        )
        args.re_x = read_float("Número de Reynolds local Re_x", args.re_x)

    try:
        validate_parameters(
            initial_guess=args.s0,
            delta_eta=args.deta,
            eta_max=args.eta_max,
            tolerance=args.tol,
            max_iterations=args.max_iter,
            reynolds_x=args.re_x,
        )

        result = solve_shooting_method(
            initial_guess=args.s0,
            delta_eta=args.deta,
            eta_max=args.eta_max,
            tolerance=args.tol,
            max_iterations=args.max_iter,
        )
    except (ValueError, RuntimeError, OverflowError) as error:
        print(f"Erro: {error}", file=sys.stderr)
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)

    solution_file = save_solution_file(result.solution, args.outdir)
    log_file = save_shooting_log(result.history, args.outdir)

    plot_files: list[Path] = []
    if not args.sem_plots:
        plot_files = generate_plots(result.solution, args.outdir)

    eta_99 = find_eta_99(result.solution.eta, result.solution.fp)
    skin_friction = 2.0 * result.s / math.sqrt(args.re_x)

    print("\nResultados finais")
    print("-----------------")
    print(f"Convergência             = {'sim' if result.converged else 'não'}")
    print(f"f''(0) convergido        = {result.s:.12f}")
    print(f"Iterações do tiro        = {result.iterations}")
    print(f"Erro final               = {result.error:.12e}")
    print(f"f'(eta_max)              = {result.solution.fp[-1]:.12f}")

    if eta_99 is not None:
        relative_difference = 100.0 * (eta_99 - 4.92) / 4.92

        print(f"eta_99                   = {eta_99:.12f}")
        print(f"C_delta = eta_99         = {eta_99:.12f}")
        print(
            "Diferença para C_delta=4.92 "
            f"= {relative_difference:.6f}%"
        )
    else:
        print(
            "eta_99                   = não encontrado "
            "no domínio integrado"
        )

    print(f"C_f para Re_x informado  = {skin_friction:.12e}")

    print("\nArquivos gerados")
    print("----------------")
    print(solution_file)
    print(log_file)

    for file_path in plot_files:
        print(file_path)

    if not result.converged:
        print(
            "\nAviso: o limite máximo de iterações foi atingido "
            "antes do critério de convergência.",
            file=sys.stderr,
        )
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
