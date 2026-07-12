#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PPC6 - Conducao de calor bidimensional em uma aleta retangular.

Este programa resolve numericamente o problema de conducao de calor em regime
permanente em uma aleta retangular de secao transversal constante. A base da
aleta possui temperatura prescrita e as demais superficies trocam calor por
conveccao com o ambiente.

Equacao governante:

    d2T/dx2 + d2T/dy2 = 0

Condicoes de contorno:

    x = 0       -> T = T_b                         (Dirichlet)
    x = L       -> -k dT/dn = h (T - T_inf)        (Robin)
    y = 0       -> -k dT/dn = h (T - T_inf)        (Robin)
    y = H       -> -k dT/dn = h (T - T_inf)        (Robin)

Metodos implementados:

    1. Montagem do sistema linear por diferencas finitas;
    2. Eliminacao de Gauss com pivoteamento parcial;
    3. Metodo de Liebmann, isto e, Gauss-Seidel sem relaxacao;
    4. Metodo de Liebmann com relaxacao, isto e, SOR;
    5. Estudo do fator de relaxacao omega;
    6. Estudo de refinamento de malha;
    7. Comparacao da linha central com a solucao analitica unidimensional.

O codigo foi escrito de forma propositalmente didatica e comentada. A solucao
numerica, a montagem do sistema e os metodos iterativos foram implementados com
estruturas basicas de Python. Nao e utilizado numpy.linalg.solve ou qualquer
rotina pronta de solucao de sistemas lineares.

Autor: Mateus Leal Silva
Disciplina: Calculo Numerico Aplicado
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


# =============================================================================
# Estruturas de dados
# =============================================================================

@dataclass
class ProblemData:
    """Agrupa os parametros fisicos e numericos do problema.

    A vantagem de usar uma estrutura como esta e evitar passar muitas variaveis
    soltas entre as funcoes. Isso reduz erros de ordem de argumentos e torna o
    codigo mais legivel.
    """

    length: float              # comprimento da aleta, L [m]
    thickness: float           # espessura/altura da aleta, H [m]
    conductivity: float        # condutividade termica, k [W/(m.K)]
    convection: float          # coeficiente convectivo, h [W/(m2.K)]
    base_temperature: float    # temperatura prescrita na base, T_b [graus C]
    ambient_temperature: float # temperatura ambiente, T_inf [graus C]
    nx: int                    # numero de nos na direcao x
    ny: int                    # numero de nos na direcao y
    tolerance: float           # tolerancia de convergencia dos metodos iterativos
    omega: float               # fator de relaxacao do metodo SOR
    max_iterations: int        # limite maximo de iteracoes


@dataclass
class IterativeResult:
    """Armazena a resposta de um metodo iterativo."""

    temperature: list[list[float]]  # matriz T[j][i]
    iterations: int                 # numero de iteracoes executadas
    final_error: float              # erro final calculado
    elapsed_time: float             # tempo de execucao [s]
    converged: bool                 # indica se o criterio de convergencia foi atendido


@dataclass
class DirectResult:
    """Armazena a resposta da eliminacao de Gauss."""

    temperature: list[list[float]]  # matriz T[j][i]
    elapsed_time: float             # tempo de execucao [s]


# =============================================================================
# Funcoes auxiliares de malha e indexacao
# =============================================================================

def global_index(i: int, j: int, nx: int) -> int:
    """Converte o indice bidimensional (i, j) para indice global p.

    O programa utiliza a numeracao por linhas horizontais, com j variando de
    baixo para cima e i variando da esquerda para a direita:

        p = j*nx + i

    Como Python usa indices iniciando em zero, o primeiro no possui p = 0.
    Para apresentar em relatorios, basta usar p + 1.
    """
    return j * nx + i


def inverse_index(p: int, nx: int) -> tuple[int, int]:
    """Converte o indice global p de volta para os indices (i, j)."""
    j = p // nx
    i = p % nx
    return i, j


def mesh_spacing(data: ProblemData) -> tuple[float, float]:
    """Calcula os espacamentos uniformes da malha.

    Com nx nos na direcao x, existem nx - 1 intervalos. Assim:

        dx = L/(nx - 1)

    Analogamente:

        dy = H/(ny - 1)
    """
    dx = data.length / (data.nx - 1)
    dy = data.thickness / (data.ny - 1)
    return dx, dy


def create_coordinate_vectors(data: ProblemData) -> tuple[list[float], list[float]]:
    """Gera os vetores de coordenadas x e y da malha."""
    dx, dy = mesh_spacing(data)
    x_values = [i * dx for i in range(data.nx)]
    y_values = [j * dy for j in range(data.ny)]
    return x_values, y_values


def create_initial_temperature_field(data: ProblemData) -> list[list[float]]:
    """Cria um campo inicial para os metodos iterativos.

    A base e iniciada diretamente com T_b. Os demais pontos sao iniciados com
    uma interpolacao linear simples entre T_b e T_inf ao longo de x. Isso nao e
    obrigatorio, mas normalmente reduz um pouco o numero de iteracoes quando
    comparado com iniciar todo o dominio em T_inf.
    """
    temperature: list[list[float]] = []

    for j in range(data.ny):
        row: list[float] = []
        for i in range(data.nx):
            if i == 0:
                row.append(data.base_temperature)
            else:
                fraction = i / (data.nx - 1)
                value = (
                    data.base_temperature
                    + fraction * (data.ambient_temperature - data.base_temperature)
                )
                row.append(value)
        temperature.append(row)

    return temperature


# =============================================================================
# Montagem do sistema linear [A]{T} = {b}
# =============================================================================

def equation_coefficients_for_node(
    i: int,
    j: int,
    data: ProblemData,
) -> tuple[dict[int, float], float]:
    """Monta a equacao discreta de um unico no da malha.

    Retorna um dicionario de coeficientes da linha da matriz A e o termo
    independente correspondente da matriz b.

    A equacao de Laplace discretizada para um no interno e:

        (T_E - 2T_P + T_W)/dx2 + (T_N - 2T_P + T_S)/dy2 = 0

    reorganizando para diagonal positiva:

        (2/dx2 + 2/dy2)T_P
        - (1/dx2)T_E - (1/dx2)T_W
        - (1/dy2)T_N - (1/dy2)T_S = 0

    Nos contornos convectivos, utiliza-se um no ficticio para aplicar a condicao
    de Robin:

        -k dT/dn = h(T - T_inf)

    A substituicao do no ficticio gera os coeficientes modificados das bordas.
    """
    nx = data.nx
    ny = data.ny
    dx, dy = mesh_spacing(data)

    # Coeficientes basicos das derivadas segundas.
    ax = 1.0 / (dx * dx)
    ay = 1.0 / (dy * dy)

    # Coeficientes associados aos termos convectivos que aparecem apos a
    # eliminacao dos nos ficticios.
    bx = data.convection / (data.conductivity * dx)
    by = data.convection / (data.conductivity * dy)

    p = global_index(i, j, nx)

    # -------------------------------------------------------------------------
    # Condicao de Dirichlet na base da aleta: T = T_b.
    # A linha da matriz fica simplesmente: 1*T_p = T_b.
    # -------------------------------------------------------------------------
    if i == 0:
        return {p: 1.0}, data.base_temperature

    coefficients: dict[int, float] = {}
    rhs = 0.0

    def add_coefficient(node_i: int, node_j: int, value: float) -> None:
        """Soma um coeficiente na posicao correta da linha da matriz."""
        q = global_index(node_i, node_j, nx)
        coefficients[q] = coefficients.get(q, 0.0) + value

    # -------------------------------------------------------------------------
    # Contribuicao da direcao x.
    # -------------------------------------------------------------------------
    if i == nx - 1:
        # Extremidade livre x = L: contorno convectivo.
        # O no ficticio externo e eliminado usando a condicao de Robin.
        add_coefficient(i, j, 2.0 * ax + 2.0 * bx)
        add_coefficient(i - 1, j, -2.0 * ax)
        rhs += 2.0 * bx * data.ambient_temperature
    else:
        # No nao esta na extremidade direita; usa vizinhos leste e oeste.
        add_coefficient(i, j, 2.0 * ax)
        add_coefficient(i - 1, j, -ax)
        add_coefficient(i + 1, j, -ax)

    # -------------------------------------------------------------------------
    # Contribuicao da direcao y.
    # -------------------------------------------------------------------------
    if j == 0:
        # Superficie inferior: contorno convectivo.
        add_coefficient(i, j, 2.0 * ay + 2.0 * by)
        add_coefficient(i, j + 1, -2.0 * ay)
        rhs += 2.0 * by * data.ambient_temperature
    elif j == ny - 1:
        # Superficie superior: contorno convectivo.
        add_coefficient(i, j, 2.0 * ay + 2.0 * by)
        add_coefficient(i, j - 1, -2.0 * ay)
        rhs += 2.0 * by * data.ambient_temperature
    else:
        # No nao esta nas superficies superior ou inferior.
        add_coefficient(i, j, 2.0 * ay)
        add_coefficient(i, j - 1, -ay)
        add_cofficient_north = -ay
        add_coefficient(i, j + 1, add_cofficient_north)

    return coefficients, rhs


def build_linear_system(data: ProblemData) -> tuple[list[list[float]], list[float]]:
    """Monta a matriz densa A e o vetor b do sistema linear.

    O armazenamento denso e usado para fins didaticos e para permitir a
    implementacao direta da eliminacao de Gauss. Para malhas muito refinadas,
    a matriz possui muitos zeros e o ideal seria usar armazenamento esparso.
    """
    number_of_nodes = data.nx * data.ny

    # Cria matriz A inicialmente preenchida por zeros.
    matrix = [[0.0 for _ in range(number_of_nodes)] for _ in range(number_of_nodes)]

    # Cria vetor b inicialmente preenchido por zeros.
    rhs = [0.0 for _ in range(number_of_nodes)]

    for j in range(data.ny):
        for i in range(data.nx):
            p = global_index(i, j, data.nx)
            coefficients, independent_term = equation_coefficients_for_node(i, j, data)

            for q, coefficient in coefficients.items():
                matrix[p][q] = coefficient

            rhs[p] = independent_term

    return matrix, rhs


# =============================================================================
# Eliminacao de Gauss com pivoteamento parcial
# =============================================================================

def gaussian_elimination(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Resolve [A]{x} = {b} por eliminacao de Gauss.

    Esta implementacao utiliza pivoteamento parcial. Em cada coluna, procura-se
    a linha com o maior valor absoluto do coeficiente pivo. Isso melhora a
    estabilidade numerica e evita divisao por valores muito pequenos quando uma
    troca de linhas resolve o problema.
    """
    n = len(rhs)

    # Copias sao criadas para nao modificar os objetos originais.
    a = [row[:] for row in matrix]
    b = rhs[:]

    # -------------------------------------------------------------------------
    # Etapa de eliminacao progressiva: transforma A em matriz triangular superior.
    # -------------------------------------------------------------------------
    for pivot_column in range(n - 1):
        # Procura o melhor pivo na coluna atual.
        pivot_row = pivot_column
        pivot_abs_value = abs(a[pivot_column][pivot_column])

        for candidate_row in range(pivot_column + 1, n):
            candidate_abs_value = abs(a[candidate_row][pivot_column])
            if candidate_abs_value > pivot_abs_value:
                pivot_abs_value = candidate_abs_value
                pivot_row = candidate_row

        if pivot_abs_value < 1.0e-20:
            raise ZeroDivisionError(
                "Pivo numericamente nulo encontrado na eliminacao de Gauss."
            )

        # Troca as linhas, se necessario.
        if pivot_row != pivot_column:
            a[pivot_column], a[pivot_row] = a[pivot_row], a[pivot_column]
            b[pivot_column], b[pivot_row] = b[pivot_row], b[pivot_column]

        pivot = a[pivot_column][pivot_column]

        # Zera os coeficientes abaixo do pivo.
        for row in range(pivot_column + 1, n):
            factor = a[row][pivot_column] / pivot

            # O coeficiente abaixo do pivo e definido explicitamente como zero.
            a[row][pivot_column] = 0.0

            # Atualiza os demais coeficientes da linha.
            for column in range(pivot_column + 1, n):
                a[row][column] -= factor * a[pivot_column][column]

            # Atualiza tambem o lado direito.
            b[row] -= factor * b[pivot_column]

    # -------------------------------------------------------------------------
    # Etapa de substituicao regressiva: resolve de baixo para cima.
    # -------------------------------------------------------------------------
    x = [0.0 for _ in range(n)]

    for row in range(n - 1, -1, -1):
        summation = 0.0

        for column in range(row + 1, n):
            summation += a[row][column] * x[column]

        diagonal = a[row][row]
        if abs(diagonal) < 1.0e-20:
            raise ZeroDivisionError(
                "Diagonal numericamente nula encontrada na substituicao regressiva."
            )

        x[row] = (b[row] - summation) / diagonal

    return x


def vector_to_temperature_field(vector: list[float], nx: int, ny: int) -> list[list[float]]:
    """Converte a solucao vetorial para uma matriz T[j][i]."""
    temperature = [[0.0 for _ in range(nx)] for _ in range(ny)]

    for p, value in enumerate(vector):
        i, j = inverse_index(p, nx)
        temperature[j][i] = value

    return temperature


def solve_by_gauss(data: ProblemData) -> DirectResult:
    """Monta e resolve o sistema linear pelo metodo direto de Gauss."""
    start = time.perf_counter()
    matrix, rhs = build_linear_system(data)
    solution_vector = gaussian_elimination(matrix, rhs)
    elapsed = time.perf_counter() - start

    temperature = vector_to_temperature_field(solution_vector, data.nx, data.ny)
    return DirectResult(temperature=temperature, elapsed_time=elapsed)


# =============================================================================
# Metodos iterativos de Liebmann e SOR
# =============================================================================

def local_update_value(
    i: int,
    j: int,
    temperature: list[list[float]],
    data: ProblemData,
) -> float:
    """Calcula o valor de T no no (i, j) usando a equacao local.

    Esta funcao usa a mesma discretizacao da matriz global, mas sem montar a
    matriz. Ela e conveniente para os metodos de Gauss-Seidel e SOR.

    A equacao geral de uma linha pode ser escrita como:

        a_P T_P - soma(a_vizinho T_vizinho) = b

    Logo:

        T_P = [soma(a_vizinho T_vizinho) + b]/a_P
    """
    # A base prescrita nao deve ser atualizada iterativamente.
    if i == 0:
        return data.base_temperature

    coefficients, rhs = equation_coefficients_for_node(i, j, data)
    p = global_index(i, j, data.nx)
    diagonal = coefficients[p]

    neighbor_sum = 0.0

    for q, coefficient in coefficients.items():
        if q == p:
            continue

        neighbor_i, neighbor_j = inverse_index(q, data.nx)

        # Como os coeficientes dos vizinhos aparecem negativos na matriz, eles
        # sao levados para o outro lado da equacao com sinal positivo.
        neighbor_sum += (-coefficient) * temperature[neighbor_j][neighbor_i]

    return (neighbor_sum + rhs) / diagonal


def solve_by_liebmann(data: ProblemData, omega: float) -> IterativeResult:
    """Resolve o problema por Liebmann/Gauss-Seidel com relaxacao opcional.

    Quando omega = 1, tem-se o metodo de Liebmann sem relaxacao.
    Quando 1 < omega < 2, tem-se sobre-relaxacao sucessiva (SOR).
    Quando 0 < omega < 1, tem-se sub-relaxacao.

    A atualizacao adotada e:

        T_novo = (1 - omega) T_velho + omega T_calculado

    em que T_calculado e o valor obtido diretamente pela equacao de diferencas.
    """
    temperature = create_initial_temperature_field(data)

    start = time.perf_counter()
    final_error = math.inf
    converged = False

    for iteration in range(1, data.max_iterations + 1):
        maximum_change = 0.0
        maximum_reference = 1.0e-12

        # Varre todos os nos da malha. A varredura linha a linha caracteriza o
        # uso imediato dos valores ja atualizados, como no Gauss-Seidel.
        for j in range(data.ny):
            for i in range(data.nx):
                old_value = temperature[j][i]

                calculated_value = local_update_value(i, j, temperature, data)
                new_value = (1.0 - omega) * old_value + omega * calculated_value

                # Garante novamente a condicao de Dirichlet da base.
                if i == 0:
                    new_value = data.base_temperature

                temperature[j][i] = new_value

                change = abs(new_value - old_value)
                maximum_change = max(maximum_change, change)
                maximum_reference = max(maximum_reference, abs(new_value))

        final_error = maximum_change / maximum_reference

        if final_error < data.tolerance:
            converged = True
            break

    elapsed = time.perf_counter() - start

    return IterativeResult(
        temperature=temperature,
        iterations=iteration,
        final_error=final_error,
        elapsed_time=elapsed,
        converged=converged,
    )


# =============================================================================
# Comparacoes, pos-processamento e solucao analitica 1D
# =============================================================================

def max_abs_difference(
    first: list[list[float]],
    second: list[list[float]],
) -> float:
    """Calcula a maior diferenca absoluta entre dois campos de temperatura."""
    maximum = 0.0

    for row_a, row_b in zip(first, second):
        for value_a, value_b in zip(row_a, row_b):
            maximum = max(maximum, abs(value_a - value_b))

    return maximum


def interpolate_temperature_at_y(
    temperature: list[list[float]],
    y_target: float,
    data: ProblemData,
) -> list[float]:
    """Extrai T(x, y_target) por interpolacao linear na direcao y."""
    _, y_values = create_coordinate_vectors(data)

    if y_target <= y_values[0]:
        return temperature[0][:]

    if y_target >= y_values[-1]:
        return temperature[-1][:]

    for j in range(1, data.ny):
        if y_values[j - 1] <= y_target <= y_values[j]:
            y0 = y_values[j - 1]
            y1 = y_values[j]
            fraction = (y_target - y0) / (y1 - y0)

            line = []
            for i in range(data.nx):
                value = temperature[j - 1][i] * (1.0 - fraction) + temperature[j][i] * fraction
                line.append(value)

            return line

    # Por seguranca, retorna a ultima linha se algo numericamente inesperado ocorrer.
    return temperature[-1][:]


def analytical_1d_temperature(x: float, data: ProblemData) -> float:
    """Calcula a solucao analitica 1D classica da aleta com ponta convectiva.

    A forma adimensional e:

        theta(x) = [cosh(m(L-x)) + h/(m k) sinh(m(L-x))]
                   /[cosh(mL) + h/(m k) sinh(mL)]

    em que:

        theta = (T - T_inf)/(T_b - T_inf)
        m = sqrt(h P/(k A_c))

    Para compatibilizar com o modelo bidimensional usado aqui, considera-se
    largura unitaria na direcao fora do plano e conveccao pelas faces superior e
    inferior. Assim:

        P = 2 * largura_unitaria
        A_c = H * largura_unitaria

    Logo P/A_c = 2/H.
    """
    unit_width = 1.0
    perimeter = 2.0 * unit_width
    area = data.thickness * unit_width

    m = math.sqrt(data.convection * perimeter / (data.conductivity * area))

    if m == 0.0:
        return data.base_temperature

    numerator = (
        math.cosh(m * (data.length - x))
        + data.convection / (m * data.conductivity)
        * math.sinh(m * (data.length - x))
    )

    denominator = (
        math.cosh(m * data.length)
        + data.convection / (m * data.conductivity)
        * math.sinh(m * data.length)
    )

    theta = numerator / denominator
    temperature = data.ambient_temperature + theta * (
        data.base_temperature - data.ambient_temperature
    )

    return temperature


def central_line_comparison(
    temperature: list[list[float]],
    data: ProblemData,
) -> tuple[list[tuple[float, float, float, float]], float]:
    """Compara T numerico na linha central com a solucao analitica 1D.

    Retorna uma lista com:

        x, T_numerico, T_analitico, erro_percentual

    e tambem o erro percentual medio.
    """
    x_values, _ = create_coordinate_vectors(data)
    y_mid = 0.5 * data.thickness
    numerical_line = interpolate_temperature_at_y(temperature, y_mid, data)

    rows: list[tuple[float, float, float, float]] = []
    errors = []

    for x, numerical_temperature in zip(x_values, numerical_line):
        analytical_temperature = analytical_1d_temperature(x, data)

        theta_numerical = (
            (numerical_temperature - data.ambient_temperature)
            / (data.base_temperature - data.ambient_temperature)
        )
        theta_analytical = (
            (analytical_temperature - data.ambient_temperature)
            / (data.base_temperature - data.ambient_temperature)
        )

        denominator = max(abs(theta_analytical), 1.0e-12)
        percent_error = abs(theta_numerical - theta_analytical) / denominator * 100.0

        rows.append((x, numerical_temperature, analytical_temperature, percent_error))
        errors.append(percent_error)

    mean_error = sum(errors) / len(errors)
    return rows, mean_error


# =============================================================================
# Estudos numericos exigidos no PPC6
# =============================================================================

def relaxation_study(data: ProblemData, omegas: list[float]) -> list[tuple[float, int, float, float, bool]]:
    """Executa o estudo do efeito do fator de relaxacao omega."""
    results = []

    for omega in omegas:
        result = solve_by_liebmann(data, omega=omega)
        results.append(
            (
                omega,
                result.iterations,
                result.final_error,
                result.elapsed_time,
                result.converged,
            )
        )

    return results


def mesh_refinement_study(
    base_data: ProblemData,
    meshes: list[tuple[int, int]],
) -> list[tuple[int, int, int, float, float, float, bool]]:
    """Executa o estudo de refinamento de malha.

    Para cada malha, resolve-se o problema com SOR e calcula-se o erro percentual
    medio da linha central em relacao ao modelo analitico 1D.
    """
    results = []

    for nx, ny in meshes:
        data = ProblemData(
            length=base_data.length,
            thickness=base_data.thickness,
            conductivity=base_data.conductivity,
            convection=base_data.convection,
            base_temperature=base_data.base_temperature,
            ambient_temperature=base_data.ambient_temperature,
            nx=nx,
            ny=ny,
            tolerance=base_data.tolerance,
            omega=base_data.omega,
            max_iterations=base_data.max_iterations,
        )

        result = solve_by_liebmann(data, omega=data.omega)
        _, mean_error = central_line_comparison(result.temperature, data)

        results.append(
            (
                nx,
                ny,
                result.iterations,
                result.final_error,
                result.elapsed_time,
                mean_error,
                result.converged,
            )
        )

    return results


# =============================================================================
# Escrita de arquivos de saida
# =============================================================================

def save_temperature_field(
    filename: Path,
    temperature: list[list[float]],
    data: ProblemData,
) -> None:
    """Salva x, y e T para todos os nos da malha."""
    x_values, y_values = create_coordinate_vectors(data)

    with filename.open("w", encoding="utf-8") as file:
        file.write("# x y T\n")

        for j, y in enumerate(y_values):
            for i, x in enumerate(x_values):
                file.write(f"{x:.12e} {y:.12e} {temperature[j][i]:.12e}\n")

            # Linha em branco separa as linhas horizontais da malha.
            file.write("\n")


def save_central_line(
    filename: Path,
    rows: list[tuple[float, float, float, float]],
) -> None:
    """Salva a comparacao da linha central em arquivo de texto."""
    with filename.open("w", encoding="utf-8") as file:
        file.write("# x T_numerico T_analitico erro_percentual\n")

        for x, numerical, analytical, error in rows:
            file.write(
                f"{x:.12e} {numerical:.12e} {analytical:.12e} {error:.12e}\n"
            )


def save_relaxation_study(
    filename: Path,
    rows: list[tuple[float, int, float, float, bool]],
) -> None:
    """Salva os resultados do estudo de omega."""
    with filename.open("w", encoding="utf-8") as file:
        file.write("# omega iteracoes erro_final tempo_s convergiu\n")

        for omega, iterations, error, elapsed, converged in rows:
            file.write(
                f"{omega:.6f} {iterations:d} {error:.12e} "
                f"{elapsed:.12e} {int(converged):d}\n"
            )


def save_mesh_refinement_study(
    filename: Path,
    rows: list[tuple[int, int, int, float, float, float, bool]],
) -> None:
    """Salva os resultados do estudo de refinamento de malha."""
    with filename.open("w", encoding="utf-8") as file:
        file.write("# nx ny iteracoes erro_final tempo_s erro_medio_percentual convergiu\n")

        for nx, ny, iterations, error, elapsed, mean_error, converged in rows:
            file.write(
                f"{nx:d} {ny:d} {iterations:d} {error:.12e} "
                f"{elapsed:.12e} {mean_error:.12e} {int(converged):d}\n"
            )


def save_summary(
    filename: Path,
    data: ProblemData,
    gauss: DirectResult,
    liebmann: IterativeResult,
    sor: IterativeResult,
    central_line_error: float,
) -> None:
    """Salva um resumo dos resultados principais."""
    maximum_difference_liebmann = max_abs_difference(gauss.temperature, liebmann.temperature)
    maximum_difference_sor = max_abs_difference(gauss.temperature, sor.temperature)

    with filename.open("w", encoding="utf-8") as file:
        file.write("Resumo do PPC6 - Aleta retangular\n")
        file.write("=================================\n\n")
        file.write(f"L = {data.length:.6g} m\n")
        file.write(f"H = {data.thickness:.6g} m\n")
        file.write(f"k = {data.conductivity:.6g} W/(m.K)\n")
        file.write(f"h = {data.convection:.6g} W/(m2.K)\n")
        file.write(f"T_b = {data.base_temperature:.6g} graus C\n")
        file.write(f"T_inf = {data.ambient_temperature:.6g} graus C\n")
        file.write(f"nx = {data.nx}, ny = {data.ny}\n")
        file.write(f"tolerancia = {data.tolerance:.3e}\n")
        file.write(f"omega = {data.omega:.6g}\n\n")

        file.write("Resultados dos metodos\n")
        file.write("----------------------\n")
        file.write(f"Gauss: tempo = {gauss.elapsed_time:.6e} s\n")
        file.write(
            f"Liebmann: iteracoes = {liebmann.iterations}, "
            f"erro = {liebmann.final_error:.6e}, "
            f"tempo = {liebmann.elapsed_time:.6e} s, "
            f"convergiu = {liebmann.converged}\n"
        )
        file.write(
            f"SOR: iteracoes = {sor.iterations}, "
            f"erro = {sor.final_error:.6e}, "
            f"tempo = {sor.elapsed_time:.6e} s, "
            f"convergiu = {sor.converged}\n"
        )
        file.write(f"Diferenca maxima Gauss x Liebmann = {maximum_difference_liebmann:.6e} graus C\n")
        file.write(f"Diferenca maxima Gauss x SOR = {maximum_difference_sor:.6e} graus C\n")
        file.write(f"Erro percentual medio da linha central = {central_line_error:.6e}%\n")


# =============================================================================
# Graficos
# =============================================================================

def generate_plots(
    temperature: list[list[float]],
    central_line_rows: list[tuple[float, float, float, float]],
    data: ProblemData,
    output_dir: Path,
) -> list[Path]:
    """Gera os graficos exigidos no enunciado.

    A biblioteca matplotlib e usada apenas para visualizacao. Se ela nao estiver
    instalada, os arquivos de dados continuam sendo gerados normalmente.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Aviso: matplotlib nao encontrado. Os graficos nao foram gerados.")
        return []

    x_values, y_values = create_coordinate_vectors(data)

    # Matrizes X e Y no mesmo formato de temperature[j][i].
    x_grid = [[x for x in x_values] for _ in y_values]
    y_grid = [[y for _ in x_values] for y in y_values]

    generated_files: list[Path] = []

    # -------------------------------------------------------------------------
    # Mapa de temperatura bidimensional.
    # -------------------------------------------------------------------------
    map_path = output_dir / "mapa_temperatura.png"
    plt.figure()
    image = plt.imshow(
        temperature,
        extent=[0.0, data.length, 0.0, data.thickness],
        origin="lower",
        aspect="auto",
    )
    plt.colorbar(image, label="Temperatura [graus C]")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Mapa de temperatura bidimensional")
    plt.tight_layout()
    plt.savefig(map_path, dpi=200)
    plt.close()
    generated_files.append(map_path)

    # -------------------------------------------------------------------------
    # Curvas de nivel ou contornos isotermicos.
    # -------------------------------------------------------------------------
    contour_path = output_dir / "contornos_isotermicos.png"
    plt.figure()
    contour = plt.contourf(x_grid, y_grid, temperature, levels=20)
    plt.colorbar(contour, label="Temperatura [graus C]")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Contornos isotermicos")
    plt.tight_layout()
    plt.savefig(contour_path, dpi=200)
    plt.close()
    generated_files.append(contour_path)

    # -------------------------------------------------------------------------
    # Distribuicao de temperatura na linha central.
    # -------------------------------------------------------------------------
    line_path = output_dir / "linha_central.png"
    x_line = [row[0] for row in central_line_rows]
    numerical_line = [row[1] for row in central_line_rows]
    analytical_line = [row[2] for row in central_line_rows]

    plt.figure()
    plt.plot(x_line, numerical_line, marker="o", label="Diferencas finitas 2D")
    plt.plot(x_line, analytical_line, linestyle="--", label="Analitica 1D")
    plt.xlabel("x [m]")
    plt.ylabel("Temperatura [graus C]")
    plt.title("Linha central da aleta")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(line_path, dpi=200)
    plt.close()
    generated_files.append(line_path)

    return generated_files


# =============================================================================
# Interface de entrada
# =============================================================================

def read_float(prompt: str, default: float) -> float:
    """Le um numero real do terminal, aceitando virgula ou ponto decimal."""
    text = input(f"{prompt} [{default}]: ").strip()

    if not text:
        return default

    return float(text.replace(",", "."))


def read_int(prompt: str, default: int) -> int:
    """Le um numero inteiro do terminal."""
    text = input(f"{prompt} [{default}]: ").strip()

    if not text:
        return default

    return int(text)


def interactive_input(default_data: ProblemData) -> ProblemData:
    """Solicita os dados do problema ao usuario."""
    print("PPC6 - Aleta retangular por diferencas finitas")
    print("Pressione Enter para manter o valor padrao indicado.\n")

    length = read_float("Comprimento da aleta L [m]", default_data.length)
    thickness = read_float("Espessura/altura da aleta H [m]", default_data.thickness)
    conductivity = read_float("Condutividade termica k [W/(m.K)]", default_data.conductivity)
    convection = read_float("Coeficiente convectivo h [W/(m2.K)]", default_data.convection)
    base_temperature = read_float("Temperatura da base T_b [graus C]", default_data.base_temperature)
    ambient_temperature = read_float("Temperatura ambiente T_inf [graus C]", default_data.ambient_temperature)
    nx = read_int("Numero de nos na direcao x", default_data.nx)
    ny = read_int("Numero de nos na direcao y", default_data.ny)
    tolerance = read_float("Tolerancia de convergencia", default_data.tolerance)
    omega = read_float("Fator de relaxacao omega", default_data.omega)
    max_iterations = read_int("Numero maximo de iteracoes", default_data.max_iterations)

    return ProblemData(
        length=length,
        thickness=thickness,
        conductivity=conductivity,
        convection=convection,
        base_temperature=base_temperature,
        ambient_temperature=ambient_temperature,
        nx=nx,
        ny=ny,
        tolerance=tolerance,
        omega=omega,
        max_iterations=max_iterations,
    )


def validate_data(data: ProblemData) -> None:
    """Verifica se os parametros fornecidos sao fisicamente e numericamente validos."""
    if data.length <= 0.0:
        raise ValueError("O comprimento L deve ser positivo.")
    if data.thickness <= 0.0:
        raise ValueError("A espessura H deve ser positiva.")
    if data.conductivity <= 0.0:
        raise ValueError("A condutividade termica k deve ser positiva.")
    if data.convection < 0.0:
        raise ValueError("O coeficiente convectivo h nao pode ser negativo.")
    if data.nx < 3:
        raise ValueError("Use pelo menos 3 nos na direcao x.")
    if data.ny < 3:
        raise ValueError("Use pelo menos 3 nos na direcao y.")
    if data.tolerance <= 0.0:
        raise ValueError("A tolerancia deve ser positiva.")
    if not (0.0 < data.omega < 2.0):
        raise ValueError("Para SOR, recomenda-se 0 < omega < 2.")
    if data.max_iterations <= 0:
        raise ValueError("O numero maximo de iteracoes deve ser positivo.")


def build_argument_parser() -> argparse.ArgumentParser:
    """Cria a interface de linha de comando."""
    parser = argparse.ArgumentParser(
        description="Resolve a aleta retangular 2D por diferencas finitas."
    )

    parser.add_argument("--L", type=float, default=0.12, help="comprimento da aleta [m]")
    parser.add_argument("--H", type=float, default=0.02, help="espessura/altura da aleta [m]")
    parser.add_argument("--k", type=float, default=205.0, help="condutividade termica [W/(m.K)]")
    parser.add_argument("--h", type=float, default=25.0, help="coeficiente convectivo [W/(m2.K)]")
    parser.add_argument("--Tb", type=float, default=100.0, help="temperatura da base [graus C]")
    parser.add_argument("--Tinf", type=float, default=25.0, help="temperatura ambiente [graus C]")
    parser.add_argument("--nx", type=int, default=21, help="numero de nos na direcao x")
    parser.add_argument("--ny", type=int, default=9, help="numero de nos na direcao y")
    parser.add_argument("--tol", type=float, default=1.0e-8, help="tolerancia iterativa")
    parser.add_argument("--omega", type=float, default=1.60, help="fator de relaxacao SOR")
    parser.add_argument("--max-iter", type=int, default=100000, help="maximo de iteracoes")
    parser.add_argument("--outdir", type=Path, default=Path("resultados_ppc6"), help="pasta de saida")
    parser.add_argument("--interativo", action="store_true", help="solicita os parametros pelo terminal")
    parser.add_argument("--sem-plots", action="store_true", help="nao gera figuras PNG")
    parser.add_argument("--pular-estudos", action="store_true", help="nao executa estudos de omega e malha")

    return parser


# =============================================================================
# Funcao principal
# =============================================================================

def main() -> int:
    """Executa o fluxo completo do PPC6."""
    parser = build_argument_parser()
    args = parser.parse_args()

    data = ProblemData(
        length=args.L,
        thickness=args.H,
        conductivity=args.k,
        convection=args.h,
        base_temperature=args.Tb,
        ambient_temperature=args.Tinf,
        nx=args.nx,
        ny=args.ny,
        tolerance=args.tol,
        omega=args.omega,
        max_iterations=args.max_iter,
    )

    if args.interativo:
        data = interactive_input(data)

    validate_data(data)
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("\nMontando e resolvendo o problema da aleta retangular...")
    print(f"Malha: {data.nx} x {data.ny} = {data.nx * data.ny} nos")

    # -------------------------------------------------------------------------
    # Metodo direto: eliminacao de Gauss.
    # -------------------------------------------------------------------------
    gauss = solve_by_gauss(data)
    print(f"Gauss concluido em {gauss.elapsed_time:.6e} s")

    # -------------------------------------------------------------------------
    # Liebmann sem relaxacao: omega = 1.
    # -------------------------------------------------------------------------
    liebmann = solve_by_liebmann(data, omega=1.0)
    print(
        "Liebmann: "
        f"iteracoes = {liebmann.iterations}, "
        f"erro = {liebmann.final_error:.6e}, "
        f"tempo = {liebmann.elapsed_time:.6e} s"
    )

    # -------------------------------------------------------------------------
    # Liebmann com relaxacao: omega definido pelo usuario.
    # -------------------------------------------------------------------------
    sor = solve_by_liebmann(data, omega=data.omega)
    print(
        "SOR: "
        f"iteracoes = {sor.iterations}, "
        f"erro = {sor.final_error:.6e}, "
        f"tempo = {sor.elapsed_time:.6e} s"
    )

    # -------------------------------------------------------------------------
    # Comparacao com a solucao analitica 1D na linha central.
    # -------------------------------------------------------------------------
    central_rows, mean_error = central_line_comparison(sor.temperature, data)
    print(f"Erro percentual medio na linha central = {mean_error:.6f}%")

    # -------------------------------------------------------------------------
    # Escrita dos arquivos principais.
    # -------------------------------------------------------------------------
    save_temperature_field(args.outdir / "temperature_gauss.dat", gauss.temperature, data)
    save_temperature_field(args.outdir / "temperature_liebmann.dat", liebmann.temperature, data)
    save_temperature_field(args.outdir / "temperature_sor.dat", sor.temperature, data)
    save_central_line(args.outdir / "linha_central.dat", central_rows)
    save_summary(args.outdir / "resumo_resultados.txt", data, gauss, liebmann, sor, mean_error)

    # -------------------------------------------------------------------------
    # Estudos adicionais solicitados no enunciado.
    # -------------------------------------------------------------------------
    if not args.pular_estudos:
        omegas = [0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 1.9]
        omega_rows = relaxation_study(data, omegas)
        save_relaxation_study(args.outdir / "estudo_relaxacao.dat", omega_rows)

        # Malhas impares ajudam a existir uma linha central exatamente em y = H/2.
        mesh_rows = mesh_refinement_study(
            data,
            meshes=[(11, 5), (21, 9), (31, 13)],
        )
        save_mesh_refinement_study(args.outdir / "estudo_malha.dat", mesh_rows)

    # -------------------------------------------------------------------------
    # Figuras.
    # -------------------------------------------------------------------------
    if not args.sem_plots:
        generated = generate_plots(sor.temperature, central_rows, data, args.outdir)
        if generated:
            print("Graficos gerados:")
            for path in generated:
                print(f"  {path}")

    print("\nArquivos salvos em:", args.outdir)
    print("Diferenca maxima Gauss x Liebmann:", f"{max_abs_difference(gauss.temperature, liebmann.temperature):.6e}")
    print("Diferenca maxima Gauss x SOR:", f"{max_abs_difference(gauss.temperature, sor.temperature):.6e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
