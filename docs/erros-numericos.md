# Erros numéricos

## 1. Erro de truncamento

Surge ao substituir uma operação matemática exata por uma aproximação finita. Exemplos:

- truncar uma série de Taylor;
- aproximar derivadas por diferenças finitas;
- integrar uma EDO com passo finito;
- representar o infinito de Blasius por um `eta_max` finito.

A redução do passo ou o refinamento da malha geralmente diminui esse erro, mas aumenta o custo computacional.

## 2. Erro de arredondamento

Computadores armazenam números reais com precisão finita. Subtrações entre valores próximos, matrizes mal condicionadas e iterações muito longas podem amplificar o arredondamento.

## 3. Erro iterativo

Métodos como Bairstow, Newton–Raphson, Liebmann e SOR são interrompidos antes da solução matemática exata. A tolerância deve ser compatível com:

- precisão dos dados;
- erro de discretização;
- custo computacional;
- grandeza física analisada.

## 4. Erro de modelagem

É a diferença entre o fenômeno real e as hipóteses do modelo. Exemplos:

- propriedades constantes;
- geometria unidimensional;
- escoamento laminar;
- ausência de geração interna;
- domínio infinito aproximado por domínio finito.

Refinar a malha não corrige um modelo físico inadequado.

## 5. Erro absoluto e relativo

```math
E_a = |x_{num}-x_{ref}|
```

```math
E_r = \frac{|x_{num}-x_{ref}|}{|x_{ref}|}
```

O erro relativo não deve ser usado diretamente quando o valor de referência é zero ou muito próximo de zero.

## 6. Norma de erro

Para um vetor de diferenças `e`:

```math
\|e\|_\infty = \max_i |e_i|
```

```math
\|e\|_2 = \sqrt{\sum_i e_i^2}
```

A norma infinita destaca o pior ponto. A norma 2 resume o comportamento global.

## 7. Estudo de ordem

Para um método com erro aproximado `E(h) ≈ C h^p`, a ordem observada pode ser estimada por:

```math
p \approx \frac{\ln(E_h/E_{h/2})}{\ln 2}
```

Esse teste é especialmente relevante no PPC1 e em estudos de refinamento dos PPCs 3 e 6.
