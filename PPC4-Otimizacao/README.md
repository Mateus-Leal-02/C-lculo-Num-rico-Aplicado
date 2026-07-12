# PPC4 — Otimização Bidimensional sem Restrições

## Problema

Maximizar:

```math
f(x,y)=2xy+2x-x^2-2y^2
```

partindo de:

```math
(x_0,y_0)=(-2,3).
```

O ótimo analítico é:

```math
(x^*,y^*)=(2,1),\qquad f^*=2.
```

## Métodos implementados

### Aclive máximo

Usa a direção do gradiente e realiza uma linha de busca a cada iteração.

### Gradientes conjugados de Fletcher–Reeves

Combina o gradiente atual com a direção anterior:

```math
\beta_k=\frac{g_{k+1}^Tg_{k+1}}{g_k^Tg_k}.
```

## Arquivos

| Arquivo | Descrição |
|---|---|
| `Otimizacao.py` | código principal existente |
| `docs/APC4.pdf` | desenvolvimento manual do aclive máximo |
| `docs/PPC4.pdf` | relatório dos dois métodos |
| `exercicios/exercicio-resolvido.md` | duas iterações do aclive máximo |
| `desafio.md` | aplicação a uma função não quadrática |

## Dependências

```bash
pip install numpy matplotlib
```

## Como executar

```bash
python Otimizacao.py
```

Com outro ponto inicial:

```bash
python Otimizacao.py --x0 -1 --y0 4
```

## Arquivos de saída

- `output1.dat`: histórico do aclive máximo;
- `output2.dat`: histórico de Fletcher–Reeves;
- `function.dat`: malha da função;
- `ppc4_caminhos.png`: trajetórias sobre curvas de nível.

## Validação

- norma do gradiente próxima de zero;
- valor final próximo de `2`;
- coordenadas próximas de `(2,1)`;
- Hessiana negativa definida, confirmando máximo estrito.

## Limitações

- a linha de busca atual explora a estrutura quadrática específica;
- gradientes conjugados não lineares podem exigir reinicialização;
- tolerância muito rígida pode ser incompatível com a precisão numérica.

## Material complementar

- [Exercício resolvido](exercicios/exercicio-resolvido.md)
- [Desafio](desafio.md)
- [APC4](docs/APC4.pdf)
- [PPC4](docs/PPC4.pdf)
