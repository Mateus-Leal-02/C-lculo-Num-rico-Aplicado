# PPC2 — Raízes de Polinômios pelo Método de Bairstow

## Sobre o projeto

O código determina raízes reais e pares complexos conjugados por meio da divisão iterativa do polinômio por fatores quadráticos reais:

```math
D(x)=x^2-rx-s.
```

O método é aplicado ao polinômio característico do sistema massa–mola–amortecedor de dois graus de liberdade:

```math
P(\lambda)=2\lambda^4+5\lambda^3+12\lambda^2+8\lambda+8.
```

## Etapas do algoritmo

1. escolher `r0` e `s0`;
2. calcular os coeficientes da divisão sintética `b`;
3. calcular os coeficientes auxiliares `c`;
4. resolver o sistema para `Δr` e `Δs`;
5. atualizar `r` e `s`;
6. extrair duas raízes do fator convergido;
7. deflacionar o polinômio;
8. repetir até grau 2 ou 1.

## Arquivos

| Arquivo | Descrição |
|---|---|
| `Bairstow.py` | implementação principal já existente |
| `docs/APC2.pdf` | formulação do polinômio característico |
| `docs/PPC2.pdf` | relatório da implementação e validação |
| `exercicios/exercicio-resolvido.md` | uma atualização completa de `(r,s)` |
| `desafio.md` | estudo das bacias de convergência |

## Dependências

```bash
pip install numpy matplotlib
```

## Como executar

```bash
python Bairstow.py
```

A geração do fractal pode exigir mais tempo que as demais análises.

## Validação

- polinômio de grau 7 construído a partir de raízes conhecidas;
- reconstrução do polinômio a partir das raízes calculadas;
- verificação do resíduo `|P(raiz)|`;
- interpretação da parte real dos autovalores do sistema dinâmico.

## Critério de convergência

A atualização pode ser encerrada quando os incrementos relativos de `r` e `s` ficam abaixo da tolerância. Também é recomendável conferir os restos `b0` e `b1`.

## Limitações

- sensibilidade aos chutes iniciais;
- possível degeneração do sistema de correção;
- acúmulo de erro durante a deflação;
- dificuldades com raízes múltiplas ou muito próximas.

## Material complementar

- [Exercício resolvido](exercicios/exercicio-resolvido.md)
- [Desafio](desafio.md)
- [APC2](docs/APC2.pdf)
- [PPC2](docs/PPC2.pdf)
