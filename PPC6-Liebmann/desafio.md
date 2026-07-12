# Desafio — fator ótimo de relaxação

Determine experimentalmente o valor de `omega` que minimiza o número de iterações do SOR.

## Malhas mínimas

```text
11 x 5
21 x 9
41 x 17
81 x 33
```

## Faixa de busca

```text
omega = 1.00, 1.05, 1.10, ..., 1.95
```

## Registre

- convergiu ou divergiu;
- número de iterações;
- tempo;
- erro final;
- resíduo discreto máximo;
- diferença para a solução de Gauss em uma malha pequena.

## Entregas

1. gráfico de iterações versus `omega` para cada malha;
2. tabela do melhor `omega` observado;
3. análise da tendência com o refinamento;
4. discussão sobre estabilidade quando `omega` se aproxima de `2`;
5. comparação entre critério por incremento e por resíduo.

<details>
<summary>Dica</summary>

O valor ótimo depende da geometria, da malha e das condições de contorno. Não assuma que `omega = 1,6` continuará ótimo para todas as discretizações.

</details>
