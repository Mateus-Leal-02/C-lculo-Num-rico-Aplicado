# Desafio — Newton–Raphson versus secante

Substitua a atualização do parâmetro de tiro por um método da secante:

```math
s_{k+1}=s_k-E(s_k)
\frac{s_k-s_{k-1}}{E(s_k)-E(s_{k-1})}.
```

## Compare

- número de integrações completas da equação de Blasius;
- número de iterações externas;
- erro final;
- sensibilidade aos chutes;
- falhas por denominador pequeno;
- tempo de execução.

## Casos mínimos

```text
(0,20; 0,40)
(0,30; 0,35)
(0,10; 0,60)
```

## Discussão

Newton com derivada por diferença finita também exige integrações adicionais para estimar `E'(s)`. Portanto, compare o número total de avaliações de `E`, e não apenas o número de iterações externas.
