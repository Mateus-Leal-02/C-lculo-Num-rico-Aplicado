# Desafio — ordem observada do RK4

Implemente um estudo com:

```text
St = 1
Re_s = 0
t_final = 5
h = 0.8, 0.4, 0.2, 0.1, 0.05
```

## Entregas

1. erro máximo em cada malha temporal;
2. tabela com `h`, erro e ordem observada;
3. gráfico `log(erro)` versus `log(h)`;
4. regressão linear para estimar a inclinação;
5. discussão sobre o ponto em que o arredondamento impede nova redução do erro.

## Resultado esperado

No intervalo em que o erro de truncamento domina, a inclinação deve se aproximar de `4`.

<details>
<summary>Dica</summary>

Use a solução exata `v*(t*) = 1 - exp(-t*/St)` e calcule a ordem entre refinamentos sucessivos por `p = ln(Eh/Eh2)/ln(2)`.

</details>
