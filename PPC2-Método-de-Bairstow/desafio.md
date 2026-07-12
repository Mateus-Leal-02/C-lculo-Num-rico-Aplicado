# Desafio — robustez aos chutes iniciais

Para o polinômio da APC2, avalie uma grade de chutes:

```text
r0 em [-3, 3]
s0 em [-3, 3]
```

## Registre para cada ponto

- convergiu ou não;
- número de iterações;
- primeiro par de raízes encontrado;
- valor final dos restos `b0` e `b1`.

## Entregas

1. mapa de convergência;
2. mapa do número de iterações;
3. identificação de regiões problemáticas;
4. comparação entre tolerâncias `1e-6`, `1e-8` e `1e-10`;
5. proposta de estratégia automática para reiniciar chutes que falham.

<details>
<summary>Dica</summary>

Além de colorir pelo número de iterações, use uma cor separada para pontos em que o determinante do sistema de correção fica próximo de zero.

</details>
