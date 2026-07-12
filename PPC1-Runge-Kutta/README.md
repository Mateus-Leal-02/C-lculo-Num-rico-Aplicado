# PPC1 — Sedimentação de Esfera por Runge–Kutta de 4ª Ordem

## Sobre o projeto

O programa resolve o movimento vertical adimensional de uma partícula esférica em um fluido para baixo número de Reynolds, incluindo a correção quadrática de Oseen.

## Modelo matemático

```math
St\frac{dv^*}{dt^*}=1-v^*-\frac{3}{8}Re_s(v^*)^2
```

com:

```math
v^*(0)=0.
```

Na forma de PVI:

```math
\frac{dv^*}{dt^*}=\frac{1}{St}\left(1-v^*-\frac{3}{8}Re_s(v^*)^2\right).
```

## Método numérico

O método RK4 calcula quatro inclinações por passo:

```math
k_1=f(t_i,v_i),
```

```math
k_2=f\left(t_i+\frac{h}{2},v_i+\frac{h}{2}k_1\right),
```

```math
k_3=f\left(t_i+\frac{h}{2},v_i+\frac{h}{2}k_2\right),
```

```math
k_4=f(t_i+h,v_i+hk_3),
```

```math
v_{i+1}=v_i+\frac{h}{6}(k_1+2k_2+2k_3+k_4).
```

## Arquivos

| Arquivo | Descrição |
|---|---|
| `rk4_sedimentation.py` | código principal já existente no repositório |
| `docs/APC1-PPC1.pdf` | relatório analítico e computacional |
| `exercicios/exercicio-resolvido.md` | primeiro passo do RK4 calculado manualmente |
| `desafio.md` | estudo da ordem de convergência |

## Entradas principais

| Parâmetro | Significado |
|---|---|
| `St` | número de Stokes |
| `Re_s` | número de Reynolds da partícula |
| `h` | passo de tempo adimensional |
| `t_max` | tempo final adimensional |

## Saídas esperadas

- vetor de tempo adimensional;
- velocidade adimensional numérica;
- solução analítica quando `Re_s = 0`;
- erro absoluto;
- gráficos da trajetória e do erro.

## Como executar

```bash
python rk4_sedimentation.py
```

## Validação

Para `Re_s = 0`, a equação é linear e possui solução:

```math
v^*(t^*)=1-e^{-t^*/St}.
```

O erro deve diminuir aproximadamente com `h^4` no regime assintótico do RK4.

## Fontes de erro

- passo de integração finito;
- arredondamento;
- tempo final insuficiente para observar a velocidade terminal;
- uso do modelo de Oseen fora de sua faixa de validade.

## Material complementar

- [Exercício resolvido](exercicios/exercicio-resolvido.md)
- [Desafio](desafio.md)
- [Relatório](docs/APC1-PPC1.pdf)
