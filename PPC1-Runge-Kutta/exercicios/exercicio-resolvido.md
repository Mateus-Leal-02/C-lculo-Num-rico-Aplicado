# Exercício resolvido — um passo de RK4

Considere o caso linear:

```math
St=1,\qquad Re_s=0,\qquad v^*(0)=0
```

com passo `h = 0,5`. A EDO fica:

```math
\frac{dv^*}{dt^*}=1-v^*.
```

## Estágios

```math
k_1=1-0=1
```

```math
k_2=1-\left(0+\frac{0,5}{2}k_1\right)=1-0,25=0,75
```

```math
k_3=1-\left(0+\frac{0,5}{2}k_2\right)=1-0,1875=0,8125
```

```math
k_4=1-(0+0,5k_3)=1-0,40625=0,59375
```

## Atualização

```math
v_1^*=0+\frac{0,5}{6}[1+2(0,75)+2(0,8125)+0,59375]
```

```math
\boxed{v_1^*=0,39322917}
```

A solução exata em `t*=0,5` é:

```math
v_{ex}^*=1-e^{-0,5}=0,39346934.
```

Logo:

```math
E_a=|0,39322917-0,39346934|\approx2,40\times10^{-4}.
```

O cálculo confirma que um único passo relativamente grande já produz boa aproximação, mas a ordem do método deve ser verificada com vários valores de `h`.
