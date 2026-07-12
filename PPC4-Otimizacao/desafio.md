# Desafio — função de Rosenbrock

Adapte os métodos para minimizar:

```math
F(x,y)=(1-x)^2+100(y-x^2)^2
```

partindo de `(-1,2; 1)`.

## Requisitos

1. trocar maximização por minimização;
2. implementar linha de busca numérica;
3. comparar descida máxima e Fletcher–Reeves não linear;
4. registrar número de iterações e avaliações da função;
5. mostrar o caminho sobre curvas de nível;
6. testar pelo menos três pontos iniciais.

## Questões para discutir

- por que o vale estreito dificulta a descida máxima?
- quando é necessário reiniciar a direção conjugada?
- o menor número de iterações também corresponde ao menor custo total?

<details>
<summary>Referência</summary>

O mínimo global é `F(1,1)=0`, mas o desafio está na trajetória e na robustez, não apenas no ponto final.

</details>
