# Roteiro de estudos

## Etapa 1 — Ler o problema

Antes do código, identifique:

1. fenômeno físico ou função matemática;
2. incógnita;
3. domínio;
4. condições iniciais ou de contorno;
5. parâmetros e unidades;
6. resultado de referência disponível.

## Etapa 2 — Acompanhar a progressão dos PPCs

### PPC1 — EDO e integração temporal

Estude a transformação do modelo de sedimentação em um PVI e acompanhe os quatro estágios do RK4.

### PPC2 — Raízes de polinômios

Observe a divisão sintética quadrática, a correção de Newton em dois parâmetros e a deflação do polinômio.

### PPC3 — Sistema tridiagonal em um problema transiente

Relacione a discretização implícita com a estrutura tridiagonal e execute manualmente o algoritmo de Thomas.

### PPC4 — Otimização

Calcule gradiente e Hessiana, interprete a linha de busca e compare direções sucessivas.

### PPC5 — Problema de valor de contorno

Entenda como o Método do Tiro transforma o PVC de Blasius em vários PVIs resolvidos por RK4.

### PPC6 — EDP elíptica

Classifique os nós, derive as equações de diferenças e compare solução direta com métodos iterativos.

## Etapa 3 — Reproduzir o exercício resolvido

Faça o cálculo sem executar o programa. Só depois compare com o arquivo em `exercicios/`.

## Etapa 4 — Executar e variar parâmetros

Altere apenas um parâmetro por vez. Registre:

- valor adotado;
- resultado;
- erro;
- número de iterações;
- tempo de execução;
- interpretação física.

## Etapa 5 — Resolver o desafio

Cada `desafio.md` propõe uma extensão que exige análise, não apenas alteração superficial do código.
