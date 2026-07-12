# Desafio — esquema de Crank–Nicolson

Implemente uma segunda opção temporal usando Crank–Nicolson.

## Requisitos

1. manter o mesmo modelo físico e as mesmas condições de contorno;
2. montar o sistema tridiagonal correspondente;
3. resolver com a mesma rotina TDMA;
4. comparar com o esquema totalmente implícito;
5. executar estudo com diferentes `dt` e `N`.

## Comparações mínimas

- erro máximo contra a solução analítica sem geração;
- temperatura máxima no caso com geração;
- tempo de execução;
- sensibilidade ao passo temporal;
- presença de oscilações para passos grandes.

## Resultado esperado

Crank–Nicolson possui segunda ordem no tempo, enquanto o esquema implícito simples possui primeira ordem temporal. Entretanto, passos excessivamente grandes podem produzir comportamento oscilatório em certos problemas.
