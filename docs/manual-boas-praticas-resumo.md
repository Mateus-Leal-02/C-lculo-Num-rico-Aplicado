# Resumo das boas práticas

## Repositório

- manter o repositório público e autoexplicativo;
- usar uma pasta por PPC;
- manter um README geral como ponto de entrada;
- manter um README específico em cada projeto;
- separar código, entradas, saídas, gráficos e relatórios;
- não depender de caminhos absolutos do computador do autor.

## Código

- comentar detalhadamente a formulação numérica;
- usar nomes de variáveis que representem a grandeza calculada;
- registrar unidades físicas;
- validar entradas;
- evitar duplicação de lógica;
- dividir códigos muito longos em módulos coerentes;
- remover variáveis não utilizadas e soluções improvisadas sem justificativa.

## README interno

Deve conter:

1. objetivo;
2. formulação matemática;
3. algoritmo;
4. dicionário das variáveis;
5. entradas e saídas;
6. dependências;
7. execução;
8. validação;
9. limitações;
10. referências específicas.

## Reprodutibilidade

Uma terceira pessoa deve conseguir:

1. instalar dependências;
2. executar o programa;
3. gerar as saídas;
4. comparar com o resultado de referência;
5. entender o critério de convergência.
