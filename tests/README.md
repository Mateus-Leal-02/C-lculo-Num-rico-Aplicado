# Estratégia futura de testes

Os códigos atuais possuem interfaces diferentes. Antes de criar testes automatizados, recomenda-se transformar cada arquivo principal em funções importáveis sem executar o programa durante o `import`.

## Testes mínimos sugeridos

- PPC1: erro do RK4 para `Re_s = 0` abaixo de uma tolerância definida;
- PPC2: raízes reconstruindo o polinômio original;
- PPC3: TDMA comparado com `numpy.linalg.solve` em um sistema pequeno;
- PPC4: solução convergindo para `(2, 1)`;
- PPC5: `f''(0)` próximo de `0.332057336`;
- PPC6: Gauss, Liebmann e SOR produzindo campos equivalentes dentro da tolerância.

A inclusão de `pytest` deve ser feita somente quando esses testes forem implementados; por isso ele não aparece no `requirements.txt` atual.
