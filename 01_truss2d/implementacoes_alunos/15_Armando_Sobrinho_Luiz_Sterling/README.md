# Portico 2D com efeito P-Delta

Implementacao: Armando Sobrinho e Luiz Sterling

Origem do trabalho: modulo `01_truss2d`, a partir do programa de trelica 2D em Fortran 77.

## Descricao

Esta implementacao adapta o programa inicial de trelica 2D para uma formulacao de portico plano 2D. O elemento original de barra axial, com dois graus de liberdade por no, foi substituido por um elemento de portico plano com tres graus de liberdade por no: deslocamento horizontal `ux`, deslocamento vertical `uy` e rotacao `rz`.

O programa monta a matriz de rigidez elastica de um elemento de portico plano de dois nos, com contribuicoes axial e flexional. Tambem foi adicionada uma rotina iterativa para considerar o efeito `P-Delta` por meio da matriz de rigidez geometrica calculada a partir do esforco normal de cada elemento.

## Formulacao resumida

Para cada no:

```text
gdl = [ux, uy, rz]
```

Para cada elemento:

```text
ue = [u1, v1, theta1, u2, v2, theta2]
```

A rigidez total do elemento em coordenadas locais e:

```text
kt = ke + kg
```

onde `ke` e a matriz elastica do elemento de portico plano e `kg` e a matriz geometrica associada ao esforco normal axial `N`.

A transformacao local-global e feita por:

```text
Ke_global = T^T * kt * T
```

O efeito `P-Delta` e resolvido iterativamente:

1. inicia-se com `N = 0`;
2. monta-se a rigidez global;
3. resolve-se o sistema;
4. atualizam-se os esforcos normais;
5. repete-se ate a convergencia dos deslocamentos.

## Arquivos

```text
src/portico2d_pdelta.f      versao principal em formato fixo
src/portico2d_pdelta.f90    versao em formato livre para consulta
exemplos/portico2d_pdelta.dat
exemplos/portico2d_pdelta.out
```

## Como compilar

No Windows com MSYS2/gfortran:

```powershell
$env:PATH='C:\msys64\mingw64\bin;' + $env:PATH
gfortran -Wall src\portico2d_pdelta.f -o portico2d_pdelta.exe
```

## Como rodar

O programa le o arquivo `portico2d_pdelta.dat` no diretorio atual. Uma forma simples de executar o exemplo e:

```powershell
Copy-Item exemplos\portico2d_pdelta.dat .
.\portico2d_pdelta.exe
```

A saida sera gerada em:

```text
portico2d_pdelta.out
```

## Formato do arquivo de entrada

```text
NN NE NM NBC NLOAD
METHOD
IPDELTA MAXITER TOL
nos:        id x y
materiais:  id E A I
elementos:  id no1 no2 matid
restricoes: no dof valor
cargas:     no dof valor
```

Graus de liberdade:

```text
1 = ux
2 = uy
3 = rz
```

Opcoes:

```text
METHOD = 1  penalidade
METHOD = 2  eliminacao de linhas/colunas

IPDELTA = 0 analise linear
IPDELTA = 1 analise P-Delta
```

## Exemplo de validacao

O exemplo incluido representa um portico plano com tres elementos e engastamento apenas no no 1. A analise usa `P-Delta` ligado:

```text
IPDELTA = 1
MAXITER = 30
TOL = 1.0D-8
```

Na ultima verificacao, o programa convergiu em 3 iteracoes. A comparacao com resultados de referencia fornecidos pelo usuario apresentou erros relativos inferiores a 0.75% para os deslocamentos e rotacoes avaliados.

## Referencias

- Cook et al. - Concepts and Applications of Finite Element Analysis.
- Bathe - Finite Element Procedures.
- Zienkiewicz e Taylor - The Finite Element Method.

