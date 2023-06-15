#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MAX 100

double mod(double x) {

    if (x < 0)
        return x * -1;
    return x;
}

void auxCholesky(int ordem, double matriz[MAX][MAX], double matrizCholesky[MAX][MAX]) {
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            // Calcular os elementos da Diagonal Principal
            if (i == j) {
                if (i == 0) {
                    matrizCholesky[i][j] = sqrt(matriz[i][j]);
                } else {
                    double soma = 0;

                    for (int k = 0; k < i; k++) {
                        soma += pow(matrizCholesky[i][k], 2);
                    }
                    matrizCholesky[i][j] = sqrt(matriz[i][j] - soma);
                }
            } else {
                // calcular os elementos fora da diagonal principal
                if (j < i) {
                    double soma = 0;
                    for (int k = 0; k < j; k++) {
                        soma += matrizCholesky[i][k] * matrizCholesky[j][k];
                    }
                    matrizCholesky[i][j] = (matriz[i][j] - soma) / matrizCholesky[j][j];
                } else {
                    matrizCholesky[i][j] = 0;
                }
            }
        }
    }
}

// Função para calcular o determinante
double determinante(int ordem, double matriz[MAX][MAX]) {
    double det = 0;
    double cofator[MAX][MAX];
    int sinal = 1;

    // Caso base para matriz de ordem 1
    if (ordem == 1) {
        return matriz[0][0];
    }

    // Loop para percorrer a primeira linha da matriz
    for (int i = 0; i < ordem; i++) {
        int submatrizLinha = 0;
        int submatrizColuna = 0;

        // Criação da matriz cofatora
        for (int linha = 1; linha < ordem; linha++) {
            for (int coluna = 0; coluna < ordem; coluna++) {
                if (coluna != i) {
                    cofator[submatrizLinha][submatrizColuna] = matriz[linha][coluna];
                    submatrizColuna++;

                    // Incrementa a coluna da submatriz
                    if (submatrizColuna == ordem - 1) {
                        submatrizColuna = 0;
                        submatrizLinha++;
                    }
                }
            }
        }

        // Chamada recursiva para calcular o determinante da submatriz
        det += sinal * matriz[0][i] * determinante(ordem - 1, cofator);

        // Alterna o sinal para o próximo cofator
        sinal = -sinal;
    }

    return det;
}

// Função para convergência dos AKs

int convergenciaAKZero(int ordem, double matriz[MAX][MAX]) {

    for (int i = 1; i <= ordem; i++) {
        // teste
        printf("Determinante A K=%d: %.4lf\n", i, determinante(i, matriz));

        if (!determinante(i, matriz)) { // se devolver 0

            // teste:
            printf("\nNao converge\n");
            return 0; // não converge
        }
    }
    // teste:
    printf("\nConverge!\n");
    return 1;
}

// Função para calcular a matriz triangular inferior

int convergenciaAKMaior(int ordem, double matriz[MAX][MAX]) {
    for (int i = 1; i <= ordem; i++) {
        // teste
        printf("Determinante A K=%d: %.4lf", i, determinante(i, matriz));

        if (determinante(i, matriz) < 0) { // se devolver 0

            // teste:
            printf("\nNao converge\n");
            return 0; // não converge
        }
    }
    // teste:
    printf("\nConverge!\n");
    return 1;
}

// Função para calcular o sistema de matriz triangular inferior

void SistemaTriangularInferior(int ordem, double coeficientes[MAX][MAX], double vetorIndInf[], double *vetorSolInf) {
    vetorSolInf[0] = vetorIndInf[0] / coeficientes[0][0];

    for (int i = 1; i < ordem; i++) {
        double soma = 0;
        for (int j = 0; j < i; j++) {
            soma += coeficientes[i][j] * vetorSolInf[j];
        }
        vetorSolInf[i] = (vetorIndInf[i] - soma) / coeficientes[i][i];
    }
}

// Função para calcular o sistema de matriz triangular superior

void SistemaTriangularSuperior(int ordem, double coeficientes[MAX][MAX], double vetorIndSup[], double *vetorSolSup) {
    vetorSolSup[ordem - 1] = vetorIndSup[ordem - 1] / coeficientes[ordem - 1][ordem - 1];
    // Mostrar o vetor Independente
    printf("\nVetor Independente: \n");
    for (int i = 0; i < ordem; i++) {
        printf("%.4lf\n", vetorIndSup[i]);
    }

    for (int i = ordem - 2; i >= 0; i--) {
        double soma = 0;
        for (int j = i + 1; j <= ordem - 1; j++) {
            soma += coeficientes[i][j] * vetorSolSup[j];
        }
        vetorSolSup[i] = (vetorIndSup[i] - soma) / coeficientes[i][i];
    }
}

// Método de decomposição LU
void DecomposicaoLU(int ordem, double coeficientes[MAX][MAX], double *vetorInd, double *vetorSol) {

    // VERIFICAR CONVERGENCIA

    if (!convergenciaAKZero(ordem, coeficientes)) { // se falso
        printf("Nao converge");
        return;
    }

    // se converge:

    double matrizU[MAX][MAX];
    double matrizL[MAX][MAX];
    double somaU, somaL;
    double vetorY[MAX];

    // incializacao da matriz L:
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            if (i == j) {
                matrizL[i][j] = 1;

            } else {

                matrizL[i][j] = 0;
            }
        }
    }
    // incializacao da matriz U:
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            matrizU[i][j] = 0;
        }
    }
    // teste Matriz L:
    printf("\nMatriz L:\n");
    for (int i = 0; i < ordem; i++) {
        printf("\n");
        for (int j = 0; j < ordem; j++) {
            // teste:
            printf("\t%.4f", matrizL[i][j]);
        }
    }
    // teste Matriz U:
    printf("\nMatriz U:\n");
    for (int i = 0; i < ordem; i++) {
        printf("\n");
        for (int j = 0; j < ordem; j++) {
            printf("\t%.4f", matrizU[i][j]);
        }
    }

    // obter matrizes U e L

    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            if (i <= j) { // calcule U
                for (int k = 0; k <= i - 1; k++) {
                    somaU = matrizL[i][k] * matrizU[k][j];
                }

                matrizU[i][j] = coeficientes[i][j] - somaU;

            } else if (i > j) { // calcule L
                for (int k = 0; k <= j - 1; k++) {
                    somaL = matrizL[i][k] * matrizU[k][j];
                }
                matrizL[i][j] = (coeficientes[i][j] - somaL) / matrizU[j][j];
            }
        }
    }

    // teste Matriz L:
    printf("\nMatriz L:\n");
    for (int i = 0; i < ordem; i++) {
        printf("\n");
        for (int j = 0; j < ordem; j++) {
            // teste:
            printf("\t%.4f", matrizL[i][j]);
        }
    }
    printf("\n");
    // teste Matriz U:
    printf("\nMatriz U:\n");
    for (int i = 0; i < ordem; i++) {
        printf("\n");
        for (int j = 0; j < ordem; j++) {
            printf("\t%.4f", matrizU[i][j]);
        }
    }
    printf("\n");

    // solucionar Ly = B

    SistemaTriangularInferior(ordem, matrizL, vetorInd, vetorY);
    // solucionar Ux = Y

    SistemaTriangularSuperior(ordem, matrizU, vetorY, vetorSol);
}

// Rotina para resolver sistema por Gauss Jordan
void GaussJordan(int ordem, double coeficientes[MAX][MAX], double vetorIndJordan[], double *vetorSolJordan) {
    // Verifica se o sistema converge
    if (convergenciaAKZero(ordem, coeficientes) == 0) {
        printf("Sistema nao converge!\n");
        return;
    }
    // Faz a matriz expandida
    double coeficientesExp[MAX][MAX + 1];
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j <= ordem; j++) {
            if (j == ordem) {
                coeficientesExp[i][j] = vetorIndJordan[i];
            } else {
                coeficientesExp[i][j] = coeficientes[i][j];
            }
        }
    }
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            if (i == j) {
                continue;
            } else {
                double fator = coeficientesExp[j][i] / coeficientesExp[i][i];
                for (int k = 0; k <= ordem; k++) {
                    coeficientesExp[j][k] -= fator * coeficientesExp[i][k];
                }
            }
        }
    }
    printf("Matriz escalonada: \n");
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem + 1; j++) {
            printf("%.4f ", coeficientesExp[i][j]);
        }
        printf("\n");
    }
    // Depois de feitas as operações com a Matriz, o vetor solução é a última coluna da matriz
    for (int i = 0; i < ordem; i++) {
        vetorSolJordan[i] = coeficientesExp[i][ordem] / coeficientesExp[i][i];
    }
}

// Rotina para resolver sistema por Cholesky
void Cholesky(int ordem, double coeficientes[MAX][MAX], double vetorIndCholesky[], double *vetorSolCholesky) {
    double matrizL[MAX][MAX];
    auxCholesky(ordem, coeficientes, matrizL);
    printf("Matriz L: \n");
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            printf("%.4f ", matrizL[i][j]);
        }
        printf("\n");
    }
    // Resolvendo o sistema Ly = b utilizando a função de matriz triangular inferior
    double vetorY[MAX];
    SistemaTriangularInferior(ordem, matrizL, vetorIndCholesky, vetorY);
    // Fazer a matriz transposta de L
    double matrizLT[MAX][MAX];
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            matrizLT[i][j] = matrizL[j][i];
        }
    }
    // Resolvendo o sistema L^t x = y utilizando a função de matriz triangular superior
    SistemaTriangularSuperior(ordem, matrizLT, vetorY, vetorSolCholesky);
}

// Rotina para resolver sistema utilizando Gauss Compacto
void GaussCompacto(int ordem, double coeficientes[MAX][MAX], double vetorIndCompacto[], double *vetorSolCompacto) {
    // Verifica se o sistema converge
    if (convergenciaAKZero(ordem, coeficientes) == 0) {
        printf("Sistema nao converge!\n");
        return;
    }
    // Faz a matriz expandida
    double coeficientesExp[MAX][MAX + 1];
    // Montar a matriz expandidda com a matriz de coeficientes e o vetor independente
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j <= ordem; j++) {
            if (j == ordem) {
                coeficientesExp[i][j] = vetorIndCompacto[i];
            } else {
                coeficientesExp[i][j] = coeficientes[i][j];
            }
        }
    }

    double matrizLU[MAX][MAX];
    // Determinar os elementos L e U da matrizLU
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem + 1; j++) {
            double soma = 0;
            if (i <= j) {
                for (int k = 0; k <= i - 1; k++) {
                    soma += matrizLU[i][k] * matrizLU[k][j];
                }
                matrizLU[i][j] = coeficientesExp[i][j] - soma;
            } else {
                for (int k = 0; k <= j - 1; k++) {
                    soma += matrizLU[i][k] * matrizLU[k][j];
                }
                matrizLU[i][j] = (coeficientesExp[i][j] - soma) / matrizLU[j][j];
            }
        }
    }
    // Zera a parte inferior da matrizLU
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < i; j++) {
            matrizLU[i][j] = 0;
        }
    }
    // Imprimir a matriz LU
    printf("Matriz LU: \n");
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem + 1; j++) {
            printf("%.4f ", matrizLU[i][j]);
        }
        printf("\n");
    }
    // modifica Vetor independente
    for (int i = 0; i < ordem; i++) {
        vetorIndCompacto[i] = matrizLU[i][ordem];
        printf("Adicionado: %.4f\n, ao vetor independente ", matrizLU[i][ordem]);
    }
    printf("Vetor independente modificado: \n");
    for (int i = 0; i < ordem; i++) {
        printf("\n%.4f\n ", vetorIndCompacto[i]);
    }
    // Resolvendo o sistema Ly = b utilizando a função de matriz triangular superior
    SistemaTriangularSuperior(ordem, matrizLU, vetorIndCompacto, vetorSolCompacto);
}

void MatrizInversa(int ordem, double matriz[MAX][MAX], double matrizInv[MAX][MAX]) {
    system("clear");
    system("cls");
    printf("Qual método deseja utilizar para calcular a matriz inversa?\n");
    printf("1 - Gauss Compacto\n");
    printf("2 - Decomp. LU\n");
    int opcao;
    scanf("%d", &opcao);
    switch (opcao) {
        case 1:;
            // Calcular matriz inversa utilizando Gauss Compacto, utilizando n vetores indepentes (n = ordem)

            for (int i = 0; i < ordem; i++) {
                double vetorInv[MAX];
                double vetorIndepnete[MAX];
                for (int j = 0; j < ordem; j++) {
                    if (j == i) {
                        vetorIndepnete[j] = 1;
                    } else {
                        vetorIndepnete[j] = 0;
                    }
                }
                GaussCompacto(ordem, matriz, vetorIndepnete, vetorInv);
                for (int j = 0; j < ordem; j++) {
                    matrizInv[j][i] = vetorInv[j];
                }
            }

            break;
        case 2:;
            for (int i = 0; i < ordem; i++) {
                double vetorInv[MAX];
                double vetorIndepnete[MAX];
                for (int j = 0; j < ordem; j++) {
                    if (j == i) {
                        vetorIndepnete[j] = 1;
                    } else {
                        vetorIndepnete[j] = 0;
                    }
                }
                DecomposicaoLU(ordem, matriz, vetorIndepnete, vetorInv);
                for (int j = 0; j < ordem; j++) {
                    matrizInv[j][i] = vetorInv[j];
                }
            }
            break;
        default:
            printf("Opcao invalida!\n");
            break;
    }
}

// Rotina utilizando método de Jacobi

int Jacobi(int ord, double matriz[MAX][MAX], double *vet_indep, double *vet_apr, double erro, int max_it, double *sol, int *n_iter) {

    enum {
        reqa,
        reqb,
        reqc,
        suc
    };

    int check = 0;

    for (int i = 0; i < ord; i++) {
        if (matriz[i][i] == 0)
            return reqa;
    }

    if (!determinante(ord, matriz))
        return reqb;

    double sum = 0, max = 0;

    for (int i = 0; i < ord; i++) {

        for (int j = 0; j < ord; j++) {

            if (j == i)
                ;
            else {
                sum += mod(matriz[i][j]);
            }
        }

        if (mod(sum / matriz[i][i]) > max)
            max = mod(sum / matriz[i][i]);

        sum = 0;
    }

    if (max >= 1) {

        max = 0;

        for (int i = 0; i < ord; i++) {

            for (int j = 0; j < ord; j++) {

                if (j == i)
                    j++;
                else
                    sum += mod(matriz[j][i]);
            }

            if (mod(sum / matriz[i][i]) > max)
                max = mod(sum / matriz[i][i]);

            sum = 0;
        }
    }

    if (max >= 1)
        return reqc;

    double oldsol[ord], err_at = 10;

    sum = 0;

    for (int i = 0; i < ord; i++)
        oldsol[i] = vet_apr[i];

    do {

        *n_iter += 1;

        for (int i = 0; i < ord; i++) {

            for (int j = 0; j < ord; j++) {

                if (j == i)
                    ;
                else {

                    sum += matriz[i][j] * oldsol[j];
                }
            }

            sol[i] = (vet_indep[i] - sum) / matriz[i][i];

            sum = 0;
        }

        err_at = 0;

        for (int i = 0; i < ord; i++) {
            if (mod(sol[i] - oldsol[i]) > err_at)
                err_at = (mod(sol[i] - oldsol[i])) / mod(sol[i]);
        }

        for (int i = 0; i < ord; i++)
            oldsol[i] = sol[i];

    } while ((err_at > erro) && *n_iter < max_it);
    return suc;
}

// Rotina de Gauss Seidel

int GaussSeidel(int ord, double matriz[MAX][MAX], double *vet_indep, double *vet_apr, double erro, int max_it, double *sol, int *n_iter) {

    enum {
        reqa,
        reqb,
        reqc,
        suc
    };

    int check = 0;

    for (int i = 0; i < ord; i++) {
        if (matriz[i][i] == 0)
            return reqa;
    }

    if (!determinante(ord, matriz))
        return reqb;

    double sum = 0, max = 0;

    for (int i = 0; i < ord; i++) {

        for (int j = 0; j < ord; j++) {

            if (j == i)
                ;
            else {
                sum += mod(matriz[i][j]);
            }
        }

        if (mod(sum / matriz[i][i]) > max)
            max = mod(sum / matriz[i][i]);

        sum = 0;
    }

    if (max >= 1) {

        max = 0;

        double b[ord];

        for (int i = 0; i < ord; i++) {

            for (int j = 0; j < i - 1; j++)
                sum += mod(matriz[i][j]) * b[j];
            for (int j = i + 1; j < ord; j++)
                sum += mod(matriz[i][j]);

            b[i] = sum / matriz[i][i];

            if (b[i] > max)
                max = b[i];

            sum = 0;
        }
    }

    if (max >= 1)
        return reqc;

    double oldsol[ord], err_at = 0;

    sum = 0;

    for (int i = 0; i < ord; i++)
        oldsol[i] = vet_apr[i];

    do {

        *n_iter += 1;

        for (int i = 0; i < ord; i++) {

            for (int j = 0; j < ord; j++) {

                if (j == i)
                    ;
                else {

                    if (j > i) {

                        sum += matriz[i][j] * oldsol[j];

                    } else {

                        sum += matriz[i][j] * sol[j];
                    }
                }
            }

            sol[i] = (vet_indep[i] - sum) / matriz[i][i];

            sum = 0;
        }

        err_at = 0;

        for (int i = 0; i < ord; i++) {
            if (mod(sol[i] - oldsol[i]) > err_at)
                err_at = (mod(sol[i] - oldsol[i])) / mod(sol[i]);
        }

        for (int i = 0; i < ord; i++) {
            printf("X%d = %lf", i, sol[i]);
            printf("\n");
        }
        printf("Num de it. = %d\n\n", *n_iter);

        for (int i = 0; i < ord; i++)
            oldsol[i] = sol[i];

    } while ((err_at > erro) && *n_iter < max_it);
    return suc;
}

int main() {
    int exit = 0;
    do {
        double matriz[MAX][MAX];
        int ordem;
        printf("Digite a ordem da matriz: ");
        scanf("%d", &ordem);

        printf("Digite os elementos da matriz (Linha por linha):\n");
        for (int i = 0; i < ordem; i++) {
            for (int j = 0; j < ordem; j++) {
                scanf("%lf", &matriz[i][j]);
            }
        }

        printf("Digite o metodo de resolucao desejado: \n");
        printf("1 - Determinante\n");
        printf("2 - Triangular Inferior\n");
        printf("3 - Triangular Superior\n");
        printf("4 - Decomposicao LU\n");
        printf("5 - Rotina Cholesky\n");
        printf("6 - Rotina Gauss-Compacto\n");
        printf("7 - Rotina Gauss-Jordan\n");
        printf("8 - Rotina Jacobi\n");
        printf("9 - Rotina Gauss-Seidel\n");
        printf("10 - Rotina Matriz Inversa\n");
        printf("11 - Para sair do programa\n");

        int opcao;

        scanf("%d", &opcao);

        switch (opcao) {
            case 1:
                printf("Determinante: %.4f\n", determinante(ordem, matriz));
                break;
            case 2:;
                double vetorIndInf[MAX];
                double vetorSolInf[MAX];
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndInf[i]);
                }
                SistemaTriangularInferior(ordem, matriz, vetorIndInf, vetorSolInf);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolInf[i]);
                }
                break;
            case 3:;
                double vetorIndSup[MAX];
                double vetorSolSup[MAX];
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndSup[i]);
                }
                SistemaTriangularSuperior(ordem, matriz, vetorIndSup, vetorSolSup);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolSup[i]);
                }
                break;
            case 4:;
                double vetorIndLU[MAX];
                double vetorSolLU[MAX];
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndLU[i]);
                }
                DecomposicaoLU(ordem, matriz, vetorIndLU, vetorSolLU);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolLU[i]);
                }
                break;
            case 5:;
                double vetorIndCholesky[MAX];
                double vetorSolCholesky[MAX];
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndCholesky[i]);
                }
                Cholesky(ordem, matriz, vetorIndCholesky, vetorSolCholesky);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolCholesky[i]);
                }
                break;
            case 6:;
                double vetorIndGaussCompacto[MAX];
                double vetorSolGaussCompacto[MAX];
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndGaussCompacto[i]);
                }
                GaussCompacto(ordem, matriz, vetorIndGaussCompacto, vetorSolGaussCompacto);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolGaussCompacto[i]);
                }
                break;
            case 7:;
                double vetorIndJordan[MAX];
                double vetorSolJordan[MAX];
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndJordan[i]);
                }
                GaussJordan(ordem, matriz, vetorIndJordan, vetorSolJordan);
                printf("Vetor Solcao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolJordan[i]);
                }
                break;
            case 8:;
                double vetorIndJacobi[MAX];
                double vetorSolJacobi[MAX];
                double vetorAproxJacobi[MAX];
                double erro;
                int max_it;
                int n_iterJacobi;
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndJacobi[i]);
                }
                printf("Digite o vetor de aproximacao inicial para a solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorAproxJacobi[i]);
                }
                printf("Digite a precisao desejada (erro): \n");
                scanf("%lf", &erro);
                printf("Digite o número máximo de iterações: \n");
                scanf("%d", &max_it);
                Jacobi(ordem, matriz, vetorIndJacobi, vetorAproxJacobi, erro, max_it, vetorSolJacobi, &n_iterJacobi);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolJacobi[i]);
                }
                break;
            case 9:;
                double vetorIndSeidel[MAX];
                double vetorSolSeidel[MAX];
                double vetorAproxSeidel[MAX];
                int n_iterSeidel;
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndSeidel[i]);
                }
                printf("Digite o vetor de aproximacao inicial para a solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorAproxSeidel[i]);
                }
                printf("Digite o número máximo de iterações: \n");
                scanf("%d", &max_it);
                GaussSeidel(ordem, matriz, vetorIndSeidel, vetorAproxSeidel, erro, max_it, vetorSolSeidel, &n_iterSeidel);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolSeidel[i]);
                }
                break;
            case 10:;
                double matrizInv[MAX][MAX];
                MatrizInversa(ordem, matriz, matrizInv);
                printf("Matriz Inversa: \n");
                for (int i = 0; i < ordem; i++) {
                    for (int j = 0; j < ordem; j++) {
                        printf("%.4f ", matrizInv[i][j]);
                    }
                    printf("\n");
                }
                break;
            case 11:
                exit = 1;
                break;
            default:
                printf("Opcao invalida!\n");
                break;
        }
    } while (exit != 1);

    return 0;
}