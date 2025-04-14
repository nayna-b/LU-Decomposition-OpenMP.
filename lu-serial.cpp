#include <iostream>
#include <chrono>
#include <cmath>
#include <stdlib.h>

using namespace std;

void print_pivot(int *&, int);
void print_matrix(double **&, int);
void swap_values(int *&, int, int);
void swap_rows(double *&, double *&);
void swap_interval(double **&, int, int);
void initialize_random_matrix(double **&, double **&, int);
void lu_decomposition(double **&, double **&, double **&, int *&, int);
void initialize_for_lu_decomposition(double **&, double **&, double **&, int *&, int);

double l2_norm(double);
double verify_result(int *&, double **&, double **&, double **&, int);

void print_matrix(double **&matrix, int size)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
            cout << matrix[i][j] << " ";
        cout << endl;
    }
    cout << endl;
}

void print_pivot(int *&P, int size)
{
    for (int i = 0; i < size; ++i)
        cout << P[i] << " ";
    cout << endl;
}

void initialize_random_matrix(double **&matrix, double **&matrix_clone, int size)
{
    matrix = new double *[size];
    matrix_clone = new double *[size];
    for (int i = 0; i < size; ++i)
    {
        matrix[i] = new double[size];
        matrix_clone[i] = new double[size];
    }

    unsigned int seed = 221;
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            matrix[i][j] = matrix_clone[i][j] = (rand_r(&seed) % 10000) - 5000;
        }
    }
}

void initialize_for_lu_decomposition(double **&A, double **&L, double **&U, int *&P, int size)
{
    L = new double *[size];
    U = new double *[size];
    P = new int[size];
    for (int i = 0; i < size; ++i)
    {
        L[i] = new double[size];
        U[i] = new double[size];
        P[i] = i;
        for (int j = 0; j < size; ++j)
        {
            L[i][j] = (i == j) ? 1 : 0;
            U[i][j] = 0;
        }
    }
}

void lu_decomposition(double **&A, double **&L, double **&U, int *&P, int size)
{
    initialize_for_lu_decomposition(A, L, U, P, size);

    int new_k = 0;
    double current_max = 0;
    for (int k = 0; k < size; ++k)
    {
        current_max = 0;
        for (int i = k; i < size; ++i)
        {
            if (abs(A[i][k]) > current_max)
            {
                current_max = abs(A[i][k]);
                new_k = i;
            }
        }
        if (current_max == 0)
            throw invalid_argument("Input is a singular matrix");

        if (k != new_k)
        {
            swap_values(P, k, new_k);
            swap_rows(A[k], A[new_k]);
            swap_interval(L, k, new_k);
        }

        U[k][k] = A[k][k];

        for (int i = k + 1; i < size; ++i)
        {
            L[i][k] = A[i][k] / U[k][k];
            U[k][i] = A[k][i];
        }

        for (int i = k + 1; i < size; ++i)
        {
            for (int j = k + 1; j < size; ++j)
                A[i][j] -= L[i][k] * U[k][j];
        }
    }
}

void swap_values(int *&P, int row1, int row2)
{
    int temp = P[row2];
    P[row2] = P[row1];
    P[row1] = temp;
}

void swap_rows(double *&row1, double *&row2)
{
    double *temp = row2;
    row2 = row1;
    row1 = temp;
}

void swap_interval(double **&L, int row, int new_row)
{
    for (int i = 0; i < row; ++i)
    {
        double temp = L[row][i];
        L[row][i] = L[new_row][i];
        L[new_row][i] = temp;
    }
}

double verify_result(int *&P, double **&A, double **&L, double **&U, int size)
{
    double residual = 0;
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            double diff = A[P[i]][j];
            for (int k = 0; k < size; ++k)
                diff -= L[i][k] * U[k][j];
            residual += l2_norm(diff);
        }
    }
    printf("The sum of Euclidean norms is: %e\n", residual);
    return residual;
}

double l2_norm(double value)
{
    return sqrt(value * value);
}

int main(int argc, char **argv)
{
    int matrix_size = atoi(argv[1]);
    int *P;
    double **A, **original_matrix, **L, **U;
    chrono::steady_clock clock;

    initialize_random_matrix(A, original_matrix, matrix_size);

    auto start = clock.now();
    lu_decomposition(A, L, U, P, matrix_size);
    auto end = clock.now();

    auto time_span = static_cast<chrono::duration<double>>(end - start);
    printf("LU Decomposition finished in: %.4f seconds for matrix size %d\n", time_span.count(), matrix_size);

    verify_result(P, original_matrix, L, U, matrix_size);

    return 0;
}
