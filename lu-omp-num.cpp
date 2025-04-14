#include <iostream>
#include <chrono>
#include <cmath>
#include <omp.h>
#include <stdlib.h>
#include <numa.h>
#include <unistd.h>
#include <random>

using namespace std;

// Function declarations
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

// Utility function to print matrix
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

// Utility function to print pivot vector
void print_pivot(int *&P, int size)
{
    for (int i = 0; i < size; ++i)
        cout << P[i] << " ";
    cout << endl;
}

// Initializes a matrix and its clone with random values in [-1, 1] using NUMA-aware allocation
void initialize_random_matrix(double **&matrix, double **&matrix_clone, int size)
{
    matrix = (double **)malloc(sizeof(double *) * size);
    matrix_clone = (double **)malloc(sizeof(double *) * size);

#pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        matrix[i] = (double *)numa_alloc_onnode(sizeof(double) * size, i < size / 2 ? 0 : 1);
        matrix_clone[i] = (double *)numa_alloc_onnode(sizeof(double) * size, i < size / 2 ? 0 : 1);
    }

    // Thread-safe random generation
#pragma omp parallel
    {
        struct drand48_data buffer;
        srand48_r(221 * (omp_get_thread_num() + 1), &buffer);

#pragma omp for schedule(static)
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                double val;
                drand48_r(&buffer, &val);
                matrix[i][j] = matrix_clone[i][j] = val * 2.0 - 1.0;
            }
        }
    }
}

// Allocates and initializes matrices L and U, and the pivot array P
void initialize_for_lu_decomposition(double **&A, double **&L, double **&U, int *&P, int size)
{
    L = (double **)numa_alloc_onnode(sizeof(double *) * size, 0);
    U = (double **)numa_alloc_onnode(sizeof(double *) * size, 1);
    P = new int[size];

#pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        L[i] = (double *)numa_alloc_onnode(sizeof(double) * size, i < size / 2 ? 0 : 1);
        U[i] = (double *)numa_alloc_onnode(sizeof(double) * size, i < size / 2 ? 0 : 1);
        P[i] = i;
    }

    // Fill L with identity matrix and U with zeros
#pragma omp parallel for schedule(static)
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            L[i][j] = (i == j) ? 1 : 0;
            U[i][j] = 0;
        }
    }
}

// Main LU Decomposition function with partial pivoting
void lu_decomposition(double **&A, double **&L, double **&U, int *&P, int size)
{
    initialize_for_lu_decomposition(A, L, U, P, size);
    int new_k = 0;
    double current_max = 0;

    for (int k = 0; k < size; ++k)
    {
        // Find pivot for numerical stability
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

        // Perform row swaps if needed
        if (k != new_k)
        {
            swap_values(P, k, new_k);
            swap_rows(A[k], A[new_k]);
            swap_interval(L, k, new_k);
        }

        U[k][k] = A[k][k];

        // Fill L and U
#pragma omp parallel for schedule(static)
        for (int i = k + 1; i < size; ++i)
        {
            L[i][k] = A[i][k] / U[k][k];
            U[k][i] = A[k][i];
        }

        // Gaussian elimination step
#pragma omp parallel for schedule(static)
        for (int i = k + 1; i < size; ++i)
        {
            for (int j = k + 1; j < size; ++j)
                A[i][j] -= L[i][k] * U[k][j];
        }
    }
}

// Swaps values in pivot array
void swap_values(int *&P, int row1, int row2)
{
    int tmp = P[row2];
    P[row2] = P[row1];
    P[row1] = tmp;
}

// Swaps entire rows in a matrix
void swap_rows(double *&row1, double *&row2)
{
    double *tmp = row2;
    row2 = row1;
    row1 = tmp;
}

// Swaps values within a row interval in matrix L
void swap_interval(double **&L, int row, int new_row)
{
    for (int i = 0; i < row; ++i)
    {
        double tmp = L[row][i];
        L[row][i] = L[new_row][i];
        L[new_row][i] = tmp;
    }
}

// Verifies that PA = LU by computing the residual
double verify_result(int *&P, double **&A, double **&L, double **&U, int size)
{
    double residual = 0.0;
#pragma omp parallel for reduction(+ : residual) schedule(static)
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            double lu_val = 0.0;
            for (int k = 0; k < size; ++k)
                lu_val += L[i][k] * U[k][j];
            double diff = A[P[i]][j] - lu_val;
            residual += l2_norm(diff);
        }
    }
    printf("The sum of Euclidean norms is: %e\n", residual);
    return residual;
}

// Computes l2 norm of a value
double l2_norm(double value)
{
    return sqrt(value * value);
}

// Entry point for the program: reads input, runs LU, prints timing and residual
int main(int argc, char **argv)
{
    int matrix_size = atoi(argv[1]), nworkers = atoi(argv[2]);
    int *P;
    double **A, **original_matrix, **L, **U;
    chrono::steady_clock clock;

    omp_set_num_threads(nworkers);
    initialize_random_matrix(A, original_matrix, matrix_size);

    auto start = clock.now();
    lu_decomposition(A, L, U, P, matrix_size);
    auto end = clock.now();

    auto time_span = static_cast<chrono::duration<double>>(end - start);
    printf("LU Decomposition finished in: %.4f seconds for %d workers\n", time_span.count(), nworkers);

    verify_result(P, original_matrix, L, U, matrix_size);

    // Free all NUMA-allocated and heap memory
    for (int i = 0; i < matrix_size; i++)
    {
        numa_free(A[i], matrix_size * sizeof(double));
        numa_free(original_matrix[i], matrix_size * sizeof(double));
        numa_free(L[i], matrix_size * sizeof(double));
        numa_free(U[i], matrix_size * sizeof(double));
    }
    free(A);
    free(original_matrix);
    numa_free(L, matrix_size * sizeof(double *));
    numa_free(U, matrix_size * sizeof(double *));
    delete[] P;

    return 0;
}
