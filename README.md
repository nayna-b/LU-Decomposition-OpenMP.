# 🔢 Parallel LU Decomposition using OpenMP

This project implements a **shared-memory parallel LU decomposition algorithm** with row pivoting using **OpenMP**. 
The goal is to decompose a large dense matrix `A` into the product of a lower-triangular matrix `L` and an upper-triangular matrix `U` such that:

PA = LU

Where `P` is a permutation matrix that captures row interchanges during pivoting to ensure numerical stability.

----------------------------------------------------------------------------------

## 🚀 Features

- ✅ LU decomposition with **partial pivoting**
- ✅ Parallelized with **OpenMP**
- ✅ Supports **variable matrix sizes and thread counts**
- ✅ **NUMA-aware** optimization (optional)
- ✅ Computes and prints the **L2,1 norm** of the residual matrix `PA - LU` for validation
- ✅ Performance evaluation with scaling across **1 to 32 threads**

----------------------------------------------------------------------------------

## 🛠️ How to Build

Make sure you have an OpenMP-compatible compiler (e.g., Intel `icc` or GCC).

```bash```
make

----------------------------------------------------------------------------------

## ⚙️ Usage

./lu_decompose <matrix_size> <num_threads>

----------------------------------------------------------------------------------

## 📈 Output

Execution Time of LU decomposition phase

L2,1 Norm of the residual matrix PA - LU (should be close to 0)

Optional: Efficiency graphs based on scaling experiments


