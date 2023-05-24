#include "hw1.h"
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>

namespace algebra {
    Matrix zeros(size_t n, size_t m) {
    // This function will create a n x m matrix with all elements equal to zero.
        std::vector<double> row(m, 0);
        return Matrix(n, row);
    }
    
    Matrix ones(size_t n, size_t m) {
    // This function will create a n x m matrix with all elements equal to one.
        std::vector<double> row(m, 1);
        return Matrix(n, row);
    }

    Matrix random(size_t n, size_t m, double min, double max) {
    // This function will create a n x m matrix with all elements a random number between min and max.
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<double> random_real(min, max);
        Matrix res(n, std::vector<double> (m, 0));
        for (int row=0; row<n; ++row)
            for (int col=0; col<m; ++col)
                res[row][col] = random_real(gen);
        return res;
    }

    void show(const Matrix& matrix) {
    // This function will display the matrix in a beautiful way: each element have 3 decimal places
        std::cout.precision(4);
        for (const auto &row : matrix)
            for (const auto &real_num : row)
                std::cout << real_num << ' ';
            std::cout << std::endl;
        std::cout.precision(7);
    }

    Matrix multiply(const Matrix & matrix, double c) {
    // it multiplies the matrix into the constant scalar c.
        Matrix res = matrix;
        for (auto &row : res)
            for (auto &real_num : row)
                    real_num *= c;
        return res;
    }

    Matrix multiply(const Matrix &matrix1, const Matrix& matrix2) {
    // This function multiplies the matrix1 into matrix2.
        const int n = matrix1.size();
        const int m = matrix2[0].size();
        const int mul_count = matrix1[0].size();
        Matrix res(n, std::vector<double> (m, 0.0));
        for (int row=0; row<n; ++row)
            for (int col=0; col<m; ++col) {
                double sum = 0.0;
                for (int i=0; i<mul_count; ++i) {
                    sum += matrix1[row][i] * matrix2[i][col];
                }
                res[row][col] = sum;
            }
        return res;
    }

    Matrix sum(const Matrix& matrix1, const Matrix& matrix2) {
    // this function it adds 2 matrices to each other.
        const int n = matrix1.size();
        const int m = matrix1[0].size();
        Matrix res(n, std::vector<double> (m, 0.0));
        for (int row=0; row<n; ++row)
            for (int col=0; col<m; ++col) 
                res[row][col] = matrix1[row][col] + matrix2[row][col];
        return res;
    }

    Matrix transpose(const Matrix& matrix) {
    // this function will generate the transpose matrix of the input matrix.
        const int n = matrix.size();
        const int m = matrix[0].size(); 
        Matrix res(m, std::vector<double> (n, 0.0));
        for (int row=0; row<m; ++row)
            for (int col=0; col<n; ++col)
                res[row][col] = matrix[col][row];
        return res;
    }











};