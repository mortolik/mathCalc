#pragma once

#include <vector>
#include <string>
#include <optional>
#include <complex>

struct MatrixResult {
    bool ok{false};
    std::string message;
};

using Complex = std::complex<double>;
using Matrix = std::vector<std::vector<Complex>>;

MatrixResult invertMatrix(const Matrix &input, Matrix &inverse);
MatrixResult eigenvaluesQR(const Matrix &input, std::vector<Complex> &eigs, int maxIter = 200, double tol = 1e-8);
MatrixResult controllabilityRank(const Matrix &A, const Matrix &B, int &rankOut, std::string *witness = nullptr);
MatrixResult observabilityRank(const Matrix &A, const Matrix &C, int &rankOut, std::string *witness = nullptr);
MatrixResult controllabilityMatrix(const Matrix &A, const Matrix &B, Matrix &block);
MatrixResult observabilityMatrix(const Matrix &A, const Matrix &C, Matrix &block);
MatrixResult determinant(const Matrix &input, Complex &detOut);
MatrixResult solveLinearSystem(const Matrix &A, const Matrix &C, Matrix &X, std::vector<size_t> *pivotColsOut = nullptr);
MatrixResult nullspaceVector(const Matrix &A, std::vector<Complex> &v);
int matrixRank(Matrix m, double tol = 1e-9, std::string *witness = nullptr);
bool multiply(const Matrix &A, const Matrix &B, Matrix &out);
bool isSquare(const Matrix &m);
std::optional<size_t> nRows(const Matrix &m);
std::optional<size_t> nCols(const Matrix &m);
