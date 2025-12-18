#pragma once

#include <vector>
#include <string>
#include <optional>

struct MatrixResult {
    bool ok{false};
    std::string message;
};

using Matrix = std::vector<std::vector<double>>;

MatrixResult invertMatrix(const Matrix &input, Matrix &inverse);
MatrixResult eigenvaluesQR(const Matrix &input, std::vector<double> &eigs, int maxIter = 200, double tol = 1e-8);
MatrixResult controllabilityRank(const Matrix &A, const Matrix &B, int &rankOut);
MatrixResult observabilityRank(const Matrix &A, const Matrix &C, int &rankOut);
MatrixResult determinant(const Matrix &input, double &detOut);
MatrixResult solveLinearSystem(const Matrix &A, const Matrix &C, Matrix &X, std::vector<size_t> *pivotColsOut = nullptr);
MatrixResult nullspaceVector(const Matrix &A, std::vector<double> &v);
int matrixRank(Matrix m, double tol = 1e-9);
bool multiply(const Matrix &A, const Matrix &B, Matrix &out);
bool isSquare(const Matrix &m);
std::optional<size_t> nRows(const Matrix &m);
std::optional<size_t> nCols(const Matrix &m);
