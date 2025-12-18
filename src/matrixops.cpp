#include "matrixops.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

Matrix identity(size_t n) {
    Matrix I(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        I[i][i] = 1.0;
    }
    return I;
}

void swapRows(Matrix &m, size_t i, size_t j) {
    if (i != j) {
        std::swap(m[i], m[j]);
    }
}

} // namespace

bool isSquare(const Matrix &m) {
    if (m.empty()) return false;
    const size_t cols = m.front().size();
    for (const auto &row : m) {
        if (row.size() != cols) return false;
    }
    return cols == m.size();
}

std::optional<size_t> nRows(const Matrix &m) {
    if (m.empty()) return std::nullopt;
    return m.size();
}

std::optional<size_t> nCols(const Matrix &m) {
    if (m.empty()) return std::nullopt;
    return m.front().size();
}

int matrixRank(Matrix m, double tol) {
    if (m.empty()) return 0;
    const size_t rows = m.size();
    const size_t cols = m.front().size();
    size_t r = 0;
    for (size_t c = 0; c < cols && r < rows; ++c) {
        size_t pivot = r;
        double maxVal = std::fabs(m[pivot][c]);
        for (size_t i = r + 1; i < rows; ++i) {
            double val = std::fabs(m[i][c]);
            if (val > maxVal) {
                maxVal = val;
                pivot = i;
            }
        }
        if (maxVal <= tol) continue;
        swapRows(m, pivot, r);
        const double pivotVal = m[r][c];
        for (size_t j = c; j < cols; ++j) {
            m[r][j] /= pivotVal;
        }
        for (size_t i = 0; i < rows; ++i) {
            if (i == r) continue;
            const double factor = m[i][c];
            if (std::fabs(factor) <= tol) continue;
            for (size_t j = c; j < cols; ++j) {
                m[i][j] -= factor * m[r][j];
            }
        }
        ++r;
    }
    return static_cast<int>(r);
}

bool multiply(const Matrix &A, const Matrix &B, Matrix &out) {
    if (A.empty() || B.empty()) return false;
    const size_t aRows = A.size();
    const size_t aCols = A.front().size();
    const size_t bRows = B.size();
    const size_t bCols = B.front().size();
    if (aCols != bRows) return false;
    out.assign(aRows, std::vector<double>(bCols, 0.0));
    for (size_t i = 0; i < aRows; ++i) {
        for (size_t k = 0; k < aCols; ++k) {
            const double aik = A[i][k];
            for (size_t j = 0; j < bCols; ++j) {
                out[i][j] += aik * B[k][j];
            }
        }
    }
    return true;
}

MatrixResult invertMatrix(const Matrix &input, Matrix &inverse) {
    if (!isSquare(input)) {
        return {false, "Matrix must be square"};
    }
    const size_t n = input.size();
    Matrix aug(n, std::vector<double>(2 * n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            aug[i][j] = input[i][j];
        }
        aug[i][n + i] = 1.0;
    }
    const double tol = 1e-12;
    for (size_t c = 0; c < n; ++c) {
        size_t pivot = c;
        double maxVal = std::fabs(aug[pivot][c]);
        for (size_t i = c + 1; i < n; ++i) {
            double val = std::fabs(aug[i][c]);
            if (val > maxVal) {
                maxVal = val;
                pivot = i;
            }
        }
        if (maxVal <= tol) return {false, "Matrix is singular or ill-conditioned"};
        swapRows(aug, pivot, c);
        const double pivotVal = aug[c][c];
        for (size_t j = 0; j < 2 * n; ++j) {
            aug[c][j] /= pivotVal;
        }
        for (size_t i = 0; i < n; ++i) {
            if (i == c) continue;
            const double factor = aug[i][c];
            if (std::fabs(factor) <= tol) continue;
            for (size_t j = 0; j < 2 * n; ++j) {
                aug[i][j] -= factor * aug[c][j];
            }
        }
    }
    inverse.assign(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            inverse[i][j] = aug[i][n + j];
        }
    }
    return {true, "OK"};
}

MatrixResult determinant(const Matrix &input, double &detOut) {
    if (!isSquare(input)) return {false, "Matrix must be square"};
    const size_t n = input.size();
    if (n == 0) return {false, "Matrix is empty"};
    Matrix m = input;
    double det = 1.0;
    const double tol = 1e-12;
    for (size_t c = 0; c < n; ++c) {
        size_t pivot = c;
        double maxVal = std::fabs(m[pivot][c]);
        for (size_t i = c + 1; i < n; ++i) {
            double val = std::fabs(m[i][c]);
            if (val > maxVal) {
                maxVal = val;
                pivot = i;
            }
        }
        if (maxVal <= tol) {
            detOut = 0.0;
            return {true, "OK"};
        }
        if (pivot != c) {
            swapRows(m, pivot, c);
            det = -det; // row swap flips sign
        }
        const double pivotVal = m[c][c];
        det *= pivotVal;
        for (size_t i = c + 1; i < n; ++i) {
            const double factor = m[i][c] / pivotVal;
            for (size_t j = c; j < n; ++j) {
                m[i][j] -= factor * m[c][j];
            }
        }
    }
    detOut = det;
    return {true, "OK"};
}

MatrixResult solveLinearSystem(const Matrix &A, const Matrix &C, Matrix &X, std::vector<size_t> *pivotColsOut) {
    if (!isSquare(A)) return {false, "Matrix must be square"};
    const size_t n = A.size();
    if (n == 0) return {false, "Matrix is empty"};
    const size_t rhsCols = C.empty() ? 0 : C.front().size();
    if (!C.empty()) {
        if (C.size() != n) return {false, "Right-hand side rows must match A"};
        for (const auto &row : C) {
            if (row.size() != rhsCols) return {false, "Right-hand side column counts differ"};
        }
    }
    const size_t m = rhsCols == 0 ? 1 : rhsCols; // allow zero RHS -> treat as zero vector
    Matrix aug(n, std::vector<double>(n + m, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) aug[i][j] = A[i][j];
        for (size_t j = 0; j < m; ++j) aug[i][n + j] = rhsCols == 0 ? 0.0 : C[i][j];
    }
    const double tol = 1e-12;
    // Forward elimination with partial pivoting, track pivots
    std::vector<size_t> pivotRows;
    std::vector<size_t> pivotCols;
    size_t row = 0;
    for (size_t c = 0; c < n && row < n; ++c) {
        size_t pivot = row;
        double maxVal = std::fabs(aug[pivot][c]);
        for (size_t i = row + 1; i < n; ++i) {
            double val = std::fabs(aug[i][c]);
            if (val > maxVal) {
                maxVal = val;
                pivot = i;
            }
        }
        if (maxVal <= tol) continue; // no pivot in this column
        swapRows(aug, pivot, row);
        const double pivotVal = aug[row][c];
        for (size_t i = row + 1; i < n; ++i) {
            const double factor = aug[i][c] / pivotVal;
            for (size_t j = c; j < n + m; ++j) {
                aug[i][j] -= factor * aug[row][j];
            }
        }
        pivotRows.push_back(row);
        pivotCols.push_back(c);
        ++row;
    }

    const size_t rank = pivotRows.size();

    // Check for inconsistency: rows with zero coefficients but non-zero RHS
    for (size_t r = rank; r < n; ++r) {
        bool allZero = true;
        for (size_t c = 0; c < n; ++c) {
            if (std::fabs(aug[r][c]) > tol) { allZero = false; break; }
        }
        if (allZero) {
            for (size_t k = 0; k < m; ++k) {
                if (std::fabs(aug[r][n + k]) > tol) {
                    return {false, "System inconsistent (нет решений)"};
                }
            }
        }
    }

    X.assign(n, std::vector<double>(m, 0.0));
    // Back substitution using pivot rows; free vars remain zero
    for (int idx = static_cast<int>(rank) - 1; idx >= 0; --idx) {
        const size_t r = pivotRows[static_cast<size_t>(idx)];
        const size_t c = pivotCols[static_cast<size_t>(idx)];
        for (size_t k = 0; k < m; ++k) {
            double sum = aug[r][n + k];
            for (size_t j = c + 1; j < n; ++j) {
                sum -= aug[r][j] * X[j][k];
            }
            const double denom = aug[r][c];
            if (std::fabs(denom) <= tol) return {false, "Matrix is singular or degenerate"};
            X[c][k] = sum / denom;
        }
    }

    if (pivotColsOut) {
        *pivotColsOut = pivotCols;
    }

    if (rank < n) {
        return {true, "Бесконечное число решений (свободные переменные = 0)"};
    }
    return {true, "OK"};
}

MatrixResult nullspaceVector(const Matrix &A, std::vector<double> &v) {
    const auto rowsOpt = nRows(A);
    const auto colsOpt = nCols(A);
    if (!rowsOpt || !colsOpt) return {false, "Matrix is empty"};
    const size_t m = *rowsOpt;
    const size_t n = *colsOpt;
    Matrix aug = A;
    aug.resize(m);
    for (auto &row : aug) row.resize(n + 1, 0.0); // append zero RHS

    const double tol = 1e-12;
    std::vector<size_t> pivotRows;
    std::vector<size_t> pivotCols;
    size_t row = 0;
    for (size_t c = 0; c < n && row < m; ++c) {
        size_t pivot = row;
        double maxVal = std::fabs(aug[pivot][c]);
        for (size_t i = row + 1; i < m; ++i) {
            const double val = std::fabs(aug[i][c]);
            if (val > maxVal) {
                maxVal = val;
                pivot = i;
            }
        }
        if (maxVal <= tol) continue;
        swapRows(aug, pivot, row);
        const double pivotVal = aug[row][c];
        for (size_t i = row + 1; i < m; ++i) {
            const double factor = aug[i][c] / pivotVal;
            for (size_t j = c; j < n + 1; ++j) aug[i][j] -= factor * aug[row][j];
        }
        pivotRows.push_back(row);
        pivotCols.push_back(c);
        ++row;
    }

    const size_t rank = pivotRows.size();
    if (rank == n) return {false, "Ядро тривиально"};

    // pick first free column and set it to 1, others 0
    std::vector<bool> isPivot(n, false);
    for (size_t c : pivotCols) isPivot[c] = true;
    size_t freeCol = 0;
    while (freeCol < n && isPivot[freeCol]) ++freeCol;
    if (freeCol == n) return {false, "Не найдена свободная переменная"};

    v.assign(n, 0.0);
    v[freeCol] = 1.0;

    // back substitute
    for (int idx = static_cast<int>(rank) - 1; idx >= 0; --idx) {
        const size_t r = pivotRows[static_cast<size_t>(idx)];
        const size_t c = pivotCols[static_cast<size_t>(idx)];
        double sum = 0.0;
        for (size_t j = c + 1; j < n; ++j) sum += aug[r][j] * v[j];
        const double denom = aug[r][c];
        if (std::fabs(denom) <= tol) return {false, "Matrix is singular"};
        v[c] = -sum / denom;
    }
    return {true, "OK"};
}

static bool qrDecomposition(const Matrix &A, Matrix &Q, Matrix &R) {
    const size_t m = A.size();
    if (m == 0) return false;
    const size_t n = A.front().size();
    Q.assign(m, std::vector<double>(n, 0.0));
    R.assign(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> v = A;
    const double tol = 1e-12;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            double dot = 0.0;
            for (size_t k = 0; k < m; ++k) dot += v[k][i] * Q[k][j];
            R[j][i] = dot;
            for (size_t k = 0; k < m; ++k) v[k][i] -= dot * Q[k][j];
        }
        double norm = 0.0;
        for (size_t k = 0; k < m; ++k) norm += v[k][i] * v[k][i];
        norm = std::sqrt(norm);
        if (norm <= tol) return false;
        R[i][i] = norm;
        for (size_t k = 0; k < m; ++k) Q[k][i] = v[k][i] / norm;
    }
    return true;
}

MatrixResult eigenvaluesQR(const Matrix &input, std::vector<double> &eigs, int maxIter, double tol) {
    if (!isSquare(input)) return {false, "Matrix must be square"};
    const size_t n = input.size();
    if (n == 0) return {false, "Matrix is empty"};
    Matrix A = input;
    for (int iter = 0; iter < maxIter; ++iter) {
        Matrix Q, R;
        if (!qrDecomposition(A, Q, R)) return {false, "QR failed (matrix may be singular)"};
        Matrix next;
        if (!multiply(R, Q, next)) return {false, "Multiply failed"};
        A.swap(next);
        double offDiag = 0.0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i == j) continue;
                offDiag = std::max(offDiag, std::fabs(A[i][j]));
            }
        }
        if (offDiag < tol) break;
    }
    eigs.clear();
    eigs.reserve(n);
    for (size_t i = 0; i < n; ++i) eigs.push_back(A[i][i]);
    return {true, "OK"};
}

MatrixResult controllabilityRank(const Matrix &A, const Matrix &B, int &rankOut) {
    if (!isSquare(A)) return {false, "A must be square"};
    const size_t n = A.size();
    if (B.size() != n) return {false, "B rows must match A"};
    Matrix block = B;
    Matrix power = A;
    Matrix tmp;
    for (size_t i = 1; i < n; ++i) {
        if (!multiply(power, B, tmp)) return {false, "Multiply failed"};
        for (size_t r = 0; r < n; ++r) {
            block[r].insert(block[r].end(), tmp[r].begin(), tmp[r].end());
        }
        Matrix next;
        if (!multiply(power, A, next)) return {false, "Multiply failed"};
        power.swap(next);
    }
    rankOut = matrixRank(block);
    return {true, "OK"};
}

MatrixResult observabilityRank(const Matrix &A, const Matrix &C, int &rankOut) {
    if (!isSquare(A)) return {false, "A must be square"};
    const size_t n = A.size();
    if (C.empty() || C.front().size() != n) return {false, "C cols must match A"};
    Matrix block = C;
    Matrix power = A;
    Matrix tmp;
    for (size_t i = 1; i < n; ++i) {
        if (!multiply(C, power, tmp)) return {false, "Multiply failed"};
        block.insert(block.end(), tmp.begin(), tmp.end());
        Matrix next;
        if (!multiply(power, A, next)) return {false, "Multiply failed"};
        power.swap(next);
    }
    rankOut = matrixRank(block);
    return {true, "OK"};
}

MatrixResult controllabilityMatrix(const Matrix &A, const Matrix &B, Matrix &block) {
    if (!isSquare(A)) return {false, "A must be square"};
    const size_t n = A.size();
    if (B.size() != n) return {false, "B rows must match A"};
    block = B;
    Matrix power = A;
    Matrix tmp;
    for (size_t i = 1; i < n; ++i) {
        if (!multiply(power, B, tmp)) return {false, "Multiply failed"};
        for (size_t r = 0; r < n; ++r) {
            block[r].insert(block[r].end(), tmp[r].begin(), tmp[r].end());
        }
        Matrix next;
        if (!multiply(power, A, next)) return {false, "Multiply failed"};
        power.swap(next);
    }
    return {true, "OK"};
}

MatrixResult observabilityMatrix(const Matrix &A, const Matrix &C, Matrix &block) {
    if (!isSquare(A)) return {false, "A must be square"};
    const size_t n = A.size();
    if (C.empty() || C.front().size() != n) return {false, "C cols must match A"};
    block = C;
    Matrix power = A;
    Matrix tmp;
    for (size_t i = 1; i < n; ++i) {
        if (!multiply(C, power, tmp)) return {false, "Multiply failed"};
        block.insert(block.end(), tmp.begin(), tmp.end());
        Matrix next;
        if (!multiply(power, A, next)) return {false, "Multiply failed"};
        power.swap(next);
    }
    return {true, "OK"};
}
