#include "matrixops.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <sstream>

namespace {

Matrix identity(size_t n) {
    Matrix I(n, std::vector<Complex>(n, Complex{0.0, 0.0}));
    for (size_t i = 0; i < n; ++i) {
        I[i][i] = Complex{1.0, 0.0};
    }
    return I;
}

void swapRows(Matrix &m, size_t i, size_t j) {
    if (i != j) {
        std::swap(m[i], m[j]);
    }
}

std::string formatComplex(const Complex &z, double tol = 1e-9) {
    std::ostringstream os;
    os.setf(std::ios::fixed);
    os.precision(6);
    const double r = z.real();
    const double im = z.imag();
    if (std::abs(im) <= tol) {
        os << r;
    } else if (std::abs(r) <= tol) {
        os << im << 'i';
    } else {
        os << r << (im >= 0 ? "+" : "") << im << 'i';
    }
    return os.str();
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

int matrixRank(Matrix m, double tol, std::string *witness) {
    if (m.empty()) return 0;
    if (tol <= 0) tol = 1e-9;
    const size_t rows = m.size();
    const size_t cols = m.front().size();
    const size_t maxK = std::min(rows, cols);

    auto next_combination = [](std::vector<size_t> &comb, size_t n) {
        const size_t k = comb.size();
        for (size_t i = k; i-- > 0;) {
            if (comb[i] < n - k + i) {
                ++comb[i];
                for (size_t j = i + 1; j < k; ++j) comb[j] = comb[j - 1] + 1;
                return true;
            }
        }
        return false;
    };

    Matrix sub;
    if (witness) witness->clear();
    // Check minors from largest to smallest until a non-zero determinant is found
    for (size_t k = maxK; k >= 1; --k) {
        std::vector<size_t> rowIdx(k), colIdx(k);
        std::iota(rowIdx.begin(), rowIdx.end(), 0);
        std::iota(colIdx.begin(), colIdx.end(), 0);
        bool moreRows = true;
        while (moreRows) {
            bool moreCols = true;
            while (moreCols) {
                sub.assign(k, std::vector<Complex>(k, Complex{0.0, 0.0}));
                for (size_t i = 0; i < k; ++i) {
                    for (size_t j = 0; j < k; ++j) {
                        sub[i][j] = m[rowIdx[i]][colIdx[j]];
                    }
                }
                Complex det = Complex{0.0, 0.0};
                auto res = determinant(sub, det);
                if (res.ok && std::abs(det) > tol) {
                    if (witness) {
                        std::ostringstream os;
                        os.setf(std::ios::fixed);
                        os.precision(6);
                        os << "Ненулевой минор " << k << "x" << k << " (строки: ";
                        for (size_t idx = 0; idx < k; ++idx) {
                            os << (rowIdx[idx] + 1);
                            if (idx + 1 < k) os << ',';
                        }
                        os << "; столбцы: ";
                        for (size_t idx = 0; idx < k; ++idx) {
                            os << (colIdx[idx] + 1);
                            if (idx + 1 < k) os << ',';
                        }
                        os << ") det = " << formatComplex(det);
                        *witness = os.str();
                    }
                    return static_cast<int>(k);
                }
                moreCols = next_combination(colIdx, cols);
            }
            moreRows = next_combination(rowIdx, rows);
            if (moreRows) {
                std::iota(colIdx.begin(), colIdx.end(), 0);
                moreCols = true;
            }
        }
    }
    if (witness) *witness = "Все миноры нулевые";
    return 0;
}

bool multiply(const Matrix &A, const Matrix &B, Matrix &out) {
    if (A.empty() || B.empty()) return false;
    const size_t aRows = A.size();
    const size_t aCols = A.front().size();
    const size_t bRows = B.size();
    const size_t bCols = B.front().size();
    if (aCols != bRows) return false;
    out.assign(aRows, std::vector<Complex>(bCols, Complex{0.0, 0.0}));
    for (size_t i = 0; i < aRows; ++i) {
        for (size_t k = 0; k < aCols; ++k) {
            const Complex aik = A[i][k];
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
    Matrix aug(n, std::vector<Complex>(2 * n, Complex{0.0, 0.0}));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            aug[i][j] = input[i][j];
        }
        aug[i][n + i] = Complex{1.0, 0.0};
    }
    const double tol = 1e-12;
    for (size_t c = 0; c < n; ++c) {
        size_t pivot = c;
        double maxVal = std::abs(aug[pivot][c]);
        for (size_t i = c + 1; i < n; ++i) {
            double val = std::abs(aug[i][c]);
            if (val > maxVal) {
                maxVal = val;
                pivot = i;
            }
        }
        if (maxVal <= tol) return {false, "Matrix is singular or ill-conditioned"};
        swapRows(aug, pivot, c);
        const Complex pivotVal = aug[c][c];
        for (size_t j = 0; j < 2 * n; ++j) {
            aug[c][j] /= pivotVal;
        }
        for (size_t i = 0; i < n; ++i) {
            if (i == c) continue;
            const Complex factor = aug[i][c];
            if (std::abs(factor) <= tol) continue;
            for (size_t j = 0; j < 2 * n; ++j) {
                aug[i][j] -= factor * aug[c][j];
            }
        }
    }
    inverse.assign(n, std::vector<Complex>(n, Complex{0.0, 0.0}));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            inverse[i][j] = aug[i][n + j];
        }
    }
    return {true, "OK"};
}

MatrixResult determinant(const Matrix &input, Complex &detOut) {
    if (!isSquare(input)) return {false, "Matrix must be square"};
    const size_t n = input.size();
    if (n == 0) return {false, "Matrix is empty"};
    Matrix m = input;
    Complex det = Complex{1.0, 0.0};
    const double tol = 1e-12;
    for (size_t c = 0; c < n; ++c) {
        size_t pivot = c;
        double maxVal = std::abs(m[pivot][c]);
        for (size_t i = c + 1; i < n; ++i) {
            double val = std::abs(m[i][c]);
            if (val > maxVal) {
                maxVal = val;
                pivot = i;
            }
        }
        if (maxVal <= tol) {
            detOut = Complex{0.0, 0.0};
            return {true, "OK"};
        }
        if (pivot != c) {
            swapRows(m, pivot, c);
            det = -det; // row swap flips sign
        }
        const Complex pivotVal = m[c][c];
        det *= pivotVal;
        for (size_t i = c + 1; i < n; ++i) {
            const Complex factor = m[i][c] / pivotVal;
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
    Matrix aug(n, std::vector<Complex>(n + m, Complex{0.0, 0.0}));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) aug[i][j] = A[i][j];
        for (size_t j = 0; j < m; ++j) aug[i][n + j] = rhsCols == 0 ? Complex{0.0, 0.0} : C[i][j];
    }
    const double tol = 1e-12;
    // Forward elimination with partial pivoting, track pivots
    std::vector<size_t> pivotRows;
    std::vector<size_t> pivotCols;
    size_t row = 0;
    for (size_t c = 0; c < n && row < n; ++c) {
        size_t pivot = row;
        double maxVal = std::abs(aug[pivot][c]);
        for (size_t i = row + 1; i < n; ++i) {
            double val = std::abs(aug[i][c]);
            if (val > maxVal) {
                maxVal = val;
                pivot = i;
            }
        }
        if (maxVal <= tol) continue; // no pivot in this column
        swapRows(aug, pivot, row);
        const Complex pivotVal = aug[row][c];
        for (size_t i = row + 1; i < n; ++i) {
            const Complex factor = aug[i][c] / pivotVal;
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
            if (std::abs(aug[r][c]) > tol) { allZero = false; break; }
        }
        if (allZero) {
            for (size_t k = 0; k < m; ++k) {
                if (std::abs(aug[r][n + k]) > tol) {
                    return {false, "System inconsistent (нет решений)"};
                }
            }
        }
    }

    X.assign(n, std::vector<Complex>(m, Complex{0.0, 0.0}));
    // Back substitution using pivot rows; free vars remain zero
    for (int idx = static_cast<int>(rank) - 1; idx >= 0; --idx) {
        const size_t r = pivotRows[static_cast<size_t>(idx)];
        const size_t c = pivotCols[static_cast<size_t>(idx)];
        for (size_t k = 0; k < m; ++k) {
            Complex sum = aug[r][n + k];
            for (size_t j = c + 1; j < n; ++j) {
                sum -= aug[r][j] * X[j][k];
            }
            const Complex denom = aug[r][c];
            if (std::abs(denom) <= tol) return {false, "Matrix is singular or degenerate"};
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

MatrixResult nullspaceVector(const Matrix &A, std::vector<Complex> &v) {
    const auto rowsOpt = nRows(A);
    const auto colsOpt = nCols(A);
    if (!rowsOpt || !colsOpt) return {false, "Matrix is empty"};
    const size_t m = *rowsOpt;
    const size_t n = *colsOpt;
    Matrix aug = A;
    aug.resize(m);
    for (auto &row : aug) row.resize(n + 1, Complex{0.0, 0.0}); // append zero RHS

    const double tol = 1e-12;
    std::vector<size_t> pivotRows;
    std::vector<size_t> pivotCols;
    size_t row = 0;
    for (size_t c = 0; c < n && row < m; ++c) {
        size_t pivot = row;
        double maxVal = std::abs(aug[pivot][c]);
        for (size_t i = row + 1; i < m; ++i) {
            const double val = std::abs(aug[i][c]);
            if (val > maxVal) {
                maxVal = val;
                pivot = i;
            }
        }
        if (maxVal <= tol) continue;
        swapRows(aug, pivot, row);
        const Complex pivotVal = aug[row][c];
        for (size_t i = row + 1; i < m; ++i) {
            const Complex factor = aug[i][c] / pivotVal;
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

    v.assign(n, Complex{0.0, 0.0});
    v[freeCol] = Complex{1.0, 0.0};

    // back substitute
    for (int idx = static_cast<int>(rank) - 1; idx >= 0; --idx) {
        const size_t r = pivotRows[static_cast<size_t>(idx)];
        const size_t c = pivotCols[static_cast<size_t>(idx)];
        Complex sum = Complex{0.0, 0.0};
        for (size_t j = c + 1; j < n; ++j) sum += aug[r][j] * v[j];
        const Complex denom = aug[r][c];
        if (std::abs(denom) <= tol) return {false, "Matrix is singular"};
        v[c] = -sum / denom;
    }
    return {true, "OK"};
}

static bool qrDecomposition(const Matrix &A, Matrix &Q, Matrix &R) {
    const size_t m = A.size();
    if (m == 0) return false;
    const size_t n = A.front().size();
    Q.assign(m, std::vector<Complex>(n, Complex{0.0, 0.0}));
    R.assign(n, std::vector<Complex>(n, Complex{0.0, 0.0}));
    std::vector<std::vector<Complex>> v = A;
    const double tol = 1e-12;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            Complex dot = Complex{0.0, 0.0};
            for (size_t k = 0; k < m; ++k) dot += std::conj(Q[k][j]) * v[k][i];
            R[j][i] = dot;
            for (size_t k = 0; k < m; ++k) v[k][i] -= dot * Q[k][j];
        }
        double norm = 0.0;
        for (size_t k = 0; k < m; ++k) norm += std::norm(v[k][i]);
        norm = std::sqrt(norm);
        if (norm <= tol) return false;
        R[i][i] = Complex{norm, 0.0};
        for (size_t k = 0; k < m; ++k) Q[k][i] = v[k][i] / norm;
    }
    return true;
}

MatrixResult eigenvaluesQR(const Matrix &input, std::vector<Complex> &eigs, int maxIter, double tol) {
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
                offDiag = std::max(offDiag, std::abs(A[i][j]));
            }
        }
        if (offDiag < tol) break;
    }
    eigs.clear();
    eigs.reserve(n);
    for (size_t i = 0; i < n; ++i) eigs.push_back(A[i][i]);
    return {true, "OK"};
}

MatrixResult controllabilityRank(const Matrix &A, const Matrix &B, int &rankOut, std::string *witness) {
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
    rankOut = matrixRank(block, 1e-9, witness);
    return {true, "OK"};
}

MatrixResult observabilityRank(const Matrix &A, const Matrix &C, int &rankOut, std::string *witness) {
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
    rankOut = matrixRank(block, 1e-9, witness);
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
