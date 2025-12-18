#include "mainwindow.h"

#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QMessageBox>
#include <QRegularExpression>
#include <QScrollBar>
#include <QVBoxLayout>
#include <algorithm>
#include <sstream>

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    auto *central = new QWidget(this);
    auto *layout = new QVBoxLayout;

    auto *grid = new QGridLayout;
    grid->addWidget(new QLabel("A (n x n, коэффициенты)"), 0, 0);
    inputA_ = new QPlainTextEdit;
    inputA_->setPlaceholderText("Строки через перевод строки; элементы через пробел или запятую");
    inputA_->setMinimumHeight(100);
    grid->addWidget(inputA_, 1, 0);

    grid->addWidget(new QLabel("B (n x m, вход)"), 0, 1);
    inputB_ = new QPlainTextEdit;
    inputB_->setPlaceholderText("Для проверки управляемости (по желанию)");
    inputB_->setMinimumHeight(100);
    grid->addWidget(inputB_, 1, 1);

    grid->addWidget(new QLabel("C (p x n / правая часть)"), 0, 2);
    inputC_ = new QPlainTextEdit;
    inputC_->setPlaceholderText("Для наблюдаемости или правой части; если пусто — возьмем нули");
    inputC_->setMinimumHeight(100);
    grid->addWidget(inputC_, 1, 2);

    layout->addLayout(grid);

    auto *buttonRow = new QHBoxLayout;
    auto *inverseBtn = new QPushButton("Обратная A");
    auto *detBtn = new QPushButton("Определитель A");
    auto *solveBtn = new QPushButton("Решить A·x = C");
    auto *eigBtn = new QPushButton("Собственные значения A");
    auto *ctrlBtn = new QPushButton("Управляемость?");
    auto *obsBtn = new QPushButton("Наблюдаемость?");
    auto *ctrlMatBtn = new QPushButton("Матрица управляемости");
    auto *obsMatBtn = new QPushButton("Матрица наблюдаемости");

    buttonRow->addWidget(inverseBtn);
    buttonRow->addWidget(detBtn);
    buttonRow->addWidget(solveBtn);
    buttonRow->addWidget(eigBtn);
    buttonRow->addWidget(ctrlBtn);
    buttonRow->addWidget(obsBtn);
    buttonRow->addWidget(ctrlMatBtn);
    buttonRow->addWidget(obsMatBtn);
    buttonRow->addStretch();
    layout->addLayout(buttonRow);

    output_ = new QPlainTextEdit;
    output_->setReadOnly(true);
    output_->setMinimumHeight(160);
    layout->addWidget(output_);

    connect(inverseBtn, &QPushButton::clicked, this, &MainWindow::computeInverse);
    connect(detBtn, &QPushButton::clicked, this, &MainWindow::computeDeterminant);
    connect(solveBtn, &QPushButton::clicked, this, &MainWindow::solveSystem);
    connect(eigBtn, &QPushButton::clicked, this, &MainWindow::computeEigenvalues);
    connect(ctrlBtn, &QPushButton::clicked, this, &MainWindow::checkControllability);
    connect(obsBtn, &QPushButton::clicked, this, &MainWindow::checkObservability);
    connect(ctrlMatBtn, &QPushButton::clicked, this, &MainWindow::showControllabilityMatrix);
    connect(obsMatBtn, &QPushButton::clicked, this, &MainWindow::showObservabilityMatrix);

    central->setLayout(layout);
    setCentralWidget(central);
    setWindowTitle("Матрицы — управление и наблюдаемость");
    resize(960, 640);
}

bool MainWindow::parseMatrix(const QString &text, Matrix &out, QString &err) const {
    const QString cleaned = text.trimmed();
    if (cleaned.isEmpty()) {
        err = "Input is empty";
        return false;
    }
    const auto rows = cleaned.split('\n', Qt::SkipEmptyParts);
    std::vector<std::vector<double>> data;
    int expectedCols = -1;
    bool ok = true;
    for (const auto &rowStr : rows) {
        QString rowNormalized = rowStr;
        rowNormalized.replace(',', ' ');
        rowNormalized.replace('|', ' ');
        const auto parts = rowNormalized.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
        if (parts.isEmpty()) continue;
        std::vector<double> row;
        for (const auto &p : parts) {
            const double val = p.toDouble(&ok);
            if (!ok) {
                err = QString("Cannot parse '%1'").arg(p);
                return false;
            }
            row.push_back(val);
        }
        if (expectedCols == -1) expectedCols = static_cast<int>(row.size());
        if (static_cast<int>(row.size()) != expectedCols) {
            err = "Row lengths differ";
            return false;
        }
        data.push_back(std::move(row));
    }
    if (data.empty()) {
        err = "No data";
        return false;
    }
    out = std::move(data);
    return true;
}

QString MainWindow::matrixToString(const Matrix &m) const {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(6);
    for (size_t i = 0; i < m.size(); ++i) {
        for (size_t j = 0; j < m[i].size(); ++j) {
            oss << m[i][j];
            if (j + 1 < m[i].size()) oss << ' ';
        }
        if (i + 1 < m.size()) oss << '\n';
    }
    return QString::fromStdString(oss.str());
}

void MainWindow::appendMessage(const QString &msg) {
    output_->appendPlainText(msg);
    auto *bar = output_->verticalScrollBar();
    bar->setValue(bar->maximum());
}

void MainWindow::showError(const QString &msg) {
    appendMessage("Ошибка: " + msg);
}

void MainWindow::computeInverse() {
    Matrix A, inv;
    QString err;
    if (!parseMatrix(inputA_->toPlainText(), A, err)) {
        showError(err);
        return;
    }
    auto res = invertMatrix(A, inv);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    appendMessage("Обратная матрица A:\n" + matrixToString(inv));
}

void MainWindow::computeEigenvalues() {
    Matrix A;
    QString err;
    if (!parseMatrix(inputA_->toPlainText(), A, err)) {
        showError(err);
        return;
    }
    std::vector<double> eigs;
    auto res = eigenvaluesQR(A, eigs);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(6);
    for (size_t i = 0; i < eigs.size(); ++i) {
        oss << eigs[i];
        if (i + 1 < eigs.size()) oss << ", ";
    }
    appendMessage(QString("Собственные значения: %1").arg(QString::fromStdString(oss.str())));
}

void MainWindow::computeDeterminant() {
    Matrix A;
    QString err;
    if (!parseMatrix(inputA_->toPlainText(), A, err)) {
        showError(err);
        return;
    }
    double det = 0.0;
    auto res = determinant(A, det);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(6);
    oss << det;
    appendMessage(QString("Определитель A: %1").arg(QString::fromStdString(oss.str())));
}

void MainWindow::solveSystem() {
    Matrix A, C;
    QString err;
    if (!parseMatrix(inputA_->toPlainText(), A, err)) {
        showError(err);
        return;
    }
    Matrix coeffForNull = A; // keep a copy of coefficient part for nullspace computation
    const QString cText = inputC_->toPlainText().trimmed();
    bool rhsProvided = false;
    if (!cText.isEmpty()) {
        if (!parseMatrix(inputC_->toPlainText(), C, err)) {
            showError("C: " + err);
            return;
        }
        rhsProvided = true;
    }

    // Auto-handle augmented input: if A has extra columns and C is empty, treat trailing columns as RHS
    if (!rhsProvided) {
        const auto rowsOpt = nRows(A);
        const auto colsOpt = nCols(A);
        if (rowsOpt && colsOpt && *colsOpt > *rowsOpt) {
            const size_t n = *rowsOpt;
            const size_t totalCols = *colsOpt;
            const size_t rhsCols = totalCols - n;
            Matrix Atrim(n, std::vector<double>(n, 0.0));
            Matrix Cauto(n, std::vector<double>(rhsCols, 0.0));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) Atrim[i][j] = A[i][j];
                for (size_t j = 0; j < rhsCols; ++j) Cauto[i][j] = A[i][n + j];
            }
            A.swap(Atrim);
            C.swap(Cauto);
            rhsProvided = true;
            appendMessage("Авто-разбор A|C: последние столбцы приняты как правая часть");
        }
    }

    if (!rhsProvided) {
        const size_t n = A.size();
        C.assign(n, std::vector<double>(1, 0.0));
    }

    Matrix X;
    std::vector<size_t> pivotCols;
    auto res = solveLinearSystem(A, C, X, &pivotCols);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    // Build human-readable mapping x1, x2, ...
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(6);
    const size_t nVars = X.size();
    const size_t rhsCols = X.empty() ? 0 : X.front().size();
    for (size_t k = 0; k < rhsCols; ++k) {
        if (rhsCols > 1) oss << "Правая часть " << (k + 1) << ": ";
        for (size_t i = 0; i < nVars; ++i) {
            oss << "x" << (i + 1) << " = " << X[i][k];
            if (i + 1 < nVars) oss << ", ";
        }
        if (k + 1 < rhsCols) oss << "\n";
    }
    QString status = QString::fromStdString(res.message);
    if (status.isEmpty()) status = "OK";
    QString extra;
    if (status.startsWith("Бесконечное")) {
        // Try to build one ненулевое решение: X_particular + v_null
        std::vector<double> vnull;
        auto nres = nullspaceVector(coeffForNull, vnull);
        if (nres.ok && !vnull.empty()) {
            // Form example: add null vector to first RHS solution (k=0)
            if (!X.empty() && !X.front().empty()) {
                std::ostringstream ex;
                ex.setf(std::ios::fixed);
                ex.precision(6);
                ex << "Пример ненулевого решения: ";
                for (size_t i = 0; i < X.size(); ++i) {
                    double val = X[i][0] + vnull[i];
                    ex << "x" << (i + 1) << " = " << val;
                    if (i + 1 < X.size()) ex << ", ";
                }
                extra = QString::fromStdString(ex.str());
            }
        }
    }

    QString result = QString("Решение A·x = C (%1):\n%2")
                          .arg(status)
                          .arg(QString::fromStdString(oss.str()));
    if (!extra.isEmpty()) result += "\n" + extra;
    appendMessage(result);
}

void MainWindow::checkControllability() {
    Matrix A, B;
    QString err;
    if (!parseMatrix(inputA_->toPlainText(), A, err)) {
        showError(err);
        return;
    }
    if (!parseMatrix(inputB_->toPlainText(), B, err)) {
        showError("B: " + err);
        return;
    }
    int rank = 0;
    auto res = controllabilityRank(A, B, rank);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    const int n = static_cast<int>(A.size());
    appendMessage(QString("Ранг управляемости = %1 (n = %2) -> %3")
                      .arg(rank)
                      .arg(n)
                      .arg(rank == n ? "УПРАВЛЯЕМА" : "НЕ управляемая"));
}

void MainWindow::checkObservability() {
    Matrix A, C;
    QString err;
    if (!parseMatrix(inputA_->toPlainText(), A, err)) {
        showError(err);
        return;
    }
    if (!parseMatrix(inputC_->toPlainText(), C, err)) {
        showError("C: " + err);
        return;
    }
    int rank = 0;
    auto res = observabilityRank(A, C, rank);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    const int n = static_cast<int>(A.size());
    appendMessage(QString("Ранг наблюдаемости = %1 (n = %2) -> %3")
                      .arg(rank)
                      .arg(n)
                      .arg(rank == n ? "НАБЛЮДАЕМА" : "НЕ наблюдаемая"));
}

void MainWindow::showControllabilityMatrix() {
    Matrix A, B, block;
    QString err;
    if (!parseMatrix(inputA_->toPlainText(), A, err)) {
        showError(err);
        return;
    }
    if (!parseMatrix(inputB_->toPlainText(), B, err)) {
        showError("B: " + err);
        return;
    }
    auto res = controllabilityMatrix(A, B, block);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    appendMessage("Матрица управляемости:\n" + matrixToString(block));
}

void MainWindow::showObservabilityMatrix() {
    Matrix A, C, block;
    QString err;
    if (!parseMatrix(inputA_->toPlainText(), A, err)) {
        showError(err);
        return;
    }
    if (!parseMatrix(inputC_->toPlainText(), C, err)) {
        showError("C: " + err);
        return;
    }
    auto res = observabilityMatrix(A, C, block);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    appendMessage("Матрица наблюдаемости:\n" + matrixToString(block));
}
