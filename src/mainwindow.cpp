#include "mainwindow.h"

#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QMessageBox>
#include <QRegularExpression>
#include <QScrollBar>
#include <QVBoxLayout>
#include <QTabWidget>
#include <algorithm>
#include <complex>
#include <cmath>
#include <string>
#include <sstream>
#include <sstream>

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    auto *tabs = new QTabWidget(this);

    // Tab 1: systems tools
    auto *systemTab = new QWidget;
    {
        auto *layout = new QVBoxLayout;

        auto *grid = new QGridLayout;
        grid->addWidget(new QLabel("A (n x n, коэффициенты)"), 0, 0);
        inputA_ = new QPlainTextEdit;
        inputA_->setPlaceholderText("Строки через перевод строки; элементы через пробел или запятую; поддерживаются комплексные a+bi");
        inputA_->setMinimumHeight(100);
        grid->addWidget(inputA_, 1, 0);

        grid->addWidget(new QLabel("B (n x m, вход)"), 0, 1);
        inputB_ = new QPlainTextEdit;
        inputB_->setPlaceholderText("Для проверки управляемости (по желанию), можно комплексные");
        inputB_->setMinimumHeight(100);
        grid->addWidget(inputB_, 1, 1);

        grid->addWidget(new QLabel("C (p x n / правая часть)"), 0, 2);
        inputC_ = new QPlainTextEdit;
        inputC_->setPlaceholderText("Для наблюдаемости или правой части; если пусто — возьмем нули; можно комплексные");
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

        auto *logRow = new QHBoxLayout;
        auto *clearLogBtn = new QPushButton("Очистить лог");
        logRow->addWidget(clearLogBtn);
        logRow->addStretch();
        layout->addLayout(logRow);

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
        connect(clearLogBtn, &QPushButton::clicked, output_, &QPlainTextEdit::clear);

        systemTab->setLayout(layout);
    }

    // Tab 2: matrix multiply and single-matrix ops
    auto *multiplyTab = new QWidget;
    {
        auto *layout = new QVBoxLayout;
        auto *grid = new QGridLayout;
        grid->addWidget(new QLabel("Матрица M1"), 0, 0);
        inputM1_ = new QPlainTextEdit;
        inputM1_->setPlaceholderText("Строки через перевод строки; элементы через пробел/запятую; комплексные a+bi");
        inputM1_->setMinimumHeight(120);
        grid->addWidget(inputM1_, 1, 0);

        grid->addWidget(new QLabel("Матрица M2"), 0, 1);
        inputM2_ = new QPlainTextEdit;
        inputM2_->setPlaceholderText("Строки через перевод строки; элементы через пробел/запятую; комплексные a+bi");
        inputM2_->setMinimumHeight(120);
        grid->addWidget(inputM2_, 1, 1);

        grid->addWidget(new QLabel("Матрица M3"), 0, 2);
        inputM3_ = new QPlainTextEdit;
        inputM3_->setPlaceholderText("Строки через перевод строки; элементы через пробел/запятую; комплексные a+bi");
        inputM3_->setMinimumHeight(120);
        grid->addWidget(inputM3_, 1, 2);

        layout->addLayout(grid);

        auto *selectRow = new QHBoxLayout;
        comboUnary_ = new QComboBox;
        comboUnary_->addItems({"M1", "M2", "M3"});
        selectRow->addWidget(new QLabel("Операция для матрицы:"));
        selectRow->addWidget(comboUnary_);
        selectRow->addStretch();
        layout->addLayout(selectRow);

        auto *row1 = new QHBoxLayout;
        auto *transposeBtn = new QPushButton("Транспонировать");
        auto *invBtn = new QPushButton("Обратная");
        auto *detBtn2 = new QPushButton("Определитель");
        auto *eigBtn2 = new QPushButton("Собственные значения");
        row1->addWidget(transposeBtn);
        row1->addWidget(invBtn);
        row1->addWidget(detBtn2);
        row1->addWidget(eigBtn2);
        row1->addStretch();
        layout->addLayout(row1);

        auto *rowMul = new QHBoxLayout;
        comboMulLeft_ = new QComboBox;
        comboMulLeft_->addItems({"M1", "M2", "M3"});
        comboMulRight_ = new QComboBox;
        comboMulRight_->addItems({"M1", "M2", "M3"});
        comboMulLeft_->setCurrentIndex(0);
        comboMulRight_->setCurrentIndex(1);
        auto *mulBtnSimple = new QPushButton("Умножить");
        rowMul->addWidget(new QLabel("Левая:"));
        rowMul->addWidget(comboMulLeft_);
        rowMul->addWidget(new QLabel("Правая:"));
        rowMul->addWidget(comboMulRight_);
        rowMul->addWidget(mulBtnSimple);
        rowMul->addStretch();
        layout->addLayout(rowMul);

        auto *exprRow = new QHBoxLayout;
        exprRow->addWidget(new QLabel("Выражение (M1,M2,M3,+,-,*,()):"));
        exprEdit_ = new QLineEdit;
        exprEdit_->setPlaceholderText("например: M1*M2 - M3*M2");
        exprRow->addWidget(exprEdit_, 1);
        auto *evalBtn = new QPushButton("Посчитать выражение");
        exprRow->addWidget(evalBtn);
        layout->addLayout(exprRow);

        auto *logMulRow = new QHBoxLayout;
        auto *clearMulBtn = new QPushButton("Очистить лог");
        logMulRow->addWidget(clearMulBtn);
        logMulRow->addStretch();
        layout->addLayout(logMulRow);

        outputMul_ = new QPlainTextEdit;
        outputMul_->setReadOnly(true);
        outputMul_->setMinimumHeight(200);
        layout->addWidget(outputMul_);

        connect(transposeBtn, &QPushButton::clicked, this, &MainWindow::transposeSelected);
        connect(invBtn, &QPushButton::clicked, this, &MainWindow::inverseSelected);
        connect(detBtn2, &QPushButton::clicked, this, &MainWindow::determinantSelected);
        connect(eigBtn2, &QPushButton::clicked, this, &MainWindow::eigenSelected);
        connect(evalBtn, &QPushButton::clicked, this, &MainWindow::evaluateExpression);
        connect(mulBtnSimple, &QPushButton::clicked, this, &MainWindow::multiplySelectedSimple);
        connect(clearMulBtn, &QPushButton::clicked, outputMul_, &QPlainTextEdit::clear);

        multiplyTab->setLayout(layout);
    }

    tabs->addTab(systemTab, "Системы");
    tabs->addTab(multiplyTab, "Произведение");

    setCentralWidget(tabs);
    setWindowTitle("Матрицы — управление и наблюдаемость");
    resize(1040, 720);
}

static bool parseRealToken(const QString &token, double &val) {
    QString t = token.trimmed();
    bool ok = false;
    val = t.toDouble(&ok);
    if (ok) return true;
    const int slashPos = t.indexOf('/');
    if (slashPos > 0 && slashPos < t.size() - 1) {
        bool okNum = false, okDen = false;
        double num = t.left(slashPos).toDouble(&okNum);
        double den = t.mid(slashPos + 1).toDouble(&okDen);
        if (okNum && okDen && std::fabs(den) > 1e-15) {
            val = num / den;
            return true;
        }
    }
    return false;
}

static bool parseComplexToken(const QString &token, Complex &val) {
    QString t = token.trimmed();
    if (t.startsWith('(') && t.endsWith(')')) {
        t = t.mid(1, t.size() - 2).trimmed();
    }
    t = t.replace('I', 'i');
    // No 'i' -> real
    if (!t.contains('i', Qt::CaseInsensitive)) {
        double real = 0.0;
        if (!parseRealToken(t, real)) return false;
        val = Complex{real, 0.0};
        return true;
    }

    // Handle forms: a+bi, a-bi, +bi, -bi, bi, -i, i
    QString withoutI = t;
    if (withoutI.endsWith('i', Qt::CaseInsensitive)) {
        withoutI.chop(1);
    }
    int splitPos = -1;
    for (int i = withoutI.size() - 1; i > 0; --i) {
        const QChar ch = withoutI[i];
        if (ch == '+' || ch == '-') { splitPos = i; break; }
    }

    double re = 0.0;
    double im = 0.0;
    if (splitPos == -1) {
        // pure imaginary
        QString imagStr = withoutI.trimmed();
        if (imagStr.isEmpty() || imagStr == "+") imagStr = "1";
        if (imagStr == "-") imagStr = "-1";
        if (!parseRealToken(imagStr, im)) return false;
    } else {
        QString realStr = withoutI.left(splitPos).trimmed();
        QString imagStr = withoutI.mid(splitPos).trimmed();
        if (!parseRealToken(realStr, re)) return false;
        if (imagStr.isEmpty() || imagStr == "+") imagStr = "1";
        if (imagStr == "-") imagStr = "-1";
        if (!parseRealToken(imagStr, im)) return false;
    }
    val = Complex{re, im};
    return true;
}

bool MainWindow::parseMatrix(const QString &text, Matrix &out, QString &err) const {
    const QString cleaned = text.trimmed();
    if (cleaned.isEmpty()) {
        err = "Input is empty";
        return false;
    }
    const auto rows = cleaned.split('\n', Qt::SkipEmptyParts);
    Matrix data;
    int expectedCols = -1;
    bool ok = true;
    for (const auto &rowStr : rows) {
        QString rowNormalized = rowStr;
        rowNormalized.replace(',', ' ');
        rowNormalized.replace('|', ' ');
        const auto parts = rowNormalized.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
        if (parts.isEmpty()) continue;
        std::vector<Complex> row;
        for (const auto &p : parts) {
            Complex val{0.0, 0.0};
            if (!parseComplexToken(p, val)) {
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
    auto fmt = [](const Complex &z) {
        std::ostringstream os;
        os.setf(std::ios::fixed);
        os.precision(6);
        const double r = z.real();
        const double im = z.imag();
        const double tol = 1e-9;
        if (std::abs(im) <= tol) {
            os << r;
        } else if (std::abs(r) <= tol) {
            os << im << 'i';
        } else {
            os << r << (im >= 0 ? "+" : "") << im << 'i';
        }
        return os.str();
    };

    std::ostringstream oss;
    for (size_t i = 0; i < m.size(); ++i) {
        for (size_t j = 0; j < m[i].size(); ++j) {
            oss << fmt(m[i][j]);
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
    std::vector<Complex> eigs;
    auto res = eigenvaluesQR(A, eigs);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    std::ostringstream oss;
    for (size_t i = 0; i < eigs.size(); ++i) {
        oss << matrixToString({{eigs[i]}}).toStdString();
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
    Complex det = Complex{0.0, 0.0};
    auto res = determinant(A, det);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    appendMessage(QString("Определитель A: %1").arg(matrixToString({{det}})));
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
            Matrix Atrim(n, std::vector<Complex>(n, Complex{0.0, 0.0}));
            Matrix Cauto(n, std::vector<Complex>(rhsCols, Complex{0.0, 0.0}));
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
        C.assign(n, std::vector<Complex>(1, Complex{0.0, 0.0}));
    }

    Matrix X;
    std::vector<size_t> pivotCols;
    auto res = solveLinearSystem(A, C, X, &pivotCols);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    // Build human-readable mapping x1, x2, ...
    auto fmt = [this](const Complex &z) { return matrixToString({{z}}).toStdString(); };
    std::ostringstream oss;
    const size_t nVars = X.size();
    const size_t rhsCols = X.empty() ? 0 : X.front().size();
    for (size_t k = 0; k < rhsCols; ++k) {
        if (rhsCols > 1) oss << "Правая часть " << (k + 1) << ": ";
        for (size_t i = 0; i < nVars; ++i) {
            oss << "x" << (i + 1) << " = " << fmt(X[i][k]);
            if (i + 1 < nVars) oss << ", ";
        }
        if (k + 1 < rhsCols) oss << "\n";
    }
    QString status = QString::fromStdString(res.message);
    if (status.isEmpty()) status = "OK";
    QString extra;
    if (status.startsWith("Бесконечное")) {
        // Try to build one ненулевое решение: X_particular + v_null
        std::vector<Complex> vnull;
        auto nres = nullspaceVector(coeffForNull, vnull);
        if (nres.ok && !vnull.empty()) {
            // Form example: add null vector to first RHS solution (k=0)
            if (!X.empty() && !X.front().empty()) {
                std::ostringstream ex;
                ex << "Пример ненулевого решения: ";
                for (size_t i = 0; i < X.size(); ++i) {
                    Complex val = X[i][0] + vnull[i];
                    ex << "x" << (i + 1) << " = " << matrixToString({{val}}).toStdString();
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
    std::string witness;
    auto res = controllabilityRank(A, B, rank, &witness);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    const int n = static_cast<int>(A.size());
    QString extra = witness.empty() ? QString() : QString::fromStdString(witness);
    appendMessage(QString("Ранг управляемости = %1 (n = %2) -> %3")
                      .arg(rank)
                      .arg(n)
                      .arg(rank == n ? "УПРАВЛЯЕМА" : "НЕ управляемая")
                      + (extra.isEmpty() ? QString() : QString("\n") + extra));
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
    std::string witness;
    auto res = observabilityRank(A, C, rank, &witness);
    if (!res.ok) {
        showError(QString::fromStdString(res.message));
        return;
    }
    const int n = static_cast<int>(A.size());
    QString extra = witness.empty() ? QString() : QString::fromStdString(witness);
    appendMessage(QString("Ранг наблюдаемости = %1 (n = %2) -> %3")
                      .arg(rank)
                      .arg(n)
                      .arg(rank == n ? "НАБЛЮДАЕМА" : "НЕ наблюдаемая")
                      + (extra.isEmpty() ? QString() : QString("\n") + extra));
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

void MainWindow::multiplyMatrices() {
    Matrix M1, M2, P;
    QString err;
    if (!parseMatrix(inputM1_->toPlainText(), M1, err)) {
        if (!err.isEmpty()) showError(err);
        return;
    }
    if (!parseMatrix(inputM2_->toPlainText(), M2, err)) {
        if (!err.isEmpty()) showError(err);
        return;
    }
    if (!multiply(M1, M2, P)) {
        if (nCols(M1) && nRows(M2)) {
            outputMul_->appendPlainText(QString("Ошибка: размеры не совпадают (%1x%2)·(%3x%4)")
                                            .arg(*nRows(M1))
                                            .arg(*nCols(M1))
                                            .arg(*nRows(M2))
                                            .arg(*nCols(M2)));
        } else {
            outputMul_->appendPlainText("Ошибка умножения (проверьте ввод)");
        }
        auto *bar = outputMul_->verticalScrollBar();
        bar->setValue(bar->maximum());
        return;
    }
    outputMul_->appendPlainText("M1·M2:\n" + matrixToString(P));
    auto *bar = outputMul_->verticalScrollBar();
    bar->setValue(bar->maximum());
}

// Helpers for multiply tab selections
bool MainWindow::fetchMatrixByName(const QString &name, Matrix &out, QString &err) const {
    QString src;
    if (name == "M1") src = inputM1_->toPlainText();
    else if (name == "M2") src = inputM2_->toPlainText();
    else if (name == "M3") src = inputM3_->toPlainText();
    else {
        err = "Неизвестная матрица";
        return false;
    }
    return parseMatrix(src, out, err);
}

void appendMul(QPlainTextEdit *out, const QString &msg) {
    out->appendPlainText(msg);
    auto *bar = out->verticalScrollBar();
    bar->setValue(bar->maximum());
}

static bool addOrSub(const Matrix &A, const Matrix &B, Matrix &out, bool add) {
    if (A.empty() || B.empty()) return false;
    if (A.size() != B.size() || A.front().size() != B.front().size()) return false;
    out.assign(A.size(), std::vector<Complex>(A.front().size(), Complex{0.0, 0.0}));
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            out[i][j] = A[i][j] + (add ? B[i][j] : -B[i][j]);
        }
    }
    return true;
}

void MainWindow::transposeSelected() {
    Matrix M;
    QString err;
    const QString name = comboUnary_->currentText();
    if (!fetchMatrixByName(name, M, err)) {
        appendMul(outputMul_, "Ошибка: " + err);
        return;
    }
    Matrix T(M.front().size(), std::vector<Complex>(M.size(), Complex{0.0, 0.0}));
    for (size_t i = 0; i < M.size(); ++i) {
        for (size_t j = 0; j < M[i].size(); ++j) T[j][i] = M[i][j];
    }
    appendMul(outputMul_, QString("T(%1):\n%2").arg(name).arg(matrixToString(T)));
}

void MainWindow::inverseSelected() {
    Matrix M, Inv;
    QString err;
    const QString name = comboUnary_->currentText();
    if (!fetchMatrixByName(name, M, err)) {
        appendMul(outputMul_, "Ошибка: " + err);
        return;
    }
    auto res = invertMatrix(M, Inv);
    if (!res.ok) {
        appendMul(outputMul_, "Ошибка: " + QString::fromStdString(res.message));
        return;
    }
    appendMul(outputMul_, QString("Обратная %1:\n%2").arg(name).arg(matrixToString(Inv)));
}

void MainWindow::determinantSelected() {
    Matrix M;
    QString err;
    const QString name = comboUnary_->currentText();
    if (!fetchMatrixByName(name, M, err)) {
        appendMul(outputMul_, "Ошибка: " + err);
        return;
    }
    Complex detVal = Complex{0.0, 0.0};
    auto res = determinant(M, detVal);
    if (!res.ok) {
        appendMul(outputMul_, "Ошибка: " + QString::fromStdString(res.message));
        return;
    }
    appendMul(outputMul_, QString("det(%1) = %2").arg(name).arg(matrixToString({{detVal}})));
}

void MainWindow::eigenSelected() {
    Matrix M;
    QString err;
    const QString name = comboUnary_->currentText();
    if (!fetchMatrixByName(name, M, err)) {
        appendMul(outputMul_, "Ошибка: " + err);
        return;
    }
    std::vector<Complex> eigs;
    auto res = eigenvaluesQR(M, eigs);
    if (!res.ok) {
        appendMul(outputMul_, "Ошибка: " + QString::fromStdString(res.message));
        return;
    }
    std::ostringstream oss;
    for (size_t i = 0; i < eigs.size(); ++i) {
        oss << matrixToString({{eigs[i]}}).toStdString();
        if (i + 1 < eigs.size()) oss << ", ";
    }
    appendMul(outputMul_, QString("λ(%1): %2").arg(name).arg(QString::fromStdString(oss.str())));
}

void MainWindow::multiplySelectedSimple() {
    Matrix L, R, P;
    QString err;
    const QString leftName = comboMulLeft_ ? comboMulLeft_->currentText() : "M1";
    const QString rightName = comboMulRight_ ? comboMulRight_->currentText() : "M2";
    if (!fetchMatrixByName(leftName, L, err)) {
        appendMul(outputMul_, "Ошибка: " + err);
        return;
    }
    if (!fetchMatrixByName(rightName, R, err)) {
        appendMul(outputMul_, "Ошибка: " + err);
        return;
    }
    if (!multiply(L, R, P)) {
        if (nCols(L) && nRows(R)) {
            appendMul(outputMul_, QString("Ошибка: размеры не совпадают (%1x%2)·(%3x%4)")
                                       .arg(*nRows(L))
                                       .arg(*nCols(L))
                                       .arg(*nRows(R))
                                       .arg(*nCols(R)));
        } else {
            appendMul(outputMul_, "Ошибка умножения (проверьте ввод)");
        }
        return;
    }
    appendMul(outputMul_, QString("%1·%2:\n%3")
                               .arg(leftName)
                               .arg(rightName)
                               .arg(matrixToString(P)));
}

// multiplySelected removed

void MainWindow::evaluateExpression() {
    const QString expr = exprEdit_->text().trimmed();
    if (expr.isEmpty()) {
        appendMul(outputMul_, "Введите выражение");
        return;
    }

    // Tokenize
    struct Tok { enum Type { Mat, Plus, Minus, Mul, LPar, RPar } t; QString name; };
    QList<Tok> tokens;
    int i = 0;
    while (i < expr.size()) {
        const QChar ch = expr[i];
        if (ch.isSpace()) { ++i; continue; }
        if (ch == '+') { tokens.push_back({Tok::Plus, {}}); ++i; continue; }
        if (ch == '-') { tokens.push_back({Tok::Minus, {}}); ++i; continue; }
        if (ch == '*') { tokens.push_back({Tok::Mul, {}}); ++i; continue; }
        if (ch == '(') { tokens.push_back({Tok::LPar, {}}); ++i; continue; }
        if (ch == ')') { tokens.push_back({Tok::RPar, {}}); ++i; continue; }
        if (ch.toUpper() == 'M') {
            if (i + 1 >= expr.size()) { appendMul(outputMul_, "Ошибка парсинга матрицы"); return; }
            QChar num = expr[i + 1];
            if (num != '1' && num != '2' && num != '3') { appendMul(outputMul_, "Ожидал M1/M2/M3"); return; }
            QString name = QString("M") + num;
            tokens.push_back({Tok::Mat, name});
            i += 2;
            continue;
        }
        appendMul(outputMul_, QString("Неизвестный символ: %1").arg(ch));
        return;
    }

    // Shunting-yard to RPN
    QList<Tok> output;
    QList<Tok> stack;
    auto prec = [](Tok::Type t) {
        if (t == Tok::Mul) return 2;
        if (t == Tok::Plus || t == Tok::Minus) return 1;
        return 0;
    };
    for (const auto &tk : tokens) {
        if (tk.t == Tok::Mat) {
            output.push_back(tk);
        } else if (tk.t == Tok::Plus || tk.t == Tok::Minus || tk.t == Tok::Mul) {
            while (!stack.isEmpty()) {
                Tok top = stack.back();
                if ((top.t == Tok::Plus || top.t == Tok::Minus || top.t == Tok::Mul) && prec(top.t) >= prec(tk.t)) {
                    output.push_back(top);
                    stack.pop_back();
                } else break;
            }
            stack.push_back(tk);
        } else if (tk.t == Tok::LPar) {
            stack.push_back(tk);
        } else if (tk.t == Tok::RPar) {
            bool matched = false;
            while (!stack.isEmpty()) {
                Tok top = stack.back(); stack.pop_back();
                if (top.t == Tok::LPar) { matched = true; break; }
                output.push_back(top);
            }
            if (!matched) { appendMul(outputMul_, "Скобки не сбалансированы"); return; }
        }
    }
    while (!stack.isEmpty()) {
        Tok top = stack.back(); stack.pop_back();
        if (top.t == Tok::LPar || top.t == Tok::RPar) { appendMul(outputMul_, "Скобки не сбалансированы"); return; }
        output.push_back(top);
    }

    // Eval RPN
    QList<Matrix> st;
    for (const auto &tk : output) {
        if (tk.t == Tok::Mat) {
            Matrix M; QString err;
            if (!fetchMatrixByName(tk.name, M, err)) { appendMul(outputMul_, "Ошибка: " + err); return; }
            st.push_back(std::move(M));
        } else if (tk.t == Tok::Plus || tk.t == Tok::Minus || tk.t == Tok::Mul) {
            if (st.size() < 2) { appendMul(outputMul_, "Ошибка: недостаточно операндов"); return; }
            Matrix b = st.back(); st.pop_back();
            Matrix a = st.back(); st.pop_back();
            Matrix r;
            bool ok = false;
            if (tk.t == Tok::Plus) ok = addOrSub(a, b, r, true);
            else if (tk.t == Tok::Minus) ok = addOrSub(a, b, r, false);
            else if (tk.t == Tok::Mul) ok = multiply(a, b, r);
            if (!ok) {
                appendMul(outputMul_, "Ошибка: несовместимые размеры");
                return;
            }
            st.push_back(std::move(r));
        }
    }
    if (st.size() != 1) { appendMul(outputMul_, "Ошибка: неверное выражение"); return; }
    appendMul(outputMul_, "Результат выражения:\n" + matrixToString(st.back()));
}
