#pragma once

#include <QMainWindow>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QTextEdit>
#include <QLabel>
#include <QComboBox>
#include <QTabWidget>
#include <QLineEdit>

#include "matrixops.h"

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override = default;

private slots:
    void computeInverse();
    void computeDeterminant();
    void computeEigenvalues();
    void checkControllability();
    void checkObservability();
    void solveSystem();
    void showControllabilityMatrix();
    void showObservabilityMatrix();
    void multiplyMatrices();
    void transposeSelected();
    void inverseSelected();
    void determinantSelected();
    void eigenSelected();
    void evaluateExpression();
    void multiplySelectedSimple();

private:
    QPlainTextEdit *inputA_{nullptr};
    QPlainTextEdit *inputB_{nullptr};
    QPlainTextEdit *inputC_{nullptr};
    QPlainTextEdit *output_{nullptr};

    // Multiply tab
    QPlainTextEdit *inputM1_{nullptr};
    QPlainTextEdit *inputM2_{nullptr};
    QPlainTextEdit *inputM3_{nullptr};
    QPlainTextEdit *outputMul_{nullptr};
    QComboBox *comboUnary_{nullptr};
    QComboBox *comboMulLeft_{nullptr};
    QComboBox *comboMulRight_{nullptr};
    QLineEdit *exprEdit_{nullptr};

    bool fetchMatrixByName(const QString &name, Matrix &out, QString &err) const;

    bool parseMatrix(const QString &text, Matrix &out, QString &err) const;
    QString matrixToString(const Matrix &m) const;
    void appendMessage(const QString &msg);
    void showError(const QString &msg);
};
