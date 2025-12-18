#pragma once

#include <QMainWindow>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QTextEdit>
#include <QLabel>
#include <QComboBox>

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

private:
    QPlainTextEdit *inputA_{nullptr};
    QPlainTextEdit *inputB_{nullptr};
    QPlainTextEdit *inputC_{nullptr};
    QPlainTextEdit *output_{nullptr};

    bool parseMatrix(const QString &text, Matrix &out, QString &err) const;
    QString matrixToString(const Matrix &m) const;
    void appendMessage(const QString &msg);
    void showError(const QString &msg);
};
