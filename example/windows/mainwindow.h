#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QApplication>
#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMessageBox>
#include <QtGui/QPainter>
#include <QtGui/QScreen>
#include <QMainWindow>
#include <QTimer>
#include <QTimer>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <QMouseEvent>
#include <QFileDialog>

#include "../../src/RayTracing/RayTracing.h"
#include "../../src/Window.h"
#include "../../src/Graphics/GraphicsND.h"

#define PI 3.141592653589

class MainWindow : public QMainWindow
{
private:

public:
    MainWindow(QWidget* parent = 0);
    ~MainWindow();

    int len = 1000;
    WindowPix windowPix;

    Mat<> p1, p2, R, G, B, center, direct;
    ObjectLib::ObjectTree obtree;

    QLabel* labelXYZ = new QLabel(this);
    QLabel* labelXY = new QLabel(this);
    QLabel* labelXZ = new QLabel(this);
    QLabel* labelYZ = new QLabel(this);

    void paint();
    void build();
    void showimg();

protected:
    void mousePressEvent(QMouseEvent* event);//鼠标点击事件
    void mouseMoveEvent (QMouseEvent* event);//鼠标点击事件
    void wheelEvent     (QWheelEvent* event);
    void keyPressEvent  (QKeyEvent*   event);
};

#endif // MAINWINDOW_H
