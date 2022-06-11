#include "mainwindow.h"


MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent)
{
    srand(time(NULL));

    QSize screenSize = QSize(1600, 1000);
    this->setMinimumSize(screenSize);
    this->setMaximumSize(screenSize);

    windowPix.init(len, len / 3);

    p1.zero(3);
    p2.zero(3);
    center.zero(3);
    direct.zero(3); direct = { 0,200, 0 };
    R.zero(len, len);
    G.zero(len, len);
    B.zero(len, len);

    labelXYZ->setGeometry(0, 0, len, len);
    labelXY->setGeometry(len + 1, 0, len / 3, len / 3);
    labelXZ->setGeometry(len + 1, len / 3, len / 3, len / 3);
    labelYZ->setGeometry(len + 1, 2 * len / 3, len / 3, len / 3);
    
    showimg();
}

MainWindow::~MainWindow()
{}

void MainWindow::build() {
    Mat<> p1(3), p2(3), R(1000, 1000), G(1000, 1000), B(1000, 1000), M(3, 3);

    Material* m = new Material;
    m->color = 1; m->refract = 1; m->refractRate = 1.8;

    //obtree.addStl("C:/Users/29753/Desktop/cube.stl", p1 = { 0, 0, 0 }, 100, &m);

    m = new Material;
    m->color = { 0.2, 1, 0.2 }; m->quickReflect = 1;
    //obtree.addStl("C:/Users/29753/Desktop/cube2.stl", p1 = { 0, 0, 0 }, 100, &m);
    m = new Material;
    m->color = 1; m->quickReflect = 1;
    obtree.addEllipsoid(p1 = { 0, 0, 0 }, M = { 1.0/1000, 0, 0,  0, 1.0 / 1000, 0,  0, 0, 1.0 / 1000 }, m);
    //obtree.addSphere(p1 = { 0, 0, 0 }, 1000, m);

    m = new Material;
    m->color = { 1, 1, 1 };  m->quickReflect = 1;
    obtree.addPlane(p1 = { 0, 0, 1 }, p2 = { 0, 0, -100 }, m);
    /*
    m = new Material;
    m->color = { 1, 0, 1 };  m->rediate = 1;
    obtree.addPlane(p1 = { 0, 0, -1 }, p2 = { 0, 0, 150 }, m);*/

    obtree.build();

    RayTracing::PointLight.push_back(p1 = { 500,500,-500 });
}

void MainWindow::paint() {

    RayTracing::traceRay(
        center, direct, 0.2,
        obtree,
        R, G, B, 0, 1
    );

    for (int x = 0; x < len; x++) {
        for (int y = 0; y < len; y++) {
            windowPix.imgXYZ(x, y) = 0xFF000000
                + std::min((int)(R(x, y) * 0xFF), 0xFF) * 0x10000
                + std::min((int)(G(x, y) * 0xFF), 0xFF) * 0x100
                + std::min((int)(B(x, y) * 0xFF), 0xFF);
        }
    }

    showimg();

}

void MainWindow::showimg() {
    QImage imgQXYZ((uchar*)windowPix.imgXYZ.data, len, len, QImage::Format_ARGB32);

    labelXYZ->setPixmap(QPixmap::fromImage(imgQXYZ));
    labelXYZ->show();
}

void MainWindow::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) {    //左键按下
        int x = event->x(),
            y = event->y();

        if (x <= len && y <= len) {
            build();
            paint();
        }

    }
    else if (event->button() == Qt::RightButton) {
        int x = event->x(),
            y = event->y();

        if (x <= len && y <= len) {


        }
        else {  // Save Image

            QImage imgQXYZ((uchar*)windowPix.imgXYZ.data, len, len, QImage::Format_ARGB32);

            QString filename = QFileDialog::getSaveFileName(this,
                tr("Save Image"),
                "",
                tr("*.bmp;; *.png;; *.jpg;; *.tif;; *.GIF")
            );

            if (filename.isEmpty())
            {
                return;
            }
            else
            {
                if (!(imgQXYZ.save(filename))) //保存图像
                {

                    QMessageBox::information(this,
                        tr("Failed to save the image"),
                        tr("Failed to save the image!")
                    );
                    return;
                }

                QMessageBox::information(this, tr("OK!"), tr("OK!"));
            }
        }
    }
    else if (event->button() == Qt::MiddleButton) {

    }

}

void MainWindow::mouseMoveEvent(QMouseEvent* event) {
    if (event->buttons() & Qt::LeftButton) {
        if (event->x() <= len && event->y() <= len) {
        

        }
    }
    
}

void MainWindow::wheelEvent(QWheelEvent* event)    // 滚轮事件
{
    if (event->angleDelta().y() > 0) {                    // 当滚轮远离使用者时
        double dis = norm(direct);
        mul(direct, dis + 50, normalize(direct));
    }
    else {                                     // 当滚轮向使用者方向旋转时
        double dis = norm(direct);
        mul(direct, dis - 50, normalize(direct));
    }

    paint();
}

void MainWindow::keyPressEvent(QKeyEvent* event)
{
    double l = norm(direct);
    Mat<> dir = direct, tmp(3);
    normalize(dir);

    // 普通键
    switch (event->key())
    {
    case Qt::Key_A: {
        normalize(cross(dir, dir, tmp = { 0,0,1 }));
        add(center, center, mul(dir, 50, dir));
        break;
    }
    case Qt::Key_D: {
        normalize(cross(dir, dir, tmp = { 0,0,1 }));
        sub(center, center, mul(dir, 50, dir));
        break;
    }
    case Qt::Key_W: {
        add(center, center, mul(dir, 50, dir));
        break;
    }
    case Qt::Key_S: {
        sub(center, center, mul(dir, 50, dir));
        break;
    }
    case Qt::Key_Q: center(2) += 50; break;
    case Qt::Key_E: center(2) -= 50; break;

    case Qt::Key_R: 
        normalize(direct);
        direct(0) = direct(0) * cos(10 * 2 * PI / 360) + direct(1) * sin(10 * 2 * PI / 360);
        direct(1) = direct(0) *-sin(10 * 2 * PI / 360) + direct(1) * cos(10 * 2 * PI / 360);

        normalize(direct);
        mul(direct, l, direct);
        break;
    case Qt::Key_F:
        normalize(direct);
        direct(0) = direct(0) * cos(-10 * 2 * PI / 360) + direct(1) * sin(-10 * 2 * PI / 360);
        direct(1) = direct(0) *-sin(-10 * 2 * PI / 360) + direct(1) * cos(-10 * 2 * PI / 360);

        normalize(direct);
        mul(direct, l, direct);
        break;
    case Qt::Key_T:
        normalize(direct);
        direct(0) = direct(0) * cos(10 * 2 * PI / 360) + direct(2) * sin(10 * 2 * PI / 360);
        direct(2) = direct(0) *-sin(10 * 2 * PI / 360) + direct(2) * cos(10 * 2 * PI / 360);

        normalize(direct);
        mul(direct, l, direct);
        break;
    case Qt::Key_G:
        normalize(direct);
        direct(0) = direct(0) * cos(-10 * 2 * PI / 360) + direct(2) * sin(-10 * 2 * PI / 360);
        direct(2) = direct(0) *-sin(-10 * 2 * PI / 360) + direct(2) * cos(-10 * 2 * PI / 360);

        normalize(direct);
        mul(direct, l, direct);
        break;
    case Qt::Key_Y: 
        normalize(direct);
        direct(1) = direct(1) * cos(10 * 2 * PI / 360) + direct(2) * sin(10 * 2 * PI / 360);
        direct(2) = direct(1) *-sin(10 * 2 * PI / 360) + direct(2) * cos(10 * 2 * PI / 360);

        normalize(direct);
        mul(direct, l, direct);
        break;
    case Qt::Key_H: 
        normalize(direct);
        direct(1) = direct(1) * cos(-10 * 2 * PI / 360) + direct(2) * sin(-10 * 2 * PI / 360);
        direct(2) = direct(1)  -sin(-10 * 2 * PI / 360) + direct(2) * cos(-10 * 2 * PI / 360);

        normalize(direct);
        mul(direct, l, direct);
        break;

    case Qt::Key_Space:
        R.zero();
        G.zero();
        B.zero();
        
        RayTracing::traceRay(
            center, direct, 0.2,
            obtree,
            R, G, B, 0, 10
        );

        for (int x = 0; x < len; x++) {
            for (int y = 0; y < len; y++) {
                windowPix.imgXYZ(x, y) = 0xFF000000
                    + std::min((int)(R(x, y) * 0xFF), 0xFF) * 0x10000
                    + std::min((int)(G(x, y) * 0xFF), 0xFF) * 0x100
                    + std::min((int)(B(x, y) * 0xFF), 0xFF);
            }
        }
        showimg(); return;
        break;
    }

    paint();
    return;
}