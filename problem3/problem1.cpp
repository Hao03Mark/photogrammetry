#include "problem.h"
#include <iostream>

#include <Eigen/Core>
#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include <cmath>
#include <vector>
using namespace Eigen;
using namespace ceres;

double f = 115.0;
double x_0 = 0;
double y_0 = 0;

//使用点投影系数法
Vector3d intersection_similarity_transform(const DataIntersection &data) {
    // 1.计算基线初值
    Eigen::Vector3d B;
    B = data.X_eo[1] - data.X_eo[0];
    B[1] = 0;
    B[2] = 0;

    // 2.计算像点坐标的像空间辅助坐标系坐标
    Eigen::Vector3d X1;
    X1[0] = data.p_img[0](0);
    X1[1] = data.p_img[0](1);
    X1[2] = -f;
    X1 = data.R_eo[0] * X1;

    Eigen::Vector3d X2;
    X2[0] = data.p_img[1](0);
    X2[1] = data.p_img[1](1);
    X2[2] = -f;
    X2 = data.R_eo[1] * X2;

    // 3.计算点投影系数
    double N1, N2;
    N1 = (B[0] * X2[2] - B[2] * X2[0]) / (X1[0] * X2[2] - X2[0] * X1[2]);
    N2 = (B[0] * X1[2] - B[2] * X1[0]) / (X1[0] * X2[2] - X2[0] * X1[2]);

    // 4.计算地面坐标
    Eigen::Vector3d XA;
    XA[0] = data.X_eo[0](0) + N1 * X1[0];
    XA[1] = (data.X_eo[0](1) + N1 * X1[1] + data.X_eo[1](1) + N2 * X2[1]) / 2;
    XA[2] = data.X_eo[0](2) + N1 * X1[2];

    return XA;
}

//使用严密法，这里的 x_initial 是 上面的XA
Vector3d intersection_least_squares(const DataIntersection &data, const Vector3d &X_initial) {
    // 已知值x0,y0,f,Xs1,Ys1,Zs1,phi1,omega1,kappa1,Xs2,Ys2,Zs2,phi2,omega2,kappa2
    // 观测值x1,y1,x2,y2
    // 未知数X,Y,Z
    Vector3d Xi = X_initial; // initial——最初的；初始的
    Vector3d delta;          // 改正量
    MatrixXd A(4, 3);        // 两个点，四行三列
    VectorXd l(4);           // 常数项
    Vector2d x[2];           // 像点坐标
    const double H = 2000;   // 行高2000
    bool tem;

    do {
        // 量测像点坐标，project方法即可（初始值，外方位线元素，旋转矩阵，f，x0，y0）
        x[0] = project(Xi, data.X_eo[0], data.R_eo[0], data.f, data.x_0, data.y_0);
        x[1] = project(Xi, data.X_eo[1], data.R_eo[1], data.f, data.x_0, data.y_0);

        for (int i = 0; i < 2; i++) {
            A(2 * i, 0) = f * cos(data.r_eo[i](2)) / H;
            A(2 * i, 1) = f * sin(data.r_eo[i](2)) / H;
            A(2 * i, 2) = x[i](0) / H;
            l[2 * i] = data.p_img[i](0) - x[i](0);

            A(2 * i + 1, 0) = -f * sin(data.r_eo[i](2)) / H;
            A(2 * i + 1, 1) = f * cos(data.r_eo[i](2)) / H;
            A(2 * i + 1, 2) = x[i](1) / H;
            l[2 * i + 1] = data.p_img[i](1) - x[i](1);
        }

        delta = (A.transpose() * A).inverse() * A.transpose() * l;
        Xi += delta;

        tem = fabs(delta[0]) < 1e-3 && fabs(delta[1]) < 1e-3 && fabs(delta[2]) < 1e-3;

    } while (!tem);

    return Xi;
}