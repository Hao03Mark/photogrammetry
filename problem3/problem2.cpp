#include "problem.h"
#include <iostream>

#include <Eigen/Core>
#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include <cmath>
#include <vector>
using namespace Eigen;

std::tuple<Vector3d, Matrix3d> relative_orientation(const DataRelativeOrientation &data) {
    // 1. 给定初值
    // x0，y0，f，像点坐标，用于相对定向的控制点坐标，基线长度都来自data
    Vector3d B(data.Bx, 0, 0);     // Bx确定为连续法相对定向
    Vector3d r = Vector3d::Zero(); // 采用 rodrigues 转角
    Vector2d uv(0, 0);
    double f = data.f; // 赋初值完成

    // K为临时旋转矩阵，T为临时系数矩阵A的角色
    Matrix3d K;
    MatrixXd T(6, 5);
    VectorXd l(6);
    // 2. 采用矩阵行列式作为观测值求解
    VectorXd delta(5);
    // 为了计算方便，所以引入下列名称
    double Bx, By, Bz, N1, N2;
    // 相对定向的控制点像空间辅助坐标
    Vector3d Xl, Xr;
    bool tem;
    do {
        K = rodrigues(r);
        Bx = B[0];
        By = B[1];
        Bz = B[2];
        for (int i = 0; i < data.tiepoints.size(); i++) {

            
            Xl[0] = data.tiepoints[i].first[0];// X1
            Xr[0] = data.tiepoints[i].second[0];// X2
            Xl[1] = data.tiepoints[i].first[1];// Y1
            Xr[1] = data.tiepoints[i].second[1];// Y2
            Xl[2] = -f;// Z1
            Xr[2] = -f;// Z2
            Xr = K * Xr;// 乘以（相对）旋转矩阵得像空间辅助

            N1 = (Bx * Xr[2] - Bz * Xr[0]) / (Xl[0] * Xr[2] - Xr[0] * Xl[2]);
            N2 = (Bx * Xl[2] - Bz * Xl[0]) / (Xl[0] * Xr[2] - Xl[2] * Xr[0]);

            T(i, 0) = Bx;
            T(i, 1) = -Xr[1] * Bx / Xr[2];
            T(i, 2) = -Xr[0] * Xr[1] / Xr[2] * N2;
            T(i, 3) = -(Xr[2] + Xr[1] * Xr[1] / Xr[2]) * N2;
            T(i, 4) = Xr[0] * N2;

            l[i] = N1 * Xl[1] - (N2 * Xr[1] + By);// l即为Q
        }

        delta = (T.transpose() * T).inverse() * T.transpose() * l;
        uv[0] += delta[0];
        uv[1] += delta[1];
        r[0] = r[0] + delta[2];
        r[1] = r[1] + delta[3];
        r[2] = r[2] + delta[4];
        // 相对定向元素初值
        B[1] = B[0] * uv[0];
        B[2] = B[0] * uv[1];

        tem = fabs(delta[0]) < 3e-5 && fabs(delta[1]) < 3e-5 && fabs(delta[2]) < 3e-5 && fabs(delta[3]) < 3e-5 &&
              fabs(delta[4]) < 3e-5;

    } while (!tem);

    // 3. 返回解
    Matrix3d R = rodrigues(r);
    B = Vector3d(B.x(), B.x() * uv(0), B.x() * uv(1));

    return std::tuple<Vector3d, Matrix3d>(B, R);
}

std::tuple<Vector3d, Matrix3d> relative_orientation_golden(const DataRelativeOrientation &_) {
    DataIntersection data = generate_intersection();
    Vector3d B = data.R_eo[0].transpose() * (data.X_eo[1] - data.X_eo[0]);
    Matrix3d R = data.R_eo[0].transpose() * data.R_eo[1];

    return std::make_tuple(B, R);
}
