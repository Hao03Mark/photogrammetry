#include "problem.h"

#include <Eigen/Geometry>
#include <ceres/ceres.h>
#include <ceres/rotation.h>

using namespace ceres;

DataIntersection absolute_orientation(const DataAbsoluteOrientation &data_ao) {

    DataIntersection result;

    // 1. 利用同名点数据，生成相对定向结果
    Vector3d B;
    Matrix3d Rr;
    DataRelativeOrientation data_ro;
    data_ro.Bx = data_ao.Bx;
    data_ro.tiepoints = data_ao.tiepoints;
    std::tie(B, Rr) = relative_orientation(data_ro);// 调用p2，生成相对定向结果B，Rr
    
    // 2. 根据相对定向结果和同名点，计算模型坐标，即在第一个影像坐标系下进行前方交会
    std::vector<Vector3d> model_points(data_ao.tiepoints.size());
    Vector3d temp;
    double X1, Y1, Z1, X2, Y2, Z2;
    double N1, N2;
    for (int i = 0; i < model_points.size(); i++) {
        // X1,Y1,Z1
        X1 = data_ro.tiepoints[i].first[0];
        Y1 = data_ro.tiepoints[i].first[1];
        Z1 = -data_ro.f;
		// X1,Y1,Z1
        temp[0] = data_ro.tiepoints[i].second[0];
        temp[1] = data_ro.tiepoints[i].second[1];
        temp[2] = -data_ro.f;
        temp = Rr * temp;
        X2 = temp[0];
        Y2 = temp[1];
        Z2 = temp[2];
       
        // N1,N2
        N1 = (B.x() * Z2 - B.z() * X2) / (X1 * Z2 - X2 * Z1);
        N2 = (B.x() * Z1 - B.z() * X1) / (X1 * Z2 - X2 * Z1);
      
        // 模型点坐标
        model_points[i][0] =  N1 * X1;
        model_points[i][1] =  0.5 * (N1 * Y1 + N2 * Y2 + B.y());
        model_points[i][2] =  N1 * Z1;
    } // 计算模型坐标完成

    
    // 3. 根据前方交会的模型坐标和控制点，计算旋转，缩放，平移参数
    double s = 1.0;
    Vector3d ra = Vector3d::Zero();
    Vector3d T = Vector3d::Zero();

    // 3.1 计算尺度、缩放、平移的初始值

    //先用最垃圾的方法算一把，把报告交了
    MatrixXd A(3 * model_points.size(), 7);
    VectorXd l(3 * model_points.size());
    Eigen::Matrix<double,7,1> delta;
    Vector3d Temp_l(1);
    Matrix3d K;// 临时旋转矩阵
    bool tem;// 判断限差
   

    // 非重心化
    do {
        K = rodrigues(ra);
        for (int i = 0; i < model_points.size(); i++) {

            Temp_l = data_ao.gcps[i] - s * K * model_points[i] - T;
            l[3 * i + 0] = Temp_l[0];
            l[3 * i + 1] = Temp_l[1];
            l[3 * i + 2] = Temp_l[2];
            //
            A(3 * i + 0, 0) = 1;
            A(3 * i + 0, 1) = 0;
            A(3 * i + 0, 2) = 0;
            A(3 * i + 0, 3) = model_points[i][0];
            A(3 * i + 0, 4) = -model_points[i][2];
            A(3 * i + 0, 5) = 0;
            A(3 * i + 0, 6) = -model_points[i][1];
            //
            A(3 * i + 1, 0) = 0;
            A(3 * i + 1, 1) = 1;
            A(3 * i + 1, 2) = 0;
            A(3 * i + 1, 3) = model_points[i][1];
            A(3 * i + 1, 4) = 0;
            A(3 * i + 1, 5) = -model_points[i][2];
            A(3 * i + 1, 6) = model_points[i][0];
            //
            A(3 * i + 2, 0) = 0;
            A(3 * i + 2, 1) = 0;
            A(3 * i + 2, 2) = 1;
            A(3 * i + 2, 3) = model_points[i][2];
            A(3 * i + 2, 4) = model_points[i][0];
            A(3 * i + 2, 5) = model_points[i][1];
            A(3 * i + 2, 6) = 0;
        }

        delta = (A.transpose() * A).inverse() * A.transpose() * l;
        // 重心化坐标后，平移量无需改正
        T[0] += delta[0];
        T[1] += delta[1];
        T[2] += delta[2];
        // s = s * (1 + delta[3]);
        s += delta[3];
        ra[0] += delta[4];
        ra[1] += delta[5];
        ra[2] += delta[6];

        tem = fabs(delta[4]) < 3.0e-5 && fabs(delta[5]) < 3.0e-5 && fabs(delta[6]) < 3.0e-5;

    } while (!tem);

    // 3.2 构造最小二乘优化，利用控制点和模型坐标，计算7参数
    // ceres::Problem problem;
    // std::vector<ceres::ResidualBlockId> blocks;
    // ceres solver方法还不太会，就先不用了

    // 4. 根据 7 参数信息，生成最后的外方位元素
    result.R_eo[0] = rodrigues(ra);
    result.X_eo[0] = T;
    result.R_eo[1] = rodrigues(ra) * Rr;
    result.X_eo[1] = T + s * rodrigues(ra) * B;

    return result;
}


