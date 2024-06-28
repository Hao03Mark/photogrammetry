#include "problem.h"
//problem.h用来获取数据的
#include <common.hpp>

#include <iostream>

using namespace photo;

//problem.h中有关于各个函数的说明
//生成模拟数据、实现相对定向，绝对定向、点投影系数法进行前方交会、基于共线方程的严密解法

//DataIntersection是problem.h中定义的一个结构体
//结构体还可以这样使用，蛙趣，长知识了
//返回值类型是结构体
// 生成模拟数据
DataIntersection generate_intersection() {
    DataIntersection data;

    // 行高2000米，重叠度 60% 左右
    data.X_eo[0] = Vector3d(0, 0, 2000);
    data.X_eo[1] = Vector3d(1200, 0, 2000) + Vector3d(3.3, 1.2, -3.9); // 加上一些噪声

    // 约绕 z 轴 旋转 5°
    data.r_eo[0] = Vector3d(0.001, 0.003, 1.0) * 5.1 * DEG2RAD;
    data.r_eo[1] = Vector3d(0.002, 0.001, 1.0) * 4.8 * DEG2RAD;

    //rodrigues是commom.hpp中的方法，将轴-角系统转换为旋转矩阵
    data.R_eo[0] = rodrigues(data.r_eo[0]);
    data.R_eo[1] = rodrigues(data.r_eo[1]);

    //X_golden是什么，没看懂啊
    data.p_img[0] = project(X_golden, data.X_eo[0], data.R_eo[0], data.f, data.x_0, data.y_0);
    data.p_img[1] = project(X_golden, data.X_eo[1], data.R_eo[1], data.f, data.x_0, data.y_0);

    return data;
}

//DataRelativeOrientation也是problem.h中的一个结构体
//获得了内外方位元素，得到了像点坐标、基线长度
// 生成模拟数据
DataRelativeOrientation generate_relative_orientation() {
    std::vector<Vector3d> X_ro_golden(6);
    // 在 -400 和 1600 这两个 X 坐标的地方生成两排点(用于相对定向的控制点)
    X_ro_golden[0] = Vector3d(-400, -1600, 0.0) + Vector3d(3.2, 2.1, 124.4);
    X_ro_golden[1] = Vector3d(-400, 0, 0.0) + Vector3d(2.1, 2.7, 1.4);
    X_ro_golden[2] = Vector3d(-400, 1600, 0.0) + Vector3d(5.3, 2.1, -23.2);

    X_ro_golden[3] = Vector3d(1600, -1600, 0.0) + Vector3d(3.2, 2.1, 65.4);
    X_ro_golden[4] = Vector3d(1600, 0, 0.0) + Vector3d(2.1, 2.7, 321.4);
    X_ro_golden[5] = Vector3d(1600, 1600, 0.0) + Vector3d(5.3, 2.1, 163.2);

    //这里又定义了一个结构体
    DataRelativeOrientation data;
    DataIntersection data_intersection = generate_intersection();

    // 投影到左右相对，得到同名点
    for (int i = 0; i < X_ro_golden.size(); ++i) {
        //这里的project是common中的方法，是得到像点坐标（x，y）
        // lp：left picture左像
        // rp：right picture右像
		Vector2d lp = project(X_ro_golden[i], data_intersection.X_eo[0], data_intersection.R_eo[0], data.f, data.x_0, data.y_0);
		Vector2d rp = project(X_ro_golden[i], data_intersection.X_eo[1], data_intersection.R_eo[1], data.f, data.x_0, data.y_0);
		std::pair<Vector2d, Vector2d> p = std::make_pair(lp, rp);
		data.tiepoints.push_back(p);
        // push_back：在尾部追加一个元素
    }

    data.X_ro_golden = X_ro_golden;

    // 根据世界坐标，计算基线长度
    Vector3d B = data_intersection.X_eo[1] - data_intersection.X_eo[0];
    B = data_intersection.R_eo[0].transpose() * B;
    data.Bx = B.x();
    return data;
}

//DataAbsoluteOrientation也是problem.h中的方法
//得到了Bx,By,Bz
//生成模拟数据
DataAbsoluteOrientation generate_absolute_orientation() {
    DataRelativeOrientation data_ro = generate_relative_orientation();

    DataAbsoluteOrientation data_ao;
    data_ao.gcps = data_ro.X_ro_golden;

    data_ao.tiepoints = data_ro.tiepoints;
    data_ao.Bx = data_ro.Bx;

    return data_ao;
}
