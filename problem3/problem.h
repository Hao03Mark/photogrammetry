#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <common.hpp>

using namespace photo;
using namespace Eigen;

struct DataIntersection {
    // 左右影像的外方位元素
    Vector3d X_eo[2];
    Vector3d r_eo[2]; // 采用轴-角系统定义
    Matrix3d R_eo[2]; // 采用旋转矩阵定义

    // 同名点
    Vector2d p_img[2];

    // 内方位元素
    const double f = 115.0;
    const double x_0 = 0;
    const double y_0 = 0;
};
const Vector3d X_golden = Vector3d(302.2, 3.5, 8.6);
//这里X_golden是什么？

//翻译：RelativeOrientation——相对定向
struct DataRelativeOrientation {
    // 用于相对定向的控制点
    std::vector<Vector3d> X_ro_golden;

    // 同名点，6对
    std::vector<std::pair<Vector2d, Vector2d>> tiepoints;

    // 内方位元素
    const double f = 115.0;
    const double x_0 = 0;
    const double y_0 = 0;

    // 固定基线长度
    double Bx;
};

//翻译：AbsoluteOrientation——绝对定向
struct DataAbsoluteOrientation {
    // 相对定向元素
    // 同名点，6对
    std::vector<std::pair<Vector2d, Vector2d>> tiepoints;

    // 内方位元素
    const double f = 115.0;
    const double x_0 = 0;
    const double y_0 = 0;

    // 固定基线长度
    double Bx;

    // 对应同名点的控制点坐标
    std::vector<Vector3d> gcps;
};

// 生成模拟数据，大概是行高 2000米处的正直摄影的，重叠度为 60% 左右
DataIntersection generate_intersection();
DataRelativeOrientation generate_relative_orientation();

DataAbsoluteOrientation generate_absolute_orientation();

// TODO: 实现点投影系数法进行前方交会
Vector3d intersection_similarity_transform(const DataIntersection &data);

// TODO: 实现基于共线方程的严密解法
Vector3d intersection_least_squares(const DataIntersection &data, const Vector3d &X_initial);

// TODO: 实现相对定向
std::tuple<Vector3d, Matrix3d> relative_orientation(const DataRelativeOrientation &data);
std::tuple<Vector3d, Matrix3d> relative_orientation_golden(const DataRelativeOrientation &data);

// TODO: 实现绝对定向
DataIntersection absolute_orientation(const DataAbsoluteOrientation &data);
