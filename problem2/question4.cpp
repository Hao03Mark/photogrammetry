#include <Eigen/Core>
#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include <vector>

#include <autodiff/reverse/reverse.hpp>
#include <autodiff/reverse/eigen.hpp> //头文件的先后顺序竟然也会报错


using namespace Eigen;

// 相平面坐标
std::vector<Vector2d> point2ds = {Vector2d(-86.15, -68.99), Vector2d(-53.40, 82.21), Vector2d(-14.78, -76.63),
                                  Vector2d(10.46, 64.43)};
// 控制点坐标
std::vector<Vector3d> point3ds = {Vector3d(36589.41, 25273.32, 2195.17), Vector3d(37631.08, 31324.51, 728.69),
                                  Vector3d(39100.97, 24934.98, 2386.50), Vector3d(40426.54, 30319.81, 757.31)};
// 内方位元素
double f = 153.24; //这里f单位是毫米，控制点坐标是米,像平面坐标也是毫米
double x_0 = 0.0;  //我觉得老师写的是有问题的，他并没有化成统一的单位（比如：米），我化为同一单位还会出错，可以证明胡老师是没有统一为米单位的
double y_0 = 0.0;

//胡老师没给比例尺1/m；

using namespace autodiff;

Matrix3var rodrigues(const Vector3var &angle_axis) {
    static const double kOne = 1.0;
    const var theta2 = angle_axis.squaredNorm();
    Matrix3var R;
    if (theta2 > var(std::numeric_limits<double>::epsilon())) {
        // We want to be careful to only evaluate the square root if the
        // norm of the angle_axis vector is greater than zero. Otherwise
        // we get a division by zero.
        const var theta = sqrt(theta2);
        const var wx = angle_axis[0] / theta;
        const var wy = angle_axis[1] / theta;
        const var wz = angle_axis[2] / theta;

        const var costheta = cos(theta);
        const var sintheta = sin(theta);

        R(0, 0) = costheta + wx * wx * (kOne - costheta);
        R(1, 0) = wz * sintheta + wx * wy * (kOne - costheta);
        R(2, 0) = -wy * sintheta + wx * wz * (kOne - costheta);
        R(0, 1) = wx * wy * (kOne - costheta) - wz * sintheta;
        R(1, 1) = costheta + wy * wy * (kOne - costheta);
        R(2, 1) = wx * sintheta + wy * wz * (kOne - costheta);
        R(0, 2) = wy * sintheta + wx * wz * (kOne - costheta);
        R(1, 2) = -wx * sintheta + wy * wz * (kOne - costheta);
        R(2, 2) = costheta + wz * wz * (kOne - costheta);
    } else {
        // Near zero, we switch to using the first order Taylor expansion.
        R(0, 0) = kOne;
        R(1, 0) = angle_axis[2];
        R(2, 0) = -angle_axis[1];
        R(0, 1) = -angle_axis[2];
        R(1, 1) = kOne;
        R(2, 1) = angle_axis[0];
        R(0, 2) = angle_axis[1];
        R(1, 2) = -angle_axis[0];
        R(2, 2) = kOne;
    }
    return R;
}

std::tuple<Vector3d, Matrix3d, double, Vector3d> question4() {
    Vector3d Xs = Vector3d::Zero();    // 表示Xs，Ys，Zs
    Matrix3d R = Matrix3d::Identity(); // 3*3的单位矩阵
    double sigma_0 = 0.0;              //验后单位权中误差
    Vector3d m_Xs = Vector3d::Zero();  //线元素的理论精度

    // // 终于搞懂了，重新写(看来这种高级的写法还是不行，还是使用low的方法先算出来吧)
    // // L矩阵
    // std::vector<Vector2d> L(4);
    // // V矩阵
    // std::vector<Vector2d> V(4);
    // // A矩阵
    // std::vector<Matrix<double,2,6,RowMajor> > A(4);
    // // 改正数向量
    // std::vector<double> delta(6);

    // low方法（我感觉很呆很傻，不像vector<vector> point,这样其实包含了点的属性，很好）
    //在使用Visual Studio 19 调试的时候发现，不初始化会有内存溢出报错的情况，所以只好进行初始化
    //我又想到了一种解决方案，在for循环结束后使用resize方法可以重置矩阵，应该可以，还没有尝试
    // L 8*1
    VectorXd L(2 * point3ds.size());
    // V 8*1
    VectorXd V(2 * point3ds.size());
    // A 8*6
    MatrixXd A(2 * point3ds.size(), 6);
    //防止列优先可以使用这样std::Matrix<double,Dynamic,Dynamic,Rowcal>
    //改正数向量 6*1
    VectorXd delta(6);
    //权逆阵
    MatrixXd Qx(6,6);

    // TODO 实现后方交会

    // 1 计算外方位线元素初值
    Xs[0] = (point3ds[0](0) + point3ds[1](0) + point3ds[2](0) + point3ds[3](0)) / 4;
    Xs[1] = (point3ds[0](1) + point3ds[1](1) + point3ds[2](1) + point3ds[3](1)) / 4;
    // Xs[2] = m * f;
    Xs[2] = 8000.0; //根据胡老师下面的数据给出的

    // 2 角元素初值就为单位阵
    // 我的思考：角元素初始值为0啊，怎么会是单位阵呢？旋转矩阵初值才为单位阵啊。
    double phi = 0.0, omega = 0.0, kappa = 0.0;
    // 3 迭代优化求解
    do {

        // 旋转矩阵R
        R(0, 0) = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
        R(0, 1) = -1 * cos(phi) * sin(kappa) - 1 * sin(phi) * sin(omega) * cos(kappa);
        R(0, 2) = -1 * sin(phi) * cos(omega);
        R(1, 0) = cos(omega) * sin(kappa);
        R(1, 1) = cos(omega) * cos(kappa);
        R(1, 2) = -1 * sin(omega);
        R(2, 0) = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
        R(2, 1) = -1 * sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
        R(2, 2) = cos(phi) * cos(omega);

        for (int i = 0; i < point3ds.size(); i++) {

            // 便于计算，Xbar，Ybar，Zbar
            double Xbar = R(0, 0) * (point3ds[i](0) - Xs[0]) + R(1, 0) * (point3ds[i](1) - Xs[1]) +
                          R(2, 0) * (point3ds[i](2) - Xs[2]),
                   Ybar = R(0, 1) * (point3ds[i](0) - Xs[0]) + R(1, 1) * (point3ds[i](1) - Xs[1]) +
                          R(2, 1) * (point3ds[i](2) - Xs[2]),
                   Zbar = R(0, 2) * (point3ds[i](0) - Xs[0]) + R(1, 2) * (point3ds[i](1) - Xs[1]) +
                          R(2, 2) * (point3ds[i](2) - Xs[2]);
            // 近似值（x），（y）
            double x_around = -1 * f * Xbar / Zbar, 
                y_around = -1 * f * Ybar / Zbar;

            // L矩阵赋值
            L[2 * i] = point2ds[i](0) - x_around;
            L[2 * i + 1] = point2ds[i](1) - y_around;

            // A矩阵第一行a11——a16
            A(2 * i, 0) = (R(0, 0) * f + R(0, 2) * (point2ds[i](0) - x_0)) / Zbar; // a11
            A(2 * i, 1) = (R(1, 0) * f + R(1, 2) * (point2ds[i](0) - x_0)) / Zbar; // a12
            A(2 * i, 2) = (R(2, 0) * f + R(2, 2) * (point2ds[i](0) - x_0)) / Zbar; // a13
            A(2 * i, 3) = (point2ds[i](1) - y_0) * sin(omega) -
                          ((point2ds[i](0) - x_0) / f *
                               ((point2ds[i](0) - x_0) * cos(kappa) - (point2ds[i](1) - y_0) * sin(kappa)) +
                           f * cos(kappa)) *
                              cos(omega); // a14
            A(2 * i, 4) = -f * sin(kappa) -
                          (point2ds[i](0) - x_0) / f *
                              ((point2ds[i](0) - x_0) * sin(kappa) + (point2ds[i](1) - y_0) * cos(kappa)); // a15
            A(2 * i, 5) = (point2ds[i](1) - y_0);                                                          // a16

            // A矩阵第二行a21——a26
            A(2 * i + 1, 0) = (R(0, 1) * f + R(0, 2) * (point2ds[i](1) - y_0)) / Zbar; // a21
            A(2 * i + 1, 1) = (R(1, 1) * f + R(1, 2) * (point2ds[i](1) - y_0)) / Zbar; // a22
            A(2 * i + 1, 2) = (R(2, 1) * f + R(2, 2) * (point2ds[i](1) - y_0)) / Zbar; // a23
            A(2 * i + 1, 3) = -(point2ds[i](0) - x_0) * sin(omega) -
                              ((point2ds[i](1) - y_0) / f *
                                   ((point2ds[i](0) - x_0) * cos(kappa) - (point2ds[i](1) - y_0) * sin(kappa)) -
                               f * sin(kappa)) *
                                  cos(omega); // a24
            A(2 * i + 1, 4) = -f * cos(kappa) -
                              (point2ds[i](1) - y_0) / f *
                                  ((point2ds[i](0) - x_0) * sin(kappa) + (point2ds[i](1) - y_0) * cos(kappa)); // a25
            A(2 * i + 1, 5) = -(point2ds[i](0) - x_0);                                                         // a26                                        
        }

        // 改正数向量
        delta = (A.transpose() * A).inverse() * A.transpose() * L;
        Xs[0] += delta[0];
        Xs[1] += delta[1];
        Xs[2] += delta[2];
        phi += delta[3];
        omega += delta[4];
        kappa += delta[5];

    } while (fabs(delta[3]) > 3.0e-5 || fabs(delta[4]) > 3.0e-5 || fabs(delta[5]) > 3.0e-5);

    // 迭代结束，计算误差方程
    V = A * delta - L;

    // 单位权中误差
    sigma_0 = sqrt((V.transpose() * V / (2 * point3ds.size() - 6))(0, 0));
    // 理论精度还未写

    Qx = (A.transpose() * A).inverse();

    //理论精度mi
    m_Xs[0] = sigma_0 * sqrt(Qx(0, 0));
    m_Xs[1] = sigma_0 * sqrt(Qx(1, 1));
    m_Xs[2] = sigma_0 * sqrt(Qx(2, 2));

    return std::make_tuple(Xs, R, sigma_0, m_Xs);
}

struct ProjectionResidual {
    ProjectionResidual(const Vector2d &point2d, const Vector3d &point3d) : point2d(point2d), point3d(point3d) {}

    template <typename T> bool operator()(const T *Xs, const T *r, T *residual) const {
        T R[9];
        ceres::AngleAxisToRotationMatrix(r, ceres::RowMajorAdapter3x3(R));

        T X = point3d.x() - Xs[0];
        T Y = point3d.y() - Xs[1];
        T Z = point3d.z() - Xs[2];

        T x = -f * (R[0] * X + R[3] * Y + R[6] * Z) / (R[2] * X + R[5] * Y + R[8] * Z) + x_0;
        T y = -f * (R[1] * X + R[4] * Y + R[7] * Z) / (R[2] * X + R[5] * Y + R[8] * Z) + y_0;

        residual[0] = point2d.x() - x;
        residual[1] = point2d.y() - y;

        return true;
    }

private:
    Vector2d point2d;
    Vector3d point3d;
};

std::tuple<Vector3d, Matrix3d, double, Vector3d> golden4() {
    Vector3d Xs = Vector3d::Zero();
    Matrix3d R = Matrix3d::Identity();
    double sigma_0 = 0.0;
    Vector3d m_Xs = Vector3d::Zero();

    // 计算初值
    for (const auto &point : point3ds) {
        Xs += point;
    }
    Xs = Xs / point3ds.size();
    Xs.z() = 8000.0;

    // 角元素初值就是单位阵
    auto angle_axis = Eigen::AngleAxisd(R);
    Vector3d r = angle_axis.angle() * angle_axis.axis();

    ceres::Problem problem;

    std::vector<ceres::ResidualBlockId> blocks;

    for (int i = 0; i < point3ds.size(); ++i) {
        ceres::CostFunction *cost = new ceres::AutoDiffCostFunction<ProjectionResidual, 2, 3, 3>(
            new ProjectionResidual(point2ds[i], point3ds[i]));
        ceres::ResidualBlockId block = problem.AddResidualBlock(cost, nullptr, Xs.data(), r.data());
        blocks.push_back(block);
    }

    ceres::Solver::Options op;
    op.max_num_iterations = 100;
    op.linear_solver_type = ceres::DENSE_QR;

    ceres::Solver::Summary summary;
    ceres::Solve(op, &problem, &summary);

    // 计算单位权中误差
    {
        ceres::Problem::EvaluateOptions op_eval;
        op_eval.residual_blocks = blocks;
        double cost;
        problem.Evaluate(op_eval, &cost, nullptr, nullptr, nullptr);
        cost *= 2.0; // ceres 记录的是包含 1/2 的
        sigma_0 = std::sqrt(cost / (2 * point2ds.size() - 6));
    }

    // 计算协方差矩阵
    {
        Matrix3d cov_Xs;
        std::vector<std::pair<const double *, const double *>> covariance_blocks;
        covariance_blocks.emplace_back(Xs.data(), Xs.data());
        ceres::Covariance::Options options;
        ceres::Covariance covariance(options);
        covariance.Compute(covariance_blocks, &problem);
        covariance.GetCovarianceBlock(Xs.data(), Xs.data(), cov_Xs.data());

        m_Xs.x() = sigma_0 * std::sqrt(cov_Xs(0, 0));
        m_Xs.y() = sigma_0 * std::sqrt(cov_Xs(1, 1));
        m_Xs.z() = sigma_0 * std::sqrt(cov_Xs(2, 2));
    }

    // 计算旋转矩阵
    {
        angle_axis.angle() = r.norm();
        angle_axis.axis() = r.normalized();

        R = angle_axis.toRotationMatrix();
    }
    return std::make_tuple(Xs, R, sigma_0, m_Xs);
}
