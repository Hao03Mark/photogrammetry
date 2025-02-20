#include <ceres/ceres.h>

int num_observations = 67;

// 2i 为x 2i+1 为y
double data[] = {
    0.000000e+00, 1.133898e+00, 7.500000e-02, 1.334902e+00, 1.500000e-01, 1.213546e+00, 2.250000e-01, 1.252016e+00,
    3.000000e-01, 1.392265e+00, 3.750000e-01, 1.314458e+00, 4.500000e-01, 1.472541e+00, 5.250000e-01, 1.536218e+00,
    6.000000e-01, 1.355679e+00, 6.750000e-01, 1.463566e+00, 7.500000e-01, 1.490201e+00, 8.250000e-01, 1.658699e+00,
    9.000000e-01, 1.067574e+00, 9.750000e-01, 1.464629e+00, 1.050000e+00, 1.402653e+00, 1.125000e+00, 1.713141e+00,
    1.200000e+00, 1.527021e+00, 1.275000e+00, 1.702632e+00, 1.350000e+00, 1.423899e+00, 1.425000e+00, 1.543078e+00,
    1.500000e+00, 1.664015e+00, 1.575000e+00, 1.732484e+00, 1.650000e+00, 1.543296e+00, 1.725000e+00, 1.959523e+00,
    1.800000e+00, 1.685132e+00, 1.875000e+00, 1.951791e+00, 1.950000e+00, 2.095346e+00, 2.025000e+00, 2.361460e+00,
    2.100000e+00, 2.169119e+00, 2.175000e+00, 2.061745e+00, 2.250000e+00, 2.178641e+00, 2.325000e+00, 2.104346e+00,
    2.400000e+00, 2.584470e+00, 2.475000e+00, 1.914158e+00, 2.550000e+00, 2.368375e+00, 2.625000e+00, 2.686125e+00,
    2.700000e+00, 2.712395e+00, 2.775000e+00, 2.499511e+00, 2.850000e+00, 2.558897e+00, 2.925000e+00, 2.309154e+00,
    3.000000e+00, 2.869503e+00, 3.075000e+00, 3.116645e+00, 3.150000e+00, 3.094907e+00, 3.225000e+00, 2.471759e+00,
    3.300000e+00, 3.017131e+00, 3.375000e+00, 3.232381e+00, 3.450000e+00, 2.944596e+00, 3.525000e+00, 3.385343e+00,
    3.600000e+00, 3.199826e+00, 3.675000e+00, 3.423039e+00, 3.750000e+00, 3.621552e+00, 3.825000e+00, 3.559255e+00,
    3.900000e+00, 3.530713e+00, 3.975000e+00, 3.561766e+00, 4.050000e+00, 3.544574e+00, 4.125000e+00, 3.867945e+00,
    4.200000e+00, 4.049776e+00, 4.275000e+00, 3.885601e+00, 4.350000e+00, 4.110505e+00, 4.425000e+00, 4.345320e+00,
    4.500000e+00, 4.161241e+00, 4.575000e+00, 4.363407e+00, 4.650000e+00, 4.161576e+00, 4.725000e+00, 4.619728e+00,
    4.800000e+00, 4.737410e+00, 4.875000e+00, 4.727863e+00, 4.950000e+00, 4.669206e+00
}; //此函数 y=exp(ax+b) 一定会是单调的，我发现数据并不能满足单调，所以数据可能有问题

//定义一个模板化对象来评估残余。每个观测值都会有一个残差。
struct ExponentialResidual {
    ExponentialResidual(double x, double y) : x_(x), y_(y) {}
    //构造代价函数，重载（）符号，仿函数
    template <typename T> bool operator()(const T *const a, const T *const b, T *residual) const {
        residual[0] = y_ - exp(a[0] * x_ + b[0]);
        return true;
    }

private:
    const double x_;
    const double y_;
};

std::tuple<double, double> question2() {
    // TODO: 用上面的 data 数据，拟合函数 y=e^(ax+b)
    //没办法，只好用ceres官方的方法了
    // 最小二乘的初始值使用 0 0
    double a = 0.0, b = 0.0;

    ceres::Problem problem;
    for (int i = 0; i < num_observations; ++i) {
        ceres::CostFunction *cost_function = new ceres::AutoDiffCostFunction<ExponentialResidual, 1, 1, 1>(
            new ExponentialResidual(data[2 * i], data[2 * i + 1]));
        //四个参数分别为（代价函数、核函数、待估参数a、待估参数b）
        problem.AddResidualBlock(cost_function, nullptr, &a, &b);
    }
    //配置Solver
    ceres::Solver::Options options;
    //设置最大迭代次数
    options.max_num_iterations = 66;
    //配置增量方程的解法
    options.linear_solver_type = ceres::DENSE_SCHUR;
    //是否输出到cout（不需要）
    //options.minimizer_progress_to_stdout = true;
    //创建Summary对象用于输出迭代结果
    ceres::Solver::Summary summary;
    //执行求解
    ceres::Solve(options, &problem, &summary);

    return std::make_tuple(a, b);
};

std::tuple<double, double> golden2() {
    double a = 0.0;
    double b = 0.0;

    ceres::Problem problem;
    for (int i = 0; i < num_observations; ++i) {
        ceres::CostFunction* cost = new ceres::AutoDiffCostFunction<ExponentialResidual, 1, 1, 1>(
            new ExponentialResidual(data[2 * i + 0], data[2 * i + 1]));
        problem.AddResidualBlock(cost, nullptr, &a, &b);
    }

    ceres::Solver::Options options;
    options.max_num_iterations = 25;//最大迭代次数
    options.linear_solver_type = ceres::DENSE_QR;
    //options.linear_solver_type = ceres::DENSE_SCHUR;

    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);

    return std::make_tuple(a, b);
}
