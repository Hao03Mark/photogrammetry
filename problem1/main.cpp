//#define CATCH_CONFIG_MAIN
//#include <catch2/catch.hpp> //老版本的用法
#include <catch2/catch_test_macros.hpp> //catch2 V3.3.2的用法
#include <catch2/catch_approx.hpp>      //catch2中approx的用法
#include <catch_amalgamated.hpp>
// 新版本要加的东西真的好多

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <unsupported/Eigen/EulerAngles>

#include <iostream>

using namespace Eigen;
using namespace std;
// DEG角度制，RAD弧度制
constexpr double DEG2RAD = 0.01745329251994329576923690768489; // 1度对应的弧度
constexpr double RAD2DEG = 57.295779513082320876798154814105;  // 1弧度对应角度

//计算绕 x 轴旋转 30∘的三维旋转矩阵
Matrix3d question1() {
    double angle_x = 30.0 * DEG2RAD;
    Matrix3d R;

    // TODO: 填充旋转矩阵的９个元素
    R(0, 0) = 1.0;
    R(0, 1) = 0.0;
    R(0, 2) = 0.0;
    R(1, 0) = 0.0;
    R(1, 1) = cos(angle_x);
    R(1, 2) = -sin(angle_x);
    R(2, 0) = 0.0;
    R(2, 1) = sin(angle_x);
    R(2, 2) = cos(angle_x);

    return R;
}

//计算在 ϕ−ω−κ 转角系统下 ϕ=5∘,ω=5∘,κ=89∘的旋转矩阵
Matrix3d question2() {
    double phi, omega, kappa;
    phi = 5.0 * DEG2RAD;
    omega = 5.0 * DEG2RAD;
    kappa = 89.0 * DEG2RAD;

    Matrix3d R;

    // TODO: 填充旋转矩阵的９个元素
    R(0, 0) = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
    R(0, 1) = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
    R(0, 2) = -sin(phi) * cos(omega);
    R(1, 0) = cos(omega) * sin(kappa);
    R(1, 1) = cos(omega) * cos(kappa);
    R(1, 2) = -sin(omega);
    R(2, 0) = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
    R(2, 1) = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
    R(2, 2) = cos(phi) * cos(omega);

    return R;
}

//计算上述旋转矩阵在 ω′−ϕ′−κ′ 转角下的旋转角
Vector3d question3() {
    Matrix3d R = question2();

    // omega'-phi'-kappa'
    // TODO: 填充３个角元素
    Vector3d r;
    r(0) = -atan2(R(1, 2), R(2, 2)) * RAD2DEG;
    r(1) = -asin(R(0, 2)) * RAD2DEG;
    r(2) = -atan2(R(0, 1), R(0, 0)) * RAD2DEG;

    return r;
}

Matrix3d golden1() {
    double angle_x = 30.0 * DEG2RAD;
    Matrix3d R = AngleAxisd(angle_x, Vector3d::UnitX()).toRotationMatrix();
    return R;
}

Matrix3d golden2() {
    double phi, omega, kappa;
    phi = 5.0 * DEG2RAD;
    omega = 5.0 * DEG2RAD;
    kappa = 89.0 * DEG2RAD;

    // 教材中的 phi 旋转其实不符合右手法则
    Matrix3d R_phi = AngleAxisd(-phi, Vector3d::UnitY()).toRotationMatrix();
    Matrix3d R_omega = AngleAxisd(omega, Vector3d::UnitX()).toRotationMatrix();
    Matrix3d R_kappa = AngleAxisd(kappa, Vector3d::UnitZ()).toRotationMatrix();

    Matrix3d R = R_phi * R_omega * R_kappa;
    return R;
}

Vector3d golden3() {
    Matrix3d R = golden2();

    using OPKSystem = EulerSystem<-EULER_X, EULER_Y, EULER_Z>;
    using OPK = EulerAngles<double, OPKSystem>;

    OPK opk = OPK(R);
    Vector3d r = opk.angles() * RAD2DEG;
    r(0) = -r(0);
    r(1) = -r(1);
    return r;
}

TEST_CASE("Question 1", "The rotation angle around x axis with 30 degrees") {
    Matrix3d R_student = question1();
    Matrix3d R_golden = golden1();

    // 新版本的Approx会报错。改为Catch::Approx就好了
    REQUIRE(R_student.determinant() == Catch::Approx(1.0));       // Approx表示约等于
    REQUIRE((R_student - R_golden).norm() == Catch::Approx(0.0)); // norm返回欧氏距离
}

TEST_CASE("Question 2", "The rotation using phi-omega-kappa system") {
    Matrix3d R_student = question2();
    Matrix3d R_golden = golden2();

    REQUIRE(R_student.determinant() == Catch::Approx(1.0)); //行列式是否为1
    REQUIRE((R_student - R_golden).norm() == Catch::Approx(0.0));
}

TEST_CASE("Question 3", "Conversion between pok and opk") {
    Vector3d r_student = question3();
    Vector3d r_golden = golden3();

    REQUIRE((r_student - r_golden).norm() == Catch::Approx(0.0));
}