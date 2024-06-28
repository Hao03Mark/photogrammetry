//用新版catch2的V3喽，不用V2了
//#include <catch2/catch_test_macros.hpp> //这句也可以不需要
#include <catch_amalgamated.hpp>
#include <Eigen/Dense>

using namespace Eigen;

// 返回偏导数值
double question1();
double golden1();

// 返回拟合参数 a,b
std::tuple<double, double> question2();
std::tuple<double, double> golden2();

// 返回外方位元素线元素、旋转矩阵、单位权中误差和线元素理论精度
std::tuple<Vector3d, Matrix3d, double, Vector3d> question4();
std::tuple<Vector3d, Matrix3d, double, Vector3d> golden4();

TEST_CASE("Question 1", "Compute partial derivative") {
    double student = question1();
    double golden = golden1();

    //新版本Approx要加Catch::
    REQUIRE(std::abs(student - golden) == Catch::Approx(0.0).margin(1e-6));
}

TEST_CASE("Question 2", "Least squares fitting of curve") {
    double a_student, b_student;
    std::tie(a_student, b_student) = question2();

    double a_golden, b_golden;
    std::tie(a_golden, b_golden) = golden2();

    REQUIRE((a_student - a_golden) == Catch::Approx(0.0).margin(0.001));
    REQUIRE((b_student - b_golden) == Catch::Approx(0.0).margin(0.001));
}

TEST_CASE("Question 4", "Space resection") {
    Vector3d Xs_student, Xs_golden;
    Matrix3d R_student, R_golden;
    double sigma_student, sigma_golden;
    Vector3d m_Xs_student, m_Xs_golden;

    std::tie(Xs_student, R_student, sigma_student, m_Xs_student) = question4();
    std::tie(Xs_golden, R_golden, sigma_golden, m_Xs_golden) = golden4();

    REQUIRE(R_student.determinant() == Catch::Approx(1.0).margin(0.001)); 
    REQUIRE((Xs_student - Xs_golden).norm() == Catch::Approx(0.0).margin(0.1)); 
    REQUIRE((R_student - R_golden).norm() == Catch::Approx(0.0).margin(0.001)); 
    REQUIRE((sigma_student - sigma_golden) == Catch::Approx(0.0).margin(0.1)); 
    REQUIRE((m_Xs_student - m_Xs_golden).norm() == Catch::Approx(0.0).margin(0.1));
}
