// #define CATCH_CONFIG_MAIN
// #include <catch2/catch.hpp>
#include <catch_amalgamated.hpp>
#include "problem.h"


TEST_CASE("Question 1", "Space Intersection") {
    
    DataIntersection data = generate_intersection();

    Vector3d X_initial = intersection_similarity_transform(data);
    Vector3d X_student = intersection_least_squares(data, X_initial);

    REQUIRE((X_golden - X_initial).norm() == Catch::Approx(0.0).margin(10.0));
    REQUIRE((X_golden - X_student).norm() == Catch::Approx(0.0).margin(0.1));
}

TEST_CASE("Question 2", "Relative Orientation") {
    DataRelativeOrientation data = generate_relative_orientation();

    Vector3d B_student, B_golden;
    Matrix3d R_student, R_golden;

    std::tie(B_student, R_student) = relative_orientation(data);
    std::tie(B_golden, R_golden) = relative_orientation_golden(data);

    REQUIRE((B_student - B_golden).norm() == Catch::Approx(0.0).margin(0.1));
    REQUIRE((R_student - R_golden).norm() == Catch::Approx(0.0).margin(0.001));
}

TEST_CASE("Question 3", "Absolute Orientation") {
    DataAbsoluteOrientation data = generate_absolute_orientation();

    DataIntersection data_student = absolute_orientation(data);
    DataIntersection data_golden = generate_intersection();
    
    REQUIRE((data_student.X_eo[1] - data_golden.X_eo[1]).norm() == Catch::Approx(0.0).margin(0.1));
    REQUIRE((data_student.R_eo[1] - data_golden.R_eo[1]).norm() == Catch::Approx(0.0).margin(0.1));
}
