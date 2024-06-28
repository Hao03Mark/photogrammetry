#include <Eigen/Core>
#include <iostream> // cmath已经包含在iostream中，单独使用cmath的话还会提示报错
//#include <cmath>
//#include <autodiff/forward/forward.hpp>
#include <autodiff/forward/dual.hpp>

using namespace autodiff;

double question1(){
    // 计算函数 f = exp(x * x + sqrt(y) / log(x * y)) 对 x 的偏导数
    double x = 1.0;
    double y = 3.0;
    
    return exp(x * x + sqrt(y) / log(x * y)) * (2 * x - sqrt(y)/pow(log(x * y), 2) / x);

}

dual f(dual x, dual y){
    return exp(x * x + sqrt(y) / log(x * y));
}

double golden1(){
    dual x = 1.0;
    dual y = 3.0;
    dual u = f(x,y);
    double dudx = derivative(f, wrt(x), at(x,y));
    return dudx;
}
