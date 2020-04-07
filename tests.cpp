#include "HTM.h"
#include <random>
#include <chrono>
const double tol = 1e-6;

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

SCENARIO("Testing if rays pass through triangles", "[rayIn3DTriangle]"){
  GIVEN("a triangle with vertices [1 0 0], [0 1 0], [0 0 1]"){
    vector<int> x {1,0,0};
    vector<int> y {0,1,0};
    vector<int> z {0,0,1};
    Triangle <int> tri(x,y,z,x.size());

    AND_GIVEN("the points p = [0,0,0], q = [0.5,0.5,0]"){
      vector<double> p {0,0,0};
      vector<double> q {0.5,0.5,0};
      THEN("the ray from p to q passes through the triangle"){
        REQUIRE(rayIn3DTriangle <double, int>(p,q,tri));
      }
      THEN("the ray from p to -q does not pass through the triangle"){
        vector<double> k {-0.5,-0.5,0};
        REQUIRE(!(rayIn3DTriangle <double, int>(p,k,tri)));
      }
    }
  }
}
