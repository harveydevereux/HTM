#include "../src/HTM.h"
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

SCENARIO("building an HTM of depth d"){
  GIVEN("depths of 0,1,2,3,4 and 5"){
    vector<int> depth {0,1,2,3,4,5};
    THEN("the HTM should have 8x4^(d) Trixels: [8,32,128,512,2048,8192]"){
      vector<int> theory {8,32,128,512,2048,8192};
      bool correct = true;
      for (int i = 0; i < depth.size(); i++){
        HTM htm;
        htm.build(depth[i]);
        vector< Trixel <double> > leaves = htm.leaves();
        if (leaves.size() != theory[i]){
          correct = false;
          break;
        }
      }
      REQUIRE(correct);
    }
  }
}

SCENARIO("triangle area"){
  GIVEN("a triangle with vertices [0,0,0],[1,0,0], and [0,1,0]"){
    vector<double> x {0,0,0};
    vector<double> y {0,1,0};
    vector<double> z {1,0,0};
    Trixel <double> t("",x,y,z,3);
    THEN("the are is 0.5"){
      REQUIRE(t.area()==0.5);
    }
  }
}
