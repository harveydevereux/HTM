#include "HTM.h"
#include <iostream>

int main(){
  vector<int> x {1,2,3};
  vector<int> y {4,5,6};
  vector<int> z {7,8,9};

  Trixel <int> trix("",x,y,z,x.size());
  vector< vector<int> > v = trix.getVertices();
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      cout << v[i][j] << ", ";
    }
    cout << endl;
  };

  HTM htm;
  htm.build();
  cout << htm.size() << endl;
  htm.build(2);
  cout << htm.size() << endl;
  //vector< Trixel <double> > leaves = htm.leaves();
  return 0;

}
