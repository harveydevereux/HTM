#include "../src/HTM.h"
#include <iostream>
#include <fstream>
#include <ctime>
const int nBuilds = 11;

int main(){
  vector< vector<double> > data;
  vector<double> times(nBuilds,0.0);
  for (int i = 0; i < nBuilds; i++){
    HTM htm;
    clock_t start;
    start = clock();
    htm.build(i);
    times[i] = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Ellapsed Time: " << times[i] << endl;
    vector< Trixel <double> > leaves = htm.leaves();
    vector<double> areas(leaves.size(),0.0);
    for (int j = 0; j < leaves.size(); j++){
      areas[j] = leaves[j].area();
    }
    data.push_back(areas);
  }

  ofstream output;
  output.open("areaData.txt");
  if (output.is_open()){
    for (int i = 0; i < data.size(); i++){
      for (int j = 0; j < data[i].size(); j++){
        output << data[i][j] << ",";
      }
      output << endl;
    }
  }
  output.close();
  output.open("runtimeData.txt");
  if (output.is_open()){
    for (int i = 0; i < times.size(); i++){
      output << times[i] << ",";
    }
  }
  output.close();

  for (int i = 0; i < 3; i++){
    HTM htm;
    htm.build(i);
    vector< Trixel <double> > leaves = htm.leaves();
    for (int j = 0; j < leaves.size(); j++){
      cout << leaves[j].getID() << ", ";
    }
    cout << endl;
  }
  return 0;
}
