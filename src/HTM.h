// basic implementation of a Hierachical triangular mesh
// H. Devereux, University of Warwick, 2020
//
// Indexing the Sphere with the Hierarchical Triangular Mesh:
// https://arxiv.org/pdf/cs/0701164.pdf
#include "Trixel.h"
using namespace std;
#include <iostream>
#include <fstream>
#include <string>
class HTM{
private:
  size_t depth;
  vector< Trixel <double> > Mesh; // vector of the trixels
  vector<int> leafIndices(){                                                    // returns the integer indexes which correspond
    vector<int> leaves;                                                         // to the current leaf nodes in the HTM
    vector<int> stack {0,1,2,3,4,5,6,7};                                        // the recursive search is marginally slower
    while (stack.size() > 0){                                                   // than using the fact the leaves are the
      int ptr = stack[stack.size()-1];                                          // Mesh.size()-8*pow(4,depth-1) to Mesh.size() nodes
      stack.pop_back();
      Trixel <double> trix = this->Mesh[ptr];
      if ( (trix.getChildren()).size() == 0 ){
        leaves.push_back(ptr);
        continue;
      }
      bool cont = false;
      for (int i = 0; i < (trix.getChildren()).size(); i++){
        if ((trix.getChildren())[i] == -1){
          leaves.push_back(ptr);
          cont = true;
          break;
        }
      }
      if (cont){
        continue;
      }
      else{
        for (int i = 0; i < (trix.getChildren()).size(); i++){
          stack.push_back((trix.getChildren())[i]);
        }
      }
    }
    return leaves;
  }
public:
  HTM()
  : depth(0), Mesh(vector< Trixel <double> >(0))
  {
    // depth 0 hard coded
    vector <double> v0 {0.0,0.0,1.0};
    vector <double> v1 {1.0,0.0,0.0};
    vector <double> v2 {0.0,1.0,0.0};
    vector <double> v3 {-1.0,0.0,0.0};
    vector <double> v4 {0.0,-1.0,0.0};
    vector <double> v5 {0.0,0.0,-1.0};

    Trixel <double> T1("S0",v1,v5,v2,3);
    Trixel <double> T2("N0",v1,v0,v4,3);
    Trixel <double> T3("S1",v2,v5,v3,3);
    Trixel <double> T4("N1",v4,v0,v3,3);
    Trixel <double> T5("S2",v3,v5,v4,3);
    Trixel <double> T6("N2",v3,v0,v2,3);
    Trixel <double> T7("S3",v4,v5,v1,3);
    Trixel <double> T8("N3",v2,v0,v1,3);
    vector< Trixel <double> > u {T1,T2,T3,T4,T5,T6,T7,T8};
    this->Mesh = u;
    this->depth = depth;
  };
  size_t size(){
    return this->depth;
  }
  void build(size_t depth = 0){                                                 // build the HTM from scrath
    if (depth == 0){
      return;
    }
    for (int i = 0; i < depth; i++){
      int start = this->Mesh.size()-8*pow(4,i);
      int end = this->Mesh.size();
      for (int j = start; j < end; j++){
        vector< Trixel <double> > newTrixels = subdivideTrixel <double> (this->Mesh[j]);
        vector<int> children(4,0);
        for (int k = 0; k < newTrixels.size(); k++){
          newTrixels[k].setParent(j);
          this->Mesh.push_back(newTrixels[k]);
          children[k] = this->Mesh.size()-1;
        }
        this->Mesh[j].setChildren(children);
      }
    }
    return;
  }
  vector< Trixel <double> > leaves(){                                           // return a vector of the leaf trixels
    vector< Trixel <double> > l;
    vector<int> indices = this->leafIndices();
    for (int i = 0; i < indices.size(); i++){
      l.push_back(this->Mesh[i]);
    }
    return l;
  }
  void writeLeaves(const string & filename){                                    // output a csv where each trixel is 9 floats followed
    ofstream output;                                                            // by it's name
    output.open(filename);
    if (output.is_open()){
      vector<int> leaves = this->leafIndices();
      for (int i = 0; i < leaves.size(); i++){
        Trixel <double> trix = Mesh[leaves[i]];
        vector< vector<double> > v = trix.getVertices();
        for (int j = 0; j < v.size(); j++){
          for (int k = 0; k < v[j].size(); k++){
            output << v[j][k] << ",";
          }
        }
        output << trix.getID() << endl;
      }
      output.close();
    }
    else{
      cout << "Error in reading output file: " << filename << endl;
      exit (EXIT_FAILURE);
    }
  }
};
