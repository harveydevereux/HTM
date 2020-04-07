// basic implementation of a Hierachical triangular mesh
// H. Devereux, University of Warwick, 2020
//
// Indexing the Sphere with the Hierarchical Triangular Mesh:
// https://arxiv.org/pdf/cs/0701164.pdf
#include <vector>
#include <string>
#include <cmath>
using namespace std;

template <class T>
class Triangle {
protected:
  size_t dimension;
  vector<T> x;
  vector<T> y;
  vector<T> z;
public:
  Triangle()
  {
    vector<T> v1{0,0,0};
    vector<T> v2{0,0,0};
    vector<T> v3{0,0,0};
    x = v1;
    y = v2;
    z = v3;
    dimension = 3;
  };
  Triangle(vector<T> v1,vector<T> v2,vector<T> v3,size_t dimension)
  : x(v1),y(v2),z(v3),dimension(dimension)
  {};
  vector< vector<T> > getVertices(){
    vector< vector<T> > v(3,vector<T>(this->dimension,0.0));
    for (int j = 0; j < dimension; j++){
      v[0][j] = this->x[j];
      v[1][j] = this->y[j];
      v[2][j] = this->z[j];
    }
    return v;
  }
  void setVertices(vector< vector<T> > v){
    for (int j = 0; j < dimension; j++){
      this->x[j] = v[0][j];
      this->y[j] = v[1][j];
      this->z[j] = v[2][j];
    }
  }
  size_t getDimension(){
    return this->dimension;
  }
};

template <class T>
class Trixel: public Triangle<T> {
private:
  string ID;
  Trixel <T> * parent;
  vector< Trixel <T> * > children;
public:
  Trixel()
  : ID(""), parent(nullptr), children(vector< Trixel <T> * >(4,nullptr))
  {
    vector<T> v1{0,0,0};
    vector<T> v2{0,0,0};
    vector<T> v3{0,0,0};
    this->x = v1;
    this->y = v2;
    this->z = v3;
    this->dimension = 3;
  };
  Trixel(string id, vector<T> x, vector<T> y, vector<T> z, size_t dimension, Trixel * parent = nullptr,children = vector< Trixel <T> * >(4,nullptr))
  : ID(""), parent(parent), children(children);
  {
    this->x = x;
    this->y = y;
    this->z = z;
    this->dimension = dimension;
  };
  string getID(){
    return this->ID;
  }
  void setParent(Trixel<T> trix){
    this->parent = &trix;
  }
  void setChildren(vector< Trixel<T> > trix){
    for (int i = 0; i < 4; i++){
      this->children[i] = &trix[i];
    }
  }
};

template <class T, class T2>
uint8_t rayIn3DTriangle(vector<T> p, vector<T> q, Triangle<T2> tri){
  double rayx = q[0]-p[0];
  double rayy = q[1]-p[1];
  double rayz = q[2]-p[2];
  vector< vector<T2> > vertices = tri.getVertices();

  double ux = vertices[1][0] - vertices[0][0];
  double uy = vertices[1][1] - vertices[0][1];
  double uz = vertices[1][2] - vertices[0][2];

  double vx = vertices[2][0] - vertices[0][0];
  double vy = vertices[2][1] - vertices[0][1];
  double vz = vertices[2][2] - vertices[0][2];

  double nx = uy*vz - uz*vy;
  double ny = -(ux*vz - uz*vx);
  double nz = ux*vy - uy*vx;

  double n = sqrt(rayx*rayx+rayy*rayy+rayz*rayz);

   double dirx = rayx/n;
   double diry = rayy/n;
   double dirz = rayz/n;

   double w0x = p[0] - vertices[0][0];
   double w0y = p[1] - vertices[0][1];
   double w0z = p[2] - vertices[0][2];

   double a = -(nx*w0x + ny*w0y + nz*w0z);
   double b = nx*dirx + ny*diry + nz*dirz;
   if (abs(b) < 1e-16){
     // technically if a==0 it is indeterminate
     return 0;
   }
   double r = a/b;
   if (r < 0.0){
     return 0;
   }
   double Ix = p[0] + r*dirx;
   double Iy = p[1] + r*diry;
   double Iz = p[2] + r*dirz;

   double uu = ux*ux + uy*uy + uz*uz;
   double uv = ux*vx + uy*vy + uz*vz;
   double vv = vx*vx + vy*vy + vz*vz;
   double wu = (Ix-vertices[0][0])*ux + (Iy-vertices[0][1])*uy + (Iz-vertices[0][2])*uz;
   double wv = (Ix-vertices[0][0])*vx + (Iy-vertices[0][1])*vy + (Iz-vertices[0][2])*vz;

   double D = uv*uv - uu*vv;
   double s = (uv*wv-vv*wu)/D;
   if (s < 0.0 || s > 1.0){
     return 0;
   }
   double t = (uv*wu-uu*wv)/D;
   if (t < 0.0 || (s+t) > 1.0){
     return 0;
   }
   return 1;
}

template <class T>
vector< Trixel <T> > subdivideTrixel(Trixel <T> trix){
  vector< vector<T> > v = trix.getVertices();
  size_t d = trix.getDimension();
  vector<double> w0(d,0);
  vector<double> w1(d,0);
  vector<double> w2(d,0);

  w0[0] = v[2][0]+v[1][0];
  w0[1] = v[2][1]+v[1][1];
  w0[2] = v[2][2]+v[1][2];

  w1[0] = v[2][0]+v[0][0];
  w1[1] = v[2][1]+v[0][1];
  w1[2] = v[2][2]+v[0][2];

  w2[0] = v[1][0]+v[0][0];
  w2[1] = v[1][1]+v[0][1];
  w2[2] = v[1][2]+v[0][2];

  double n0 = sqrt(w0[0]*w0[0]+w0[1]*w0[1]+w0[2]*w0[2]);
  double n1 = sqrt(w1[0]*w1[0]+w1[1]*w1[1]+w1[2]*w1[2]);
  double n2 = sqrt(w2[0]*w2[0]+w2[1]*w2[1]+w2[2]*w2[2]);

  w0[0] = w0[0]/n0;
  w0[1] = w0[1]/n0;
  w0[2] = w0[2]/n0;

  w1[0] = w1[0]/n1;
  w1[1] = w1[1]/n1;
  w1[2] = w1[2]/n1;

  w2[0] = w2[0]/n2;
  w2[1] = w2[1]/n2;
  w2[2] = w2[2]/n2;

  vector<T> v0 {v[0][0],v[0][1],v[0][2]};
  vector<T> v1 {v[1][0],v[1][1],v[1][2]};
  vector<T> v2 {v[2][0],v[2][1],v[2][2]};

  string id = trix.getID();

  Trixel <T> T1(id+"0",v0,vector<T>(w2),vector<T>(w1),d);
  Trixel <T> T2(id+"1",v1,vector<T>(w0),vector<T>(w2),d);
  Trixel <T> T3(id+"2",v2,vector<T>(w1),vector<T>(w0),d);
  Trixel <T> T4(id+"3",vector<T>(w0),vector<T>(w1),vector<T>(w2),d);

  return vector< Trixel <T> > u {T1,T2,T3,T4};
}

template <class T>
class HTM{
private:
  size_t size;
  vector< Trixel <T> > Mesh; // vector of the 8 foundational trixels
public:
  HTM()
  : size(0), Mesh(vector< Trixel <T> >(0))
  {};
  void BuildHTM(size_t depth){

  }
};
