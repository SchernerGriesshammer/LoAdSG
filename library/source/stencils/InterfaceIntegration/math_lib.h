/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer, Rainer Hartmann
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/

// ------------------------------------------------------------
// math_lib.h
//
// ------------------------------------------------------------

#ifndef MATHLIB_H_
#define MATHLIB_H_


//////////////////////////////////////////////////////////////////////
// a simple 3D vector class
//////////////////////////////////////////////////////////////////////

class D3vector {
 public:
  double x,y,z;
  D3vector(double cx, double cy, double cz) : x(cx), y(cy), z(cz) {};
  explicit D3vector(double c) : x(c), y(c), z(c) {};
  D3vector() : x(0), y(0), z(0) {};
  ~D3vector(){};
  void Print();
  void Print(ofstream *Datei);
  void operator=(const D3vector& v) { x=v.x; y=v.y; z=v.z; }
  double operator[](int i) {
    if(i==0) return x;
    if(i==1) return y;
    return z;
  }
};

inline double MIN(D3vector V) {
  if(V.x<V.z && V.x<V.y) return V.x;
  if(V.y<V.z && V.y<V.x) return V.y;
    return V.z;
}

inline double MAX(D3vector V) {
  if(V.x>V.z && V.x>V.y) return V.x;
  if(V.y>V.z && V.y>V.x) return V.y;
    return V.z;
}

/*
inline D3vector MAX(D3vector V1, D3vector V2) {
  return D3vector(MAX(V1.x,V2.x),MAX(V1.y,V2.y),MAX(V1.z,V2.z));
}

inline D3vector MIN(D3vector V1, D3vector V2) {
  return D3vector(MIN(V1.x,V2.x),MIN(V1.y,V2.y),MIN(V1.z,V2.z));
}
*/

inline D3vector operator+(const D3vector& v,const D3vector& w) {
  return D3vector(v.x+w.x,v.y+w.y,v.z+w.z);
}

inline D3vector operator*(const D3vector& v,const D3vector& w) {
  return D3vector(v.x*w.x,v.y*w.y,v.z*w.z);
}

inline D3vector operator/(const D3vector& v,const D3vector& w) {
  return D3vector(v.x/w.x,v.y/w.y,v.z/w.z);
}

inline D3vector operator-(const D3vector& v,const D3vector& w) {
  return D3vector(v.x-w.x,v.y-w.y,v.z-w.z);
}

inline D3vector operator/(const D3vector& v,const double f) {
  return D3vector(v.x/f,v.y/f,v.z/f);
}

inline D3vector operator*(const D3vector& v,const double f) {
  return D3vector(v.x*f,v.y*f,v.z*f);
}

inline D3vector operator*(const double f, const D3vector& v) {
  return D3vector(v.x*f,v.y*f,v.z*f);
}

inline bool operator<(const D3vector& v,const D3vector& w) {
  return (v.x<w.x && v.y<w.y && v.z<w.z);
}
inline bool operator==(const D3vector& v,const D3vector& w) {
  return (v.x==w.x && v.y==w.y && v.z==w.z);
}


//////////////////////////////////////////////
// Implementierung einiger Memberfunktionen
//////////////////////////////////////////////


// D3vector 
// ----------

inline void D3vector::Print() {
  cout << "Coordinate: " << x << ", " << y << ", " << z << ";";
}

inline void D3vector::Print(ofstream *Datei) {
  *Datei  << x << " " << y << " " << z;
}


#endif // MATHLIB_H_


