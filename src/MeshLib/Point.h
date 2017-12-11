//Copyright 2017 Yalin Wang and Boris Gutman
//
//Permission is hereby granted, free of charge, to any person obtaining a
//copy of this software and associated documentation files (the "Software"),
//to deal in the Software without restriction, including without limitation
//the rights to use, copy, modify, merge, publish, distribute, sublicense,
//and/or sell copies of the Software, and to permit persons to whom the
//Software is furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in
//all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
//OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
//OTHER DEALINGS IN THE SOFTWARE.
//

#ifndef _MESHLIB_POINT_H_

#define _MESHLIB_POINT_H_

//#ifndef pi
//#define pi 3.1415926535897932384626433832795
//#endif
//#ifndef M_PI
//#define M_PI 3.1415926535897932384626433832795
//#endif


#include <iostream>
#include <stdexcept>

//TODO: switch to std lib
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


namespace MeshLib{

class Point {

public:

  Point( double x, double y, double z ){ v[0] = x; v[1] = y; v[2] = z;};

  Point() { v[0] = v[1] = v[2] = 0.0;};

  //spherical coordinates
  Point(double theta, double phi){v[0] = sin(theta)*cos(phi); v[1] = sin(theta)*sin(phi); v[2] = cos(theta);};

  ~Point(){};

  double angle(Point  p);

  double SymDet(); //symmetric matrix operations
  double SymTrace();
  Point SymInverse();
  Point SymRoot();
  double SymNorm();
  double SymNorm2();

  Point SymInvRoot_derivative(Point dM);
  Point SymInverse_derivative(Point dM);
  double SymDet_derivative(Point dM);
  Point SymConjDeriv(int coord);


  double theta();
  double phi();
  double line_dist(Point p1, Point p2, Point closest, double * t);
  double line_dist(Point p1, Point p2, Point closest);

  void print();

  Point rotate(double theta, Point vector);//Rodriguez formula

  //parallel transport of the stored vector, from p_old to p_new
  Point PT(Point p_new, Point p_old, bool verbose = false);

  void EulerRotate(double alpha, double beta, double gamma) {
    if(gamma < 0.0)
      gamma += 2.0*M_PI;
    if(alpha < 0.0)
      alpha += 2.0*M_PI;

    double x = v[0];
    double y = v[1];
    double z = v[2];
    double sine = sin(gamma);
    double cose = cos(gamma);

    v[0]=cose* x - sine* y;
    v[1]=cose* y + sine* x;

    x = v[0];
    y = v[1];

    sine = sin(/*(-1.0)**/beta);
    cose = cos(/*(-1.0)**/beta);

    v[0]=cose* x - sine* z;
    v[2]=cose* z + sine* x;

    x = v[0];
    z = v[2];

    //print();

    sine = sin(alpha);
    cose = cos(alpha);

    v[0]=cose* x - sine* y;
    v[1]=cose* y + sine* x;
  }


  double & operator[]( int i ) { 
    if ( 0>i || i>=3 ) throw std::out_of_range("invalid index for point coordinates");
    return v[i]; 
  }

  double operator()( int i) const { 
    if ( 0>i || i>=3 ) throw std::out_of_range("invalid index for point coordinates");
    return v[i]; 
  }

  double operator[]( int i) const {
      if ( 0>i || i>=3 ) throw std::out_of_range("invalid index for point coordinates");
      return v[i]; 
  }

  double norm() const {
    return sqrt( fabs( v[0] * v[0] + v[1] * v[1] + v[2] * v[2] ) );
  }

  double norm2() const {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  }

  Point& operator += ( const Point & p) {
    v[0] += p(0); v[1] += p(1); v[2] += p(2);
    return *this;
  }

  Point& operator -= ( const Point & p) {
    v[0] -= p(0); v[1] -= p(1); v[2] -= p(2);
    return *this;
  }

  Point& operator *= ( const double  s) {
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
    return *this;
  }

  Point& operator /= ( const double  s) {
    v[0] /= s;
    v[1] /= s;
    v[2] /= s;
    return *this;
  }

  double operator *( const Point & p ) const {
    return v[0] * p[0] + v[1] * p[1] + v[2] * p[2];
  }

  //2x2 symmetric matrix multiplication
  Point operator %(const Point & p) const {         
    return Point(v[0] * p[0] + v[1] * p[1],  v[1] * p[0] + v[2]*p[1], 0.0);
  }

  Point   operator+( const Point & p  ) const {
    Point r( v[0] + p[0], v[1] + p[1], v[2] + p[2] );
    return r;
  }

  Point   operator-( const Point & p  ) const {
      Point r( v[0] - p[0], v[1] - p[1], v[2] - p[2] );
      return r;
  }

  Point   operator*( const double s  ) const {
      Point r( v[0] * s, v[1] * s, v[2] * s );
      return r;
  }

  Point   operator/( const double s  ) const {
      Point r( v[0] / s, v[1] / s, v[2] / s );
      return r;
  }

  Point operator^( const Point & p2) const {
    Point r( v[1] * p2[2] - v[2] * p2[1],
      v[2] * p2[0] - v[0] * p2[2],
      v[0] * p2[1] - v[1] * p2[0]);
    return r;
  }
        
  // 2x2 symmetric matrix operation 
  //Symmetric matrix conjugation,  
  //A = | v[0] v[1] |    ,  B =  | p2[0] p2[1] |  , returns B*A*B   
  //    | v[1] v[2] |            | p2[1] p2[2] |
  Point operator|(const Point & p2) const {
    Point E1(p2[0] * p2[0], 2.0*p2[0] * p2[1], p2[1] * p2[1]);
    Point E2(p2[0] * p2[1], p2[0] * p2[2] + p2[1] * p2[1], p2[1] * p2[2]);
    Point E3(p2[1] * p2[1], 2.0*p2[1] * p2[2], p2[2] * p2[2]);
    Point r(*this*E1,
      *this*E2,
      *this*E3);
    return r;
  };

  Point operator-() const {
      Point p(-v[0],-v[1],-v[2]);
      return p;
  };

  bool operator == (const Point & p) {
    return(v[0] == p(0) && v[1] == p(1) && v[2] == p(2));
  };

  //seriously?
  bool operator< (const Point & p) {
    bool returnValue = false;

    if ( v[0] < p(0) ) returnValue = true;
    else if ( v[0] == p[0] ) {

      if ( v[1] < p(1) ) returnValue = true;
      else if ( v[1] == p(1) ) {

        if ( v[2] < p(2) ) returnValue = true;
        else returnValue = false;
      }
      else returnValue = false;

    }
    else returnValue = false;

    return returnValue;
  };

protected:
  double v[3];

};



std::ostream & operator<<( std::ostream & os, const Point & p);



}//end of namespace MeshLib



#endif //_MESHLIB_POINT_H_ defiined

