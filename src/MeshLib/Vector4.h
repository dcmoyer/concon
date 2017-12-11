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

//Vertor4

#ifndef _VECTOR4_H_
#define _VECTOR4_H_

#include <cstdio>
#include "Point.h"

namespace MeshLib
{

    class Vector4 
    {
    public:
		Vector4() /*{_data[0] = _data[1] = _data[2] = _data[3] = 0.0;}*/;

        Vector4(float a);

        Vector4(float x, float y, float z, float w);

        Vector4(Point p, float w);	 //3-vector, scalar	  = 4-vector

		Vector4(Point p);				// symmetric 2x2 matrix


        Vector4& operator+=(const Vector4& v);

        Vector4& operator-=(const Vector4& v);


        Vector4& operator*=(float f);

        float operator*(Vector4 & v);//used to compute (1x4)*(4x1)

	//	Vector4 &  operator^(Vector4 & v);//used to compute 2x2 matrix product
		Vector4   operator^(Vector4  v);//used to compute 2x2 matrix product
		Vector4	 operator+(Vector4 & v);
		Vector4	 operator-(Vector4 & v);
		Vector4	 operator*(float & f);

        Vector4& operator/=(float f);

        Point point();

		Point Sym2x2();	 //returns 2x2 symmetric matrix in point form
		Vector4 Transpose2x2(); //returns transpose of 2x2 matrix

		//Vector4 LogSymMat();
		//Vector4 LogSymMat(Point & spec1, Point & spec2); //also returns spectral decomposition

		//void DLogSymMat(Point dM, Vector4 & Log, Vector4 & dLog); 
		//Vector4 DLogSymMat(Point dM, Point spec1, Point spec2); //takes pre-computed spectral decomposition
		//Vector4 DLogSymMat_num(Point dM);

		//Vector4 ExpSymMat();

		//Vector4 InvMat();
		//Vector4 DInvSymMat(Vector4 dM);

		//meant for vector conjugation of 2x2 matrices	 x^t M x
		double induced_norm(double v1, double v2);

		double Det2x2();

		double Trace2x2();

        void output();



        float& operator[](int i) ;
        const float& operator[](int i) const;

        friend inline Vector4 operator+(const Vector4& v1, const Vector4& v2);

    private:

        float _data[4];
    };

    inline Vector4 
    operator+(const Vector4& v1, const Vector4& v2)
    {

        float r[4];

        for ( int i = 0; i < 4; i ++ )
        {
            r[i] = v1[i] + v2[i];
        }

        return Vector4(r[0],r[1],r[2],r[3]);


    };




}//end of namespace MeshLib
#endif // Vector4

