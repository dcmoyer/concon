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
#include "Vector4.h"
//#include "Num_compute.h"

using namespace MeshLib;

Vector4::Vector4() 
{
    Vector4(0, 0, 0, 0);
};

Vector4::Vector4(float x, float y, float z, float w)
{
    _data[0]=x; _data[1]=y; _data[2]=z; _data[3]=w;
};

Vector4::Vector4(Point p, float w)
{
    _data[0]=(float)p(0); _data[1]=(float)p(1); _data[2]=(float)p(2); _data[3]=w;
};

Vector4::Vector4(Point p)
{
	_data[0] = (float)p(0); _data[1] = (float)p(1); _data[2] = (float)p(1); _data[3] = (float)p(2);
};

Vector4&
Vector4:: operator+=(const Vector4& v)
{
    _data[0]+=v._data[0]; _data[1]+=v._data[1]; _data[2]+=v._data[2]; _data[3]+=v._data[3];
    return *this;
};

Vector4& 
Vector4::operator-=(const Vector4& v)
{
    _data[0]-=v._data[0]; _data[1]-=v._data[1]; _data[2]-=v._data[2]; _data[3]-=v._data[3];
    return *this;
};

Vector4&
Vector4:: operator*=(float f)
{
    _data[0]*=f; _data[1]*=f; _data[2]*=f; _data[3]*=f; 
    return *this;
};

float 
Vector4::operator*(Vector4 & v)//used to compute (1x4)*(4x1)
{
    return _data[0]*v._data[0]+ _data[1]*v._data[1]+ _data[2]*v._data[2]+ _data[3]*v._data[3];
};

//Vector4	&
//Vector4::operator^(Vector4 & v)//used to compute 2x2 matix product	 , M = | v[0] v[2] |
//																		// | v[1] v[3] |
//{
//	double v1 = _data[0] * v._data[0] + _data[2] * v._data[1];
//	double v2 = _data[1] * v._data[0] + _data[3] * v._data[1];
//	double v3 = _data[0] * v._data[2] + _data[2] * v._data[3];
//	double v4 = _data[1] * v._data[2] + _data[3] * v._data[3];
//
//	Vector4 r(v1, v2, v3, v4);
//
//	return r;
//};

Vector4		 //probably much slower this way
Vector4::operator^(Vector4  v)//used to compute 2x2 matix product	 , M = | v[0] v[2] |
// | v[1] v[3] |
{
	double v1 = _data[0] * v._data[0] + _data[2] * v._data[1];
	double v2 = _data[1] * v._data[0] + _data[3] * v._data[1];
	double v3 = _data[0] * v._data[2] + _data[2] * v._data[3];
	double v4 = _data[1] * v._data[2] + _data[3] * v._data[3];

	Vector4 r(v1, v2, v3, v4);

	return r;
};



Vector4
Vector4::operator+(Vector4 & v)// adds vectors

{
	double v1 = _data[0] + v._data[0] ;
	double v2 = _data[1] + v._data[1] ;
	double v3 = _data[2] + v._data[2] ;
	double v4 = _data[3] + v._data[3] ;

	Vector4 r(v1, v2, v3, v4);

	return r;
};

Vector4
Vector4::operator-(Vector4 & v)// subtracts vectors

{
	double v1 = _data[0] - v._data[0];
	double v2 = _data[1] - v._data[1];
	double v3 = _data[2] - v._data[2];
	double v4 = _data[3] - v._data[3];

	Vector4 r(v1, v2, v3, v4);

	return r;
};

Vector4
Vector4::operator*(float & f)// scalar multiplication

{
	double v1 = _data[0] * f;
	double v2 = _data[1] * f;
	double v3 = _data[2] * f;
	double v4 = _data[3] * f;

	Vector4 r(v1, v2, v3, v4);

	return r;
};


Vector4&
Vector4::operator/=(float f)
{
    if ( f <= 1e-10 ) throw std::domain_error("divided by 0 error");
    _data[0] /= f;
	_data[1] /= f;
	_data[2] /= f;
	_data[3] /= f;
		
	return *this;
};

Point
Vector4::point()
{
    return Point(_data[0], _data[1], _data[2]);
};

Point
Vector4::Sym2x2()
{	   
	return Point(_data[0], _data[1], _data[3]);
};

Vector4
Vector4::Transpose2x2()
{
	return Vector4(_data[0], _data[2], _data[1],_data[3]);
};

double 
Vector4::Det2x2()
{
	return _data[0] * _data[3] - _data[1] * _data[2];
}

double
Vector4::Trace2x2()
{
	return _data[0] + _data[3] ;
}

/*

Vector4
Vector4::DInvSymMat(Vector4 dM)
{
	double det = _data[0] * _data[3] - _data[1] * _data[2];

	if (det == 0.0)
	{
		printf("Vector4 InvMat: 2x2 matrix NOT invertible\n");
		exit(1);
	}
//	printf("I get here, det = %f\n", det);


	Vector4 r1(_data[3] , (-1.0)*_data[1] , (-1.0)*_data[2] , _data[0] );
	Vector4 r2(dM[3] / det, (-1.0)*dM[1] / det, (-1.0)*dM[2] / det, dM[0] / det);
//	printf("I get here 3\n");
//	r2.output();
//	r2 /= det;

//	printf("I get here 4\n");

	double dDet = dM[0] * _data[3] + dM[3] * _data[0] - _data[2] * dM[1] - _data[1] * dM[2];
	r1 *= (float)(dDet/(det*det));

	r2 -= r1;

	return r2;
}

Vector4
Vector4::DLogSymMat_num(Point dM)
{
	Vector4 res(0.0, 0.0, 0.0, 0.0);
	Vector4 res1(0.0, 0.0, 0.0, 0.0);
	Vector4 res2(0.0, 0.0, 0.0, 0.0);
	Vector4 eye(1.0, 0.0, 0.0, 1.0);
	Vector4 neg_eye(-1.0, 0.0, 0.0, -1.0);
	Vector4 dM_V(dM);
	Vector4 tmp;
	Vector4 inv_tmp;
//	Vector4 copy(Sym2x2());

	int samples = 50;	   //at 100, this is 20% of the computation for Metric Mismatch

	//original
	for (int i = 0; i < samples; i++)
	{
		float s = (double)i / (double)samples;
		tmp = (*this + neg_eye);
		tmp *= s;
		tmp += eye;
		inv_tmp = tmp.InvMat();

		res += inv_tmp^dM_V^inv_tmp;	
	
	}


	res /= (float)(samples - 1);


	//simpson's rule

	//for (int i = 1; i < samples/2; i++)
	//{
	//	float s = (double)i*2 / (double)samples;
	//	Vector4 tmp = (*this + neg_eye);
	//	tmp *= s;
	//	tmp += eye;
	//	Vector4 inv_tmp = tmp.InvMat();

	//	Vector4 plus = inv_tmp^dM_V^inv_tmp;
	//	res1 += plus;

	//}

	//for (int i = 1; i <= samples / 2; i++)
	//{
	//	float s = (double)(i * 2 - 1) / (double)samples;
	//	Vector4 tmp = (*this + neg_eye);
	//	tmp *= s;
	//	tmp += eye;
	//	Vector4 inv_tmp = tmp.InvMat();

	//	Vector4 plus = inv_tmp^dM_V^inv_tmp;
	//	res2 += plus;

	//}

	//Vector4 F0 = dM_V;
	//Vector4 Fn = *this^dM_V^*this;
	//res1 *= 2.0;
	//res2 *= 4.0;

	//res = F0 + res1 + res2 + Fn;  
	//res /= (float)(3*(samples-1));

	return res;
}
*/

double 
Vector4::induced_norm(double v1, double v2)//meant for vector conjugation of 2x2 matrices	 x^t M x
{
	double x1 = _data[0] * v1 + _data[2] * v2;
	double x2 = _data[1] * v1 + _data[3] * v2;

	return v1*x1 + v2*x2;
}

/*
Vector4
Vector4::DLogSymMat(Point dM, Point spec1, Point spec2)
{
	//if (_data[1] != _data[2])
	if (fabs(_data[1] - _data[2]) / (fabs(_data[1]) + fabs(_data[2]) + 1e-8) > 1e-6)
	{
		printf("Vector4 LogSymMat: 2x2 matrix NOT symmetric\n");
		exit(1);
	}

	double l1 = spec1[0];
	double l2 = spec2[0];
	double x1 = spec1[1];
	double y1 = spec1[2];
	double x2 = spec2[1];
	double y2 = spec2[2];

	double Dl1;
	double Dl2;
	double Dx1;
	double Dy1;
	double Dx2;
	double Dy2;

	Num_compute nc;

	//bool diag = nc.D_eigen2(_data[0], _data[1], _data[3], dM[0], dM[1], dM[2],
	//	l1, l2, x1, y1, x2, y2,		Dl1, Dl2, Dx1, Dy1, Dx2, Dy2);

	bool diag = nc.eigen(_data[0], _data[1], _data[3],
		l1, l2, x1, y1, x2, y2);

	if (!diag)
	{
		printf("Vector4 LogSymMat: 2x2 matrix NOT diagonalizeable\n");
		exit(1);
	}


	if (l1 < 0.0 || l2 < 0.0)
	{
		printf("Vector4 LogSymMat: matrix NOT positive definite\n");
		exit(1);
	}

	//first, deriv of eigenvectors
	Vector4 V(x1, y1, x2, y2);
	Vector4 dM_vec(dM);

	Vector4 eye1(l1, 0.0, 0.0, l1);
	Vector4 dif1 = eye1 - *this;
	dif1 = dif1.Transpose2x2() ^ dM_vec;
	Dx1 = dif1[0] * x1 + dif1[2] * y1;
	Dy1 = dif1[1] * x1 + dif1[3] * y1;

	Vector4 eye2(l2, 0.0, 0.0, l2);
	Vector4 dif2 = eye2 - *this;
	dif2 = dif2.Transpose2x2() ^ dM_vec;
	Dx2 = dif2[0] * x2 + dif2[2] * y2;
	Dy2 = dif2[1] * x2 + dif2[3] * y2;

	Vector4 dV(Dx1, Dy1, Dx2, Dy2);


	//next, deriv of eigenvalues
	Dl1 = dM_vec.induced_norm(x1, y1);
	Dl2 = dM_vec.induced_norm(x2, y2);	

	Vector4 L(log(l1), 0.0, 0.0, log(l2));
	Vector4 dL(Dl1 / l1, 0.0, 0.0, Dl2 / l2); 

	Vector4 Vinv = V.InvMat();
	Vector4 dVinv = V.DInvSymMat(dV);	

	Vector4 dLog = dV^L^Vinv;
	dLog += V^dL^Vinv;
	dLog += V^L^dVinv; 

	return dLog;

};

void
Vector4::DLogSymMat(Point dM, Vector4 & Log, Vector4 & dLog)
{
	//if (_data[1] != _data[2])
	if (fabs(_data[1] - _data[2]) / (fabs(_data[1]) + fabs(_data[2]) + 1e-8) > 1e-6)
	{
		printf("Vector4 LogSymMat: 2x2 matrix NOT symmetric\n");
		exit(1);
	}

	double l1;
	double l2;
	double x1;
	double y1;
	double x2;
	double y2;

	double Dl1;
	double Dl2;
	double Dx1;
	double Dy1;
	double Dx2;
	double Dy2;

	Num_compute nc;

	//bool diag = nc.D_eigen2(_data[0], _data[1], _data[3], dM[0], dM[1], dM[2],
	//	l1, l2, x1, y1, x2, y2,		Dl1, Dl2, Dx1, Dy1, Dx2, Dy2);

	bool diag = nc.eigen(_data[0], _data[1], _data[3], 
		l1, l2, x1, y1, x2, y2);

	if (!diag)
	{
		printf("Vector4 LogSymMat: 2x2 matrix NOT diagonalizeable\n");
		exit(1);
	}


	if (l1 < 0.0 || l2 < 0.0)
	{
		printf("Vector4 LogSymMat: matrix NOT positive definite\n");
		exit(1);
	}

	//first, deriv of eigenvectors
	Vector4 V(x1, y1, x2, y2);	  	
	Vector4 dM_vec(dM);

	Vector4 eye1(l1, 0.0, 0.0, l1);
	Vector4 dif1 = eye1 - *this;
	dif1 = dif1.Transpose2x2() ^ dM_vec;
	Dx1 = dif1[0] * x1 + dif1[2] * y1;
	Dy1 = dif1[1] * x1 + dif1[3] * y1;

	Vector4 eye2(l2, 0.0, 0.0, l2);
	Vector4 dif2 = eye2 - *this;
	dif2 = dif2.Transpose2x2() ^ dM_vec;
	Dx2 = dif2[0] * x2 + dif2[2] * y2;
	Dy2 = dif2[1] * x2 + dif2[3] * y2;	  

	Vector4 dV(Dx1, Dy1, Dx2, Dy2);

	Dl1 = dM_vec.induced_norm(x1, y1);
	Dl2 = dM_vec.induced_norm(x2, y2); 


	Vector4 L(log(l1), 0.0, 0.0, log(l2));
	Vector4 dL(Dl1/l1, 0.0, 0.0, Dl2/l2);

//	if (V.Det2x2() < 0.0)
//		V = Vector4(x1, y1, (-1.0)*x2, (-1.0)*y2);

	Vector4 Vinv = V.InvMat();
	Vector4 dVinv = V.DInvSymMat(dV);

	//printf("dv:\n");
	//dV.output();
	//printf("V:\n");
	//V.output();

	//printf("Vinv:\n");
	//Vinv.output();

	Log = V^L^Vinv;

	dLog = dV^L^Vinv;
	dLog += V^dL^Vinv;
	dLog += V^L^dVinv;

	V.output();
	L.output();
	V.InvMat().output();

	//testing Num_compute

	//Vector4 Ltest((l1), 0.0, 0.0, (l2));
	//printf("\n\n printing reconstruction from spectral decomp\n\n:");
	//(V^Ltest^V.InvMat()).output();
	//printf("\n\n");

	// end test

//	return V^L^V.InvMat();

};


Vector4
Vector4::LogSymMat()
{
	//if (_data[1] != _data[2])
	if (fabs(_data[1] - _data[2]) / (fabs(_data[1]) + fabs(_data[2]) + 1e-8) > 1e-6)
	{
		printf("Vector4 LogSymMat: 2x2 matrix NOT symmetric\n");
		exit(1);
	}

	double l1;
	double l2;
	double x1;
	double y1;
	double x2;
	double y2;

	Num_compute nc;

	bool diag =  nc.eigen(_data[0], _data[1], _data[3], l1, l2, x1, y1, x2, y2);

	if (!diag)
	{
		printf("Vector4 ExpSymMat: 2x2 matrix NOT diagonalizeable\n");
		exit(1);
	}


	if (l1 < 0.0 || l2 < 0.0)
	{
		printf("Vector4 LogSymMat: matrix NOT positive definite\n");
		exit(1);
	}

	Vector4 V(x1, y1, x2, y2);	
	Vector4 L(log(l1), 0.0, 0.0, log(l2));
  
	//	return V^L^V.InvMat(); works for VS compiler

	Vector4 Vinvmat = V.InvMat();
	Vector4 intermid = L^Vinvmat;

	Vector4 res = V^intermid;

	return res;
	 
};

Vector4
Vector4::LogSymMat(Point & spec1, Point & spec2)
{
	//if (_data[1] != _data[2])
	if (fabs(_data[1] - _data[2]) / (fabs(_data[1]) + fabs(_data[2]) + 1e-8) > 1e-6)
	{
		printf("Vector4 LogSymMat: 2x2 matrix NOT symmetric\n");
		exit(1);
	}

	double l1;
	double l2;
	double x1;
	double y1;
	double x2;
	double y2;

	Num_compute nc;

	bool diag = nc.eigen(_data[0], _data[1], _data[3], l1, l2, x1, y1, x2, y2);

	spec1[0] = l1;
	spec1[1] = x1;
	spec1[2] = y1;

	spec2[0] = l2;
	spec2[1] = x2;
	spec2[2] = y2;

	if (!diag)
	{
		printf("Vector4 ExpSymMat: 2x2 matrix NOT diagonalizeable\n");
		exit(1);
	}


	if (l1 < 0.0 || l2 < 0.0)
	{
		printf("Vector4 LogSymMat: matrix NOT positive definite\n");
		exit(1);
	}

	Vector4 V(x1, y1, x2, y2);
	Vector4 L(log(l1), 0.0, 0.0, log(l2));

//	return V^L^V.InvMat(); works for VS compiler

	Vector4 Vinvmat = V.InvMat();
	Vector4 intermid = L^Vinvmat;

	Vector4 res = V^intermid;

	return res;

};

Vector4
Vector4::ExpSymMat()
{
	if (fabs(_data[1] - _data[2]) / (fabs(_data[1]) + fabs(_data[2]) + 1e-8) > 1e-6)
	{
		printf("Vector4 ExpSymMat: 2x2 matrix NOT symmetric\n");
		exit(1);
	}

	if (_data[0] == 0.0 &&_data[1] == 0.0 &&_data[3] == 0.0)
		return Vector4(1.0, 0.0, 0.0, 1.0);

	double l1;
	double l2;
	double x1;
	double y1;
	double x2;
	double y2;

	Num_compute nc;

	bool diag = nc.eigen(_data[0], _data[1], _data[3], l1, l2, x1, y1, x2, y2);

	if (!diag)
	{
		printf("Vector4 ExpSymMat: 2x2 matrix NOT diagonalizeable\n");
		exit(1);
	}
		

	Vector4 V(x1, y1, x2, y2);
	Vector4 L(exp(l1), 0.0, 0.0, exp(l2));

	Vector4 Vinvmat = V.InvMat();
	Vector4 intermid = L^Vinvmat;

	Vector4 res = V^intermid;

	return res;

};

Vector4
Vector4::InvMat()
{
	double det = _data[0] * _data[3] - _data[1] * _data[2];

	if (det == 0.0)
	{
		printf("Vector4 InvMat: 2x2 matrix NOT invertible\n");
		exit(1);
	}

	Vector4 r(_data[3] / det, (-1.0)*_data[1] / det, (-1.0)*_data[2] / det, _data[0] / det);

//	r.output();
//	printf("det  = %lf\n", det);
 
	return r;

}

*/
void 
Vector4::output()
{
/*
    std::cout<< "Vector4(" << _data[0] << "," << _data[1] <<
                 "," << _data[2] << "," << _data[3] << ")\n";
*/
    printf("Vector4: %g %g %g %g\n", _data[0], _data[1], _data[2], _data[3] );
};


float& 
Vector4::operator[](int i) { return _data[i];};

const float& 
Vector4::operator[](int i) const { return _data[i];};

