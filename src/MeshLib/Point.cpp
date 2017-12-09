
#include "Point.h"

using namespace MeshLib;


std::ostream & operator<<( std::ostream & os, const Point & p)
{
    os << "Point: " << p(0) << " " << p(1) << " " << p(2) << std::endl;
    return os;
};

void Point::print() {printf("Point: %g %g %g\n",v[0],v[1],v[2]);fflush(stdout);};

Point Point::rotate(double theta, Point axis)//Rodriguez formula
{
    Point result;
	Point vector = axis/axis.norm();
    double cos_t = cos(theta);
    double sin_t = sin(theta);
    result = vector * (v[0] * vector[0] + v[1] * vector[1] + v[2] * vector[2]) * (1-cos_t);
    result[0] += v[0] * cos_t;
    result[1] += v[1] * cos_t;
    result[2] += v[2] * cos_t;
    result[0] -= (v[1] * vector[2] - v[2] * vector[1]) * sin_t;
    result[1] -= (v[2] * vector[0] - v[0] * vector[2]) * sin_t;
    result[2] -= (v[0] * vector[1] - v[1] * vector[0]) * sin_t;
    return result;

}


Point Point::PT(Point p_new, Point p_old, bool verbose)//parallel transport of the stored vector, from p_old to p_new
{
	//verbose = true;

	if(verbose)
	{
		printf("in PT\n");
		printf("p old:");p_old.print();
		printf("p new:");p_new.print();
	}

	Point cur(v[0],v[1],v[2]);


	Point axis = p_old^p_new;
	double norm_a = axis.norm();
	double n_old = p_old.norm();
	double n_new = p_new.norm();
	double ratio = norm_a/(n_old*n_new);

	if(ratio < 1e-6)//1 for testing only
	{
		if(p_old*p_new > 0.0)
			return cur;
		else
			return cur*(-1.0);
	}


//	else
//		printf("ratio > 1\n");

	if(verbose)
		printf("norm axis, ratio = %lf, %lf\n",norm_a,ratio);

	axis /= norm_a;


	//Point eo_1 = axis;
	//Point en_1 = axis;

	Point eo_2 = p_old^axis;
	eo_2 /= eo_2.norm();

	Point en_2 = p_new^axis;
	en_2 /= en_2.norm();

	Point res = axis*(axis*cur) + en_2*(eo_2*cur);





	//double theta;
	//if(fabs(ratio-1.0) < 1e-8)
	//	theta = M_PI/2.0;
	//else theta = fabs(atan2(ratio,1.0-pow(ratio,2.0)));//asin(ratio);

	//if(verbose)
	//{
	//	printf("theta = %lf M_PI\n",theta/M_PI);
	//	axis.print();
	//}

	//


	//Point res = rotate(theta,axis);

	//if(verbose)
	//	printf("rotated\n");
	return res;
}

double Point::SymDet()
{
	return v[0] * v[2] - v[1] * v[1];
}

double Point::SymTrace()
{
	return v[0] + v[2];
}

Point Point::SymInverse()
{
	Point r(v[2], (-1.0)*v[1], v[0]);
	return r / SymDet();
}

Point Point::SymInverse_derivative(Point dM)
{
	Point r1(v[2], (-1.0)*v[1], v[0]);
	Point r2(dM[2], (-1.0)*dM[1], dM[0]);
	double det = SymDet();
	r2 /= det;

	double dDet = dM[0] * v[2] + dM[2] * v[0] - 2.0 * v[1] * dM[1];
	r1 *= (dDet/(det*det));

	return r2 - r1;
}

double Point::SymDet_derivative(Point dM)
{
	double dDet = dM[0] * v[2] + dM[2] * v[0] - 2.0 * v[1] * dM[1];
	return dDet;
}

Point Point::SymRoot()
{
	double s = sqrt(fabs(SymDet()));//only valid for positive determinant
	double t = sqrt(SymTrace() + 2.0*s);
	Point r(v[0]+s, v[1], v[2]+s);
	return r / t;
}

double Point::SymNorm2()
{
	return v[0] * v[0] + 2.0*v[1] * v[1] + v[2] * v[2];
}

double Point::SymNorm()
{
	return sqrt(SymNorm2());
}

Point Point::SymInvRoot_derivative(Point dM)
{
	double s = sqrt(fabs(SymDet()));
	double t = sqrt(SymTrace() + 2.0*s);
	double ds_dx = dM[0] * v[2] + dM[2] * v[0] - 2.0 * v[1] * dM[1];
	ds_dx /= (2.0*s);

	double dt_dx = (0.5*dM.SymTrace() + ds_dx) / t;

	Point res(dM[2] + ds_dx, (-1.0)*dM[1], dM[0] + ds_dx);

	res /= (s*t);

	Point res2(v[2] + s, (-1.0)*v[1], v[0] + s);

	res2 *= (dt_dx*s + ds_dx*t) / (s*s*t*t);

	return res - res2;

}

Point Point::SymConjDeriv(int coord) //derivative of | operator with respect to middle term
{
	if (coord < 0 || coord > 2)
	{
		printf("SymConjDeriv only takes coordintes between 0 and 2\n");
		exit(1);
	}

	if (coord == 0)
		return Point(v[0] * v[0], v[0] * v[1], v[1] * v[1]);

	if (coord == 1)
		return Point(2.0*v[0] * v[1], v[0] * v[2] + v[1] * v[1], 2.0*v[1] * v[2]);

//	if (coord == 2)
		return Point(v[1] * v[1], v[1] * v[2], v[2] * v[2]);


	/*Point E1(p2[0] * p2[0], 2.0*p2[0] * p2[1], p2[1] * p2[1]);
	Point E2(p2[0] * p2[1], p2[0] * p2[2] + p2[1] * p2[1], p2[1] * p2[2]);
	Point E3(p2[1] * p2[1], 2.0*p2[1] * p2[2], p2[2] * p2[2]);*/
}


double Point::angle(Point p)
{
	double res = acos(((double)(p[0]*v[0]) + (double)(p[1]*v[1]) + (double)(p[2]*v[2]))/(this->norm()*p.norm()));

	if(res < -1e-8)
		res += 2*M_PI;
	if(res > M_PI)
		res = 2*M_PI - res;
	return res;
}

double Point::theta()
{
	double x = v[0], y = v[1], z = v[2];
	double r = sqrt(x*x + y*y + z*z);
	double res;
	if(r == 0.0)
		return 0.0;
	
	else
	{
		res = acos(z/r);
		if(res > M_PI)
			res -= M_PI;
	
		return res; 
	}
}

double Point::phi()
{
	double x = v[0], y = v[1];
	double res = atan2(y,x);

	if(res < 0.0) res += 2.0*M_PI;
	return res;
}


double Point::line_dist(Point p1, Point p2, Point closest, double * t)
{
	Point cur = Point(v[0],v[1],v[2]);
	if(p1 == p2)
		p2 += Point(1.0,1.0,1.0);
	double t_loc = (p2-p1)*(cur-p1)/(p2-p1).norm2();

	closest = p1 + (p2-p1)*t_loc;

	t = &t_loc;
	return (cur-closest).norm();
}

double Point::line_dist(Point p1, Point p2, Point closest)
{
	Point cur = Point(v[0],v[1],v[2]);
	if(p1 == p2)
		p2 += Point(1.0,1.0,1.0);
	double t_loc = (p2-p1)*(cur-p1)/(p2-p1).norm2();

	closest = p1 + (p2-p1)*t_loc;

	return (cur-closest).norm();
}
	

