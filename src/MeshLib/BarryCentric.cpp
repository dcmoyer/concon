
#include "BarryCentric.h"

using namespace MeshLib;

Barry::Barry( Point * w, Point & p )
{
    Point n = (w[1] - w[0])^(w[2]-w[0]);
    double area = n.norm();
    if ( area < 1e-14 ) throw std::domain_error("area must be >1e-14");
    n/=area;

    for ( int i = 0; i < 3; i ++ )
    {
        Point a = (p-w[(i+1)%3])^(p-w[(i+2)%3]);
        v[i] = a*n/area;
    }

};

