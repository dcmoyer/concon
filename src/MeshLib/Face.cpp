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

#include "iterators.h"
#include "Face.h"
#include "HalfEdge.h"

using namespace MeshLib;


std::ostream & operator << ( std::ostream & co, const Face & f)
{
    co << "Face " << f.id() << std::endl;
    return co;
};

bool Face::include_vertex(Vertex *v)
{
    HalfEdge * he = m_halfedge;
    if ( he->target () == v || he->source () == v || he->he_next ()->target () == v )
        return true;
    return false;

}

int Face::on_sphere(Point p, double radius){
    HalfEdge * he = m_halfedge;
    for ( int i=0; i<3; i++ )
    {
        if ( he->target ()->string ().size()>0 && he->target ()->string ().find("valid=(0)") != std::string::npos )
            return 0;
    }
    int count = 0;
    for ( int i=0; i<3; i++ )
    {
        if ( he->edge ()->on_sphere(p,radius) )
            count ++;
        he = he->he_next ();
    }
    return count;
}

Edge * Face::conjunction(Face * face)
{
    HalfEdge * he = m_halfedge;
    for ( int i=0; i<3; i++ )
    {
        if ( face->include_edge (he->edge ()) )
            return he->edge ();
        he = he->he_next ();
    }
    return NULL;
}

Point Face::barycenter()
{
    HalfEdge * he = m_halfedge;
    Point p1 = he->source ()->point ();
	Point p2 = he->target ()->point ();
    Point p3 = he->he_next ()->target ()->point ();
    return ((p1 + p2 + p3)/3.0);
}


double Face::area()
{
    HalfEdge * he1 = m_halfedge;
    HalfEdge * he2 = he1->he_next ();
    Point p1 = he1->target ()->point () - he1->source ()->point ();
    Point p2 = he2->target ()->point () - he2->source ()->point ();
    return(p1^p2).norm ()/2.0;
}

void Face::set_UV_attribute(Point p1, Point p2, Point p3, char * key_c)
{

	char *rgb = new char[128];       
	std::string key(key_c);
	std::string s = Trait::getTraitValue( string(), key );
	
	sprintf(rgb,"%g %g %g %g %g %g", p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
    std::string st_rgb = std::string(rgb);
    Trait::updateTraitString(string(),key, st_rgb);

	delete [] rgb;
}

Point Face::set_UV_barrycentric_coords(Barry b)
{
	Point res;
	Point U, V;

	std::string key = std::string("UV");
	std::string  UV = MeshLib::Trait::getTraitValue(string(), key); 

//	if(UV == NULL)
//		return res;
	if(UV.length() == 0)
		return res;

	sscanf(UV.c_str(),"%lf %lf %lf %lf %lf %lf",&U[0],&V[0],&U[1],&V[1],&U[2],&V[2]);
	
	res[0] = U*b;
	res[1] = V*b;

	return res;
}

int Face::in_sphere(Point p, double radius)
{
    HalfEdge * he = m_halfedge;
    int count = 0;
    for ( int i=0; i<3; i++ )
    {
        Vertex * v = he->vertex ();
        if ( v->in_sphere (p,radius) )
            count++;
        he = he->he_next ();
    }
    return count;
}


Point Face::norm()
{
    HalfEdge * he = m_halfedge;
    Point p1 = he->target ()->point () - he->source ()->point ();
    Point p2 = he->he_next ()->target ()->point () - he->target ()->point ();
    Point n = p1 ^ p2;
    n /= n.norm ();
    return n;
}

//added by kun on 01/25/2005
//include_edge(Edge * e)
bool Face::include_edge(Edge *e)
{
    HalfEdge * he = m_halfedge;
    if ( he->edge () == e || he->he_next ()->edge () == e || he->he_prev ()->edge () == e )
        return true;
    return false;
}

bool Face::segment_cross(Point & p1, Point & p2) // returns true if the segment between the two points crosses the face
{
	Point o,v1,v2, l_p1(p1[0],p1[1],p1[2]), l_p2(p2[0],p2[1],p2[2]);
	HalfEdge * he = m_halfedge;
	in_triangle = false;

	o = he->source()->point();
	v1 = he->target ()->point ();
	v2 = he->he_next()->target ()->point ();

	Point normal = p2 - p1;

	double theta = atan2(normal[0],normal[2]);
    double phi = atan2(normal[1],(normal[0]*sin(theta) + normal[2]*cos(theta)));
	x_rotate(theta,phi, o);
	x_rotate(theta,phi, v1);
	x_rotate(theta,phi, v2);
	x_rotate(theta,phi, l_p1);
	x_rotate(theta,phi, l_p2);

	double x = x_triangle_point(o,v1,v2,l_p1[1],l_p1[2]); //if(in_triangle)printf("planes cross\n");
	if((in_triangle || aligned) && ((x <= l_p1[0] && x >= l_p2[0])||(x >= l_p1[0] && x <= l_p2[0])))
		return true;
	else 
		return false;
}

bool Face::segment_cross2(Point & p1, Point & p2, Point & inter)
{
//	printf("in seg cross 2\n");
	HalfEdge * he = m_halfedge;
	Point o = he->source()->point();
	Point v1 = he->target ()->point ();
	Point v2 = he->he_next()->target ()->point ();
	Point l1 = v1 - o;
	Point l2 = v2 - o;
	Point n = l1^l2;

	n/=n.norm();

	double denom = ((p2-p1)*n);
	if(fabs(denom)/(p2-p1).norm() < 1e-5){
    //std::cout << denom << std::endl;
    //std::cout << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
    //std::cout << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
    //std::cout << n[0] << " " << n[1] << " " << n[2] << std::endl;
    //std::cout << (p2 - p1)[0] << " " << (p2 - p1)[1] << " " << (p2 - p1)[2] << std::endl;
    //std::cout << "Points too close to orthogonal." << std::endl;
		return false;
  }

	double t = (o - p1)*n /denom;

  //NOTES:
  //this is a bug fix from 171114. Who knows how long this has been here.
  //It's technically still working without it (...it gives the correct point)
  //but without this, we could end up with intersect points OUTSIDE p1->p2,
  //which is bad.
  if(t > 1.0){
    return false;
  }

	inter = p1 + (p2 - p1)*t;

//	double A_sum = ((o-inter)^(v1-inter)).norm() + ((v1-inter)^(v2-inter)).norm() + ((v2-inter)^(o-inter)).norm();
//	double A_trig = (l1^l2).norm();

	double angle1 = l1.angle(inter-o);
	double angle2 = l2.angle(inter-o);
	double angle3 = l1.angle(l2);	

	if(fabs(angle1+angle2-angle3) > 1e-6){
		return false;
  }

	angle1 = (o-v1).angle(inter-v1);
	angle2 = (v2-v1).angle(inter-v1);
	angle3 = (v2-v1).angle(o-v1);

	if(fabs(angle1+angle2-angle3) > 1e-6){
		return false;
  }

	angle1 = (o-v2).angle(inter-v2);
	angle2 = (v1-v2).angle(inter-v2);
	angle3 = (v1-v2).angle(o-v2);	

	if(fabs(angle1+angle2-angle3) > 1e-6){
		return false;
  }

	return true;
}


bool Face::segment_cross(Point & p1, Point & p2, double * coords) // returns true if the segment between the two points crosses the face
																 // different from above: if the segment lies in the triangle's plane, it is 
																 // classified as outside
{
	Point o,v1,v2, l_p1(p1[0],p1[1],p1[2]), l_p2(p2[0],p2[1],p2[2]);
	HalfEdge * he = m_halfedge;
//	coords = new double[2];
//	in_triangle = false;

	o = he->source()->point();
	v1 = he->target ()->point ();
	v2 = he->he_next()->target ()->point ();

	/*Point *ap = new Point[3];
	int ind = 0;
	Face *f = he->face();

	 for ( FaceVertexIterator fviter(f); !fviter.end(); ++fviter )
     {
           Solid::tVertex v = *fviter;
           ap[ind] = v->point();
           ind++;
     }

	 o = ap[0];
	 v1 = ap[1];
	 v2 = ap[2];*/

	Point normal = p2 - p1;

	double theta = atan2(normal[0],normal[2]);
    double phi = atan2(normal[1],(normal[0]*sin(theta) + normal[2]*cos(theta)));
	x_rotate(theta,phi, o);
	x_rotate(theta,phi, v1);
	x_rotate(theta,phi, v2);
	x_rotate(theta,phi, l_p1);
	x_rotate(theta,phi, l_p2);

	double x = x_triangle_point(o,v1,v2,l_p1[1],l_p1[2]); //if(in_triangle)printf("planes cross\n");

	coords[0] = m_s;
	coords[1] = m_t;

	if((in_triangle /*|| aligned*/) && ((x <= l_p1[0] && x >= l_p2[0])||(x >= l_p1[0] && x <= l_p2[0])))
		return true;
	else 
		return false;
}

void Face::x_rotate(double theta, double phi, Point & o)
{
	double temp;

	temp = o[0]*cos(theta) - o[2]*sin(theta);
    o[2] =o[0]*sin(theta) + o[2]*cos(theta);
    o[0] = temp;

	temp = o[1]*cos(phi) - o[2]*sin(phi);
    o[2] = o[1]*sin(phi) + o[2]*cos(phi);
    o[1] = temp;

	temp = o[0];
	o[0] = o[2];
	o[2] = temp;
}


double Face::x_triangle_point(Point o, Point v1, Point v2, double y_cur, double z_cur)
{
    Point w1, w2, w3;
    double s = -1.0, t;
    double /*z_cur, y_cur, */y_const, z_const, x1, x2;

    aligned = false;

    y_const = y_cur - o[1];
    z_const = z_cur - o[2];

    w1 = v1 - o;
    w2 = v2 - o;
    w3 = v1 - v2;

    if ( fabs(w1[2]) > FLT_EPSILON /*!= 0.0*/ )
    {
        //x_const -= z_const*w1[1]/w1[2];

        if ( fabs(w2[1] - w2[2]*w1[1]/w1[2]) > FLT_EPSILON /*!=0.0*/ )
        {
            t = (y_const - z_const*w1[1]/w1[2])/(w2[1] - w2[2]*w1[1]/w1[2]);
            s = (z_const - t*w2[2]) /w1[2];
        }
        else
            aligned = true;

    }
    else if ( fabs(w2[2]) > FLT_EPSILON /*!= 0.0*/ )
    {
        t = z_const/w2[2];

        if ( fabs(w1[1]) > FLT_EPSILON /*!= 0.0*/ )
            s = (y_const - w2[1]*t)/w1[1];

        else
            aligned = true;
    }


    else
        aligned = true;

    if ( aligned == true )
    {
        //printf("x aligned\n");
        if ( fabs(w1[0]*w2[2] - w1[2]*w2[0]) > FLT_EPSILON/*!=*/ ) // if the triangle projects a triangle onto the x-z plane
        {
            if ( ( z_cur <= v1[2] && z_cur >= o[2]) || ( z_cur <= v1[2] && z_cur >= o[2]) )
            {
                if ( fabs(w1[2]) > FLT_EPSILON /*!= 0.0*/ )
                {
                    x1 = o[0] + (z_cur -o[2])*w1[0]/w1[2];

                    if ( (( z_cur <= v2[2] && z_cur >= o[2]) || ( z_cur <= v2[2] && z_cur >= o[2]))&&(fabs(w2[2]) > FLT_EPSILON /*!= 0.0*/) )
                        x2 = o[0] + (z_cur -o[2])*w2[0]/w2[2];

                    else
                        x2 = w2[0] + (z_cur - w2[2])*w3[0]/w3[2];
                }

                else
                {
                    x1 = v1[0];
                    x2 = o[0];
                }
            }
            else
            {
                x1 = w2[0] + (z_cur - w2[2])*w3[0]/w3[2];
                x2 = o[0] + (z_cur -o[2])*w2[0]/w2[2];
            }
        }

        else // if the triangle projects a triangle onto the x-y plane
        {
            if ( ( y_cur <= v1[1] && y_cur >= o[1]) || ( y_cur <= v1[1] && y_cur >= o[1]) )
            {
                if ( fabs(w1[1]) > FLT_EPSILON /*!= 0.0*/ )
                {
                    x1 = o[0] + (y_cur -o[1])*w1[0]/w1[1];

                    if ( (( y_cur <= v2[1] && y_cur >= o[1]) || ( y_cur <= v2[1] && y_cur >= o[1]))&&(fabs(w2[1]) > FLT_EPSILON /*!= 0.0*/) )
                        x2 = o[0] + (y_cur -o[1])*w2[0]/w2[1];

                    else
                        x2 = w2[0] + (y_cur - w2[1])*w3[0]/w3[1];
                }

                else
                {
                    x1 = v1[0];
                    x2 = o[0];
                }
            }
            else
            {
                x1 = w2[0] + (y_cur - w2[1])*w3[0]/w3[1];
                x2 = o[0] + (y_cur -o[1])*w2[0]/w2[1];
            }
        }

        return(x1 + x2)/2.0;
    }
    if ( s >= 0.0 && t >= 0.0 && s + t <= 1.0 + FLT_EPSILON ) in_triangle = true;
    else in_triangle = false;

	m_s = s;
	m_t = t;

    return o[0] + s*w1[0] + t*w2[0];

}

Point Face::project(Point p)
{
	HalfEdge * he = m_halfedge;
	Point o = he->source()->point();
	Point v1 = he->target()->point();
	Point v2 = he->he_next()->target()->point();

	Point p1 = v1 - o;
	Point p2 = v2 - o;

	p1 /= p1.norm();
	p2 -= p1*(p1*p2);
	p2 /= p2.norm();

	return Point(p1*p,p2*p,0.0);
}

Point Face::apply_projection(Point p)
{
	HalfEdge * he = m_halfedge;
	Point o = he->source()->point();
	Point v1 = he->target()->point();
	Point v2 = he->he_next()->target()->point();

	Point p1 = v1 - o;
	Point p2 = v2 - o;

//	printf("p1\t");p1.print();
//	printf("p2\t");p2.print();

	p1 /= p1.norm();
	p2 -= p1*(p1*p2);
	p2 /= p2.norm();

//	printf("p1\t");p1.print();
//	printf("p2\t");p2.print();

	return p1*p(0) + p2*p(1);
}

Point Face::apply_coords(double *coords)
{
	HalfEdge * he = m_halfedge;
	Point o = he->source()->point();
	Point v1 = he->target()->point();
	Point v2 = he->he_next()->target()->point();

	Point p1 = v1 - o;
	Point p2 = v2 - o;
	
	return o + p1*coords[0] + p2*coords[1];
}

Point Face::find_coords(Point pp, double *coords)
{
	return find_coords(pp,norm() + pp, coords);
}

Point Face::find_coords_from_origin(Point pp, double *coords)
{
	Point zero;
	return find_coords(zero, pp, coords);
}

Point Face::find_coords(Point pp1, Point pp2, double *coords)
{
	HalfEdge * he = m_halfedge;
	Point o = he->source()->point();
	Point v1 = he->target()->point();
	Point v2 = he->he_next()->target()->point();

	Point p1 = v1 - o;
	Point p2 = v2 - o;
	Point n = p1^p2;
	double dot_n = (pp2-pp1)*n;
	Point res;

	if(dot_n > 0.0)
	{
		double t_norm = (o-pp1)*n/(dot_n);

		res = pp1 + pp2*t_norm;

		double dot12 = p1*p2;
		double norm1 = p1.norm2();
		double norm2 = p2.norm2();
		double dotp1 = res*p1;
		double dotp2 = res*p2;
	
		double det = norm1*norm2-dot12*dot12;
		double s = dotp1*norm2 - dot12*dotp2;
		double t = dotp2*norm1 - dot12*dotp1;
		t /= det;
		s /= det;

		coords[0] = s;
		coords[1] = t;
	}
	else
		coords[0] = coords[1] = -1.0;
	
	return res;
}

Point Face::find_coords(Point pp1, Point pp2, Point coords)//third coord is the segment parameterization
{
	HalfEdge * he = m_halfedge;
	Point o = he->source()->point();
	Point v1 = he->target()->point();
	Point v2 = he->he_next()->target()->point();

	Point p1 = v1 - o;
	Point p2 = v2 - o;
	Point n = p1^p2;
	double dot_n = (pp2-pp1)*n;
	Point res;
	double t_norm = -1.0;

	if(dot_n > 0.0)
	{
		t_norm = (o-pp1)*n/(dot_n);

		res = pp1 + pp2*t_norm;

		double dot12 = p1*p2;
		double norm1 = p1.norm2();
		double norm2 = p2.norm2();
		double dotp1 = res*p1;
		double dotp2 = res*p2;
	
		double det = norm1*norm2-dot12*dot12;
		double s = dotp1*norm2 - dot12*dotp2;
		double t = dotp2*norm1 - dot12*dotp1;
		t /= det;
		s /= det;

		coords[0] = s;
		coords[1] = t;
		coords[2] = t_norm;
	}
	else
		coords[0] = coords[1] = coords[2] = -1.0;
	
	return res;
}

bool Face::segment_cross3(Point &p1, Point &p2, double *coords)
{
	Point cds_pt;

	find_coords(p1, p2, cds_pt);
	coords[0] = cds_pt[0];
	coords[1] = cds_pt[1];

	if(cds_pt[0] >= 0.0 && cds_pt[1] >= 0.0 && 
		(cds_pt[0] + cds_pt[1]) <= 1.0 && 
		cds_pt[2] >= 0.0 && cds_pt[2] <= 1.0)
		return true;

	return false;
}


Point Face::find_projection(Point p)
{
	HalfEdge * he = m_halfedge;
	Point o = he->source()->point();
	Point tmp = p - o;
	return this->project(tmp);
}

Barry Face::find_barrycentric_coords(Point p)
{
	Point *w = new Point[3];
	HalfEdge * he = m_halfedge;
	w[0] = he->source()->point();
	w[1] = he->target()->point();
	w[2] = he->he_next()->target()->point();

	Barry res(w,p);
	delete [] w;
	return res;
}

Barry Face::find_barrycentric_coords2(Point pp)
{	
	Point *w = new Point[3];
	HalfEdge * he = m_halfedge;
	w[0] = he->source()->point();
	w[1] = he->target()->point();
	w[2] = he->he_next()->target()->point();
  std::cout << w[1][0] << w[1][1] << w[1][2] << std::endl;

	Point n_unsc = (w[0] - w[1])^(w[2] - w[1]);
  
	Point p = pp*(w[0]*n_unsc)/(pp*(n_unsc));

	Barry res(w,p);
	delete [] w;
	return res;
}

//FIND POINT FROM BARRY
Point Face::set_barrycentric_coords(Barry b)
{
	Point *w = new Point[3];
	HalfEdge * he = m_halfedge;
	w[0] = he->source()->point();
	w[1] = he->target()->point();
	w[2] = he->he_next()->target()->point();

	Point res = w[0]*b(0) + w[1]*b(1) + w[2]*b(2);
	delete [] w;
	return res;
}

float Face::set_barrycentric_coords(Barry b, float * att, bool nearest_neighbor)
{
	int *w = new int[3];
	HalfEdge * he = m_halfedge;
	w[0] = he->source()->id2()-1;
	w[1] = he->target()->id2()-1;
	w[2] = he->he_next()->target()->id2()-1;

	float res = att[w[0]]*b(0) + att[w[1]]*b(1) + att[w[2]]*b(2);

	if (nearest_neighbor)
	{
		int near_ind = 0;
		if (b(0) < b(1)) near_ind = 1;
		if (b(0) < b(2) && b(1) < b(2)) near_ind = 2;
		res = att[w[near_ind]];
	}

	delete [] w;
	return res;
}

double Face::set_barrycentric_coords(Barry b, double * att)
{
	int *w = new int[3];
	HalfEdge * he = m_halfedge;
	w[0] = he->source()->id2()-1;
	w[1] = he->target()->id2()-1;
	w[2] = he->he_next()->target()->id2()-1;

	double res = att[w[0]]*b(0) + att[w[1]]*b(1) + att[w[2]]*b(2);
	delete [] w;
	return res;
}

Point Face::set_barrycentric_coords(Barry b, Point * att)
{
	int *w = new int[3];
	HalfEdge * he = m_halfedge;
	w[0] = he->source()->id2()-1;
	w[1] = he->target()->id2()-1;
	w[2] = he->he_next()->target()->id2()-1;

	Point res = att[w[0]]*b(0) + att[w[1]]*b(1) + att[w[2]]*b(2);
	delete [] w;
	return res;
}

Point Face::get_world_coords_from_barry(Point b){

	Point *w = new Point[3];
	HalfEdge * he = m_halfedge;
	w[0] = he->source()->point();
	w[1] = he->target()->point();
	w[2] = he->he_next()->target()->point();

	Point res = w[0]*b[0] + w[1]*b[1] + w[2]*b[2];
	delete [] w;
  return res;
}


