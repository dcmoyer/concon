#include "Vertex.h"

//THIS IS CIRCULAR BY DESIGN
//See MeshLib.h for explaination.
#include "HalfEdge.h"



namespace MeshLib {

HalfEdge *  Vertex::most_ccw_in_halfedge()  
{ 
    if ( !m_boundary )
    {
        return m_halfedge; //current half edge is the most ccw in halfedge 
    }

    HalfEdge * he = m_halfedge->ccw_rotate_about_target();

    while ( he != NULL )
    {
        m_halfedge = he;
        he = m_halfedge->ccw_rotate_about_target();
    }

    return m_halfedge;
};

HalfEdge *  Vertex::most_clw_in_halfedge()  
{ 
    if ( !m_boundary )
    {
        return most_ccw_in_halfedge()->ccw_rotate_about_target(); 
    }

    HalfEdge * he = m_halfedge->clw_rotate_about_target();

    while ( he != NULL )
    {
        m_halfedge = he;
        he = m_halfedge->clw_rotate_about_target();
    }

    return m_halfedge;
};

HalfEdge *  Vertex::most_ccw_out_halfedge()  
{ 
    if ( !m_boundary )
    {
        return most_ccw_in_halfedge()->he_sym(); 
    }

    HalfEdge * he = m_halfedge->he_next();
    HalfEdge * ne = he->ccw_rotate_about_source();

    while ( ne != NULL )
    {
        he = ne;
        ne = he->ccw_rotate_about_source();
    }

    return he;
};

HalfEdge *  Vertex::most_clw_out_halfedge()  
{ 
    if ( !m_boundary )
    {
        return most_ccw_out_halfedge()->ccw_rotate_about_source(); 
    }

    HalfEdge * he = m_halfedge->he_next();
    HalfEdge * ne = he->clw_rotate_about_source();

    while ( ne != NULL )
    {
        he = ne;
        ne = he->clw_rotate_about_source();
    }

    return he;
};

double Vertex::theta()
{
	double x = m_point[0], y = m_point[1], z = m_point[2];
	double r = sqrt(x*x + y*y + z*z);
	if(r == 0.0)
		return 0.0;
	
	else return acos(z/r);
}

double Vertex::phi()
{
	double x = m_point[0], y = m_point[1];
	double res = atan2(y,x);

	if(res < 0.0) res += 2.0*M_PI;
	return res;
}

void Vertex::color_vertex(double r, double g, double b)
{

	char *rgb = new char[128];       
	std::string key( "rgb" );
	std::string s = Trait::getTraitValue( string(), key );
//	if ( s.length() > 0 )
    {
		sprintf(rgb,"%g %g %g", r,g,b);
        std::string st_rgb = std::string(rgb);
        Trait::updateTraitString(string(),key, st_rgb);
    }
}

void Vertex::set_scalar_attribute(double att, char * key_c)
{

	char *rgb = new char[128];       
	std::string key(key_c);
	std::string s = Trait::getTraitValue( string(), key );
	
	sprintf(rgb,"%g", att);
    std::string st_rgb = std::string(rgb);
    Trait::updateTraitString(string(),key, st_rgb);

	delete [] rgb;
}

bool Vertex::get_scalar_attribute(float *att, char *key_c)
{	
	bool found = true;

	std::string key_s(key_c);
	
	std::string s = Trait::getTraitValue( string(), key_s );
	if( s.length() > 0 )
		sscanf( s.c_str(), "%g ", &att[0] );			
	else
		found = false;

	return found;
}

bool Vertex::get_vector_attribute(Point *att, char *key_c)
{
	bool found = true;
	float a, b, c;

	std::string key_s(key_c);

	std::string s = Trait::getTraitValue(string(), key_s);
	if (s.length() > 0)
	{
		sscanf(s.c_str(), "%g %g %g", &a, &b, &c);
		att[0][0] = a;
		att[0][1] = b;
		att[0][2] = c;
	}
	else
		found = false;

	return found;
}



/*
Edge *Vertex::edge_with(Vertex *v)
{
    for(VertexEdgeIterator eviter(this); !eviter.end(); ++eviter)
    {
        Edge *e = *eviter;
        if(e->other_vertex(this) == v)
            return e;
    }

    return NULL;
}
*/

std::ostream & operator << ( std::ostream & co, const Vertex & v)
{
    co << "Vertex " << v.id() << std::endl;
    return co;
};


}//end of meshlib namespace


