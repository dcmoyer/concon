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

//#include "stdafx.h"

#include "Edge.h"


//THIS IS CIRCULAR BY DESIGN
//See MeshLib.h for explaination.
#include "HalfEdge.h"
#include "Vertex.h"
#include "Face.h"

using namespace MeshLib;

bool Edge::coface(Vertex * v)
{
    HalfEdge * he0 = m_halfedge[0];
    if ( he0->face()->include_vertex(v) )
        return true;
    if ( m_halfedge[1] != NULL && m_halfedge[1]->face()->include_vertex(v) )
        return true;
    return false;
}

bool Edge::on_sphere (Point p, double radius){
    Vertex * v1 = m_halfedge[0]->source ();
    Vertex * v2 = m_halfedge[0]->target ();
    Point p1 = v1->point ();
    Point p2 = v2->point ();
    if ( ((p1-p).norm ()-radius)*((p2-p).norm ()-radius)<=0 )
        return true;
    return false;
}

void Edge::crosspoint (Point p, double radius, double & b0, double & b1)
{
    Point p0 = m_halfedge[0]->source ()->point ();
    Point p1 = m_halfedge[0]->target ()->point ();
    if ( (p-p0).norm ()<1e-15 )
    {
        b0 = radius;
        b1 = radius;
        return;
    }
    double cost = (p1-p0) * (p-p0)/(p1-p0).norm ()/(p-p0).norm ();
    double a = (p-p0).norm ();
    if ( radius*radius-a*a+a*a*cost*cost<0 )
    {
        b0 = -1;
        b1 = -1;
        return;
    }
    double t = sqrt(fabs(a*a*cost*cost-a*a+radius*radius));
    b0 = a * cost + t;
    b1 = a * cost - t;
}

bool Edge::coface(Edge * e)
{
    HalfEdge * he1 = m_halfedge[0];
    HalfEdge * he2 = m_halfedge[1];
    HalfEdge * hee1 = e->halfedge (0);
    HalfEdge * hee2 = e->halfedge (1);
    if ( he1->face ()==hee1->face () )
        return true;
    if ( he2 != NULL )
    {
        if ( he2->face () == hee1->face () ) return true;
    }
    if ( hee2 != NULL )
    {
        if ( hee2->face () == he1->face () ) return true;
    }
    if ( he2 != NULL && hee2 != NULL )
    {
        if ( he2->face () == hee2->face () ) return true;
    }
    return false;
};

Vertex * Edge::other_vertex(Vertex *v)
{
    HalfEdge * he = m_halfedge[0];
    if ( he->target () == v ) return he->source ();
    else if ( he->source () == v ) return he->target ();
    return NULL;

}
double Edge::length()
{
    HalfEdge * he = m_halfedge[0];
    Point p1 = he->target ()->point ();
    Point p2 = he->source ()->point ();
    return(p1-p2).norm ();
}

Vertex * Edge::conjunction(Edge *e)
{
    if ( (m_halfedge[0]->source()->id() == e->vertex(0)) || (m_halfedge[0]->source()->id() == e->vertex (1)) )
        return m_halfedge[0]->source();
    else if ( m_halfedge[0]->target()->id() == e->vertex(0) || m_halfedge[0]->target()->id() == e->vertex (1) )
        return m_halfedge[0]->target();
    else return NULL;

}

void Edge::get_vertices(Vertex *&v1, Vertex *&v2)
{
    v1 = m_halfedge[0]->source();
    assert(v1);
    v2 = other_vertex(v1);
    assert(v2);
}

bool Edge::include_point(Point & p)
{
	HalfEdge * he = m_halfedge[0];
	Point p1 = he->source()->point();
	Point p2 = he->target()->point();

	if( (p1 - p).norm() + (p2-p).norm() - (p1-p2).norm() <= 0.0)
		return true;

	return false;
}

Point Edge::project_segment(double t)
{
	HalfEdge * he = m_halfedge[0];
	Point p1 = he->source()->point();
	Point p2 = he->target()->point();

	Point curve_dir = p2 - p1;
	return p1 + curve_dir*t;
}
	
