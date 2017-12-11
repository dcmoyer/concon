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


#ifndef _MESHLIB_FACE_H_
#define _MESHLIB_FACE_H_

#include <iostream>
#include <string>
#include <assert.h>
#include <math.h>
#include <float.h>
#include "Edge.h"
#include "BarryCentric.h"

namespace MeshLib{


    class HalfEdge;
    class Trait;

    class Face
    {
    public:

        Face(){ m_halfedge = NULL; m_trait=NULL; m_s = m_t = 0.0;};
        ~Face(){};

        Point norm();
        HalfEdge    * & halfedge() { return m_halfedge;};
        Trait       * & trait()    { return m_trait;};
        int           & id()       { return m_id;};
        const int       id() const { return m_id;};
    //  int           & id2()      { return m_id2;};
    //  const int       id2()const { return m_id2;};
        double          area();
    Point barycenter();
    bool segment_cross(Point & p1, Point & p2);
    bool segment_cross(Point & p1, Point & p2, double * coords);
    bool segment_cross2(Point & p1, Point & p2, Point & inter);
    bool segment_cross3(Point & p1, Point & p2, double * coords);

    Point apply_coords(double *coords);//applies result of Face::segment_cross, or find_coords
    Point find_coords(Point pp1, double *coords);
    Point find_coords_from_origin(Point pp1, double *coords);
    Point find_coords(Point pp1, Point pp2, double *coords);
    Point find_coords(Point pp1, Point pp2, Point coords);
    Point project(Point p); //returns weights for projection: proj = (n1)*res(0) + (n2)*res(1)
    //o = halfedge()->source()->point();
    //v1 = halfedge()->target()->point();
    //v2 = halfedge()->he_next()->target()->point();
    //n1 = (v1-o)/||v1=0||
    //n2 = (v2-o) - n1*((v2-o)*n1); n2 /= ||n2||
    Point apply_projection(Point p); // applies result of Face::project
    Point find_projection(Point p);

    Barry find_barrycentric_coords(Point p);
    Barry find_barrycentric_coords2(Point p);
    Point find_point_from_barry(Barry b);

    //FIND POINT FROM BARRY
    Point set_barrycentric_coords(Barry p);
    float set_barrycentric_coords(Barry b, float * att, bool nearest_neighbor = false);
    double set_barrycentric_coords(Barry b, double * att);
    Point set_barrycentric_coords(Barry b, Point * att);
    Point set_UV_barrycentric_coords(Barry b);
    Point get_world_coords_from_barry(Point b);

    bool include_edge(Edge * e);
    bool include_vertex(Vertex * v);
    int in_sphere(Point p, double radius);
    int on_sphere(Point p, double radius);
    Edge * conjunction(Face * face);

    void set_UV_attribute(Point p1, Point p2, Point p3, char * key_c);

        std::string &   string()   { return m_string;};

        bool operator== (const Face & f) const
        {
            return m_id == f.id();
        };

        bool operator< (const Face & f) const
        {
            return f.id() >= m_id; 
        };
  private:
    double x_triangle_point(Point o, Point v1, Point v2, double y_cur, double z_cur);
    void x_rotate(double theta, double phi, Point & o);


    private:
        int           m_id;
        //int           m_id2;
        HalfEdge    * m_halfedge;
        std::string   m_string;   //string
        Trait       * m_trait;
    bool aligned, in_triangle;
    double m_s, m_t;

    };

    std::ostream & operator << ( std::ostream & co, const Face & f);

}//name space MeshLib

#endif //_MESHLIB_FACE_H_ defined

