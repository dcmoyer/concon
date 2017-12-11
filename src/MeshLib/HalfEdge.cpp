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



#include "HalfEdge.h"

using namespace MeshLib;


HalfEdge * HalfEdge::ccw_rotate_about_target()
{
    HalfEdge * he_dual = he_sym();
    if ( he_dual == NULL ) return NULL;

    return he_dual->he_prev();
};

HalfEdge * HalfEdge::clw_rotate_about_target()
{
    HalfEdge * he = he_next()->he_sym();
    return he;
};

HalfEdge * HalfEdge::ccw_rotate_about_source()
{

    HalfEdge * he = he_prev()->he_sym();
    return he;
};

HalfEdge * HalfEdge::clw_rotate_about_source()
{
    HalfEdge * he = he_sym();
    if ( he == NULL ) return NULL;
    return he->he_next();
};
