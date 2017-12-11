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

// chain.cpp: implementation of the Chain class.
//
//////////////////////////////////////////////////////////////////////

#include "Chain.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
namespace MeshLib{
Chain::Chain()
{

}


bool Chain::closed()
{
    return m_start_x == m_end_x && m_start_y == m_end_y;

}

bool Chain::circle(int x, int y)
{
    double chain_x = (double)m_start_x;
    double chain_y = (double)m_start_y;
    Point start_p = Point(chain_x,chain_y,0) - Point(x,y,0);
    if ( start_p.norm () == 0 ) return true;
    start_p /= start_p.norm ();
    double sum_theta = 0;
    for ( int i=0; i<m_length; i++ )
    {
        chain_x += (double)chain_dx[m_directions[i]];
        chain_y += (double)chain_dy[m_directions[i]];
        Point p = Point(chain_x, chain_y,0) - Point(x,y,0);
        if ( p.norm () == 0 ) return true;
        p/= p.norm ();
        Point thp = start_p ^ p;
        if ( thp[2]>0 )
            sum_theta += thp.norm ();
        else
            sum_theta -= thp.norm ();
        start_p = p;
    }
    if ( sum_theta > acos(-1.0) ) return true;
    else return false;
}
}
