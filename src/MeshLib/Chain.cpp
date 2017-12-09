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
