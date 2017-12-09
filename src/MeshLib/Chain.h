// Chain.h: interface for the Chain class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _MESHLIB_CHAIN_H_
#define _MESHLIB_CHAIN_H_
//#include "Solid.h"
#include <cstdlib>
#include "Point.h"

namespace MeshLib
{
    static int chain_dx[8] = {1, 1, 0, -1, -1, -1, 0, 1};
    static int chain_dy[8] = {0, -1, -1, -1, 0, 1, 1, 1};

    class Chain  
    {
    public:
        bool circle(int x, int y);
        bool closed();
        Chain();
        virtual ~Chain(){
          if ( m_directions != NULL )
            delete [] m_directions;
          //I hate whomever wrote this lib
          if( chain_dx == chain_dy ){
            //do nothing
            int x = chain_dx[0];
            x++;
          }
        };
        void set_start_x(int x) { m_start_x = x;};
        void set_start_y(int y) { m_start_y = y;};
        int get_start_x() { return m_start_x;};
        int get_start_y() { return m_start_y;};
        void initialize(int l) { m_length = l; m_directions = new unsigned char[l];};
        void copy_directions(int l, unsigned char * dir) {m_length = l; m_directions = dir;};
        unsigned char get_direction(int i) { if ( i<0 || i>=m_length ) return -1;
            else return m_directions[i];};
        void set_direction(int i, unsigned char dir) { if ( i<0 || i>= m_length ) return;
            else m_directions[i] = dir;};
        void set_end_x(int x) { m_end_x = x;};
        void set_end_y(int y) { m_end_y = y;};
        int get_end_x() { return m_end_x;};
        int get_end_y() { return m_end_y;};
        int get_length() { return m_length;};

    private:
        int m_start_x;
        int m_start_y;
        int m_end_x;
        int m_end_y;
        int m_length;
        unsigned char * m_directions;

    };
};
#endif // !defined(AFX_CHAIN_H__298771C1_E516_46F8_A8E4_187412316EFE__INCLUDED_)
