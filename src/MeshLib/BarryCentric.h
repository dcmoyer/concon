#ifndef _MESHLIB_BARRY_H_
#define _MESHLIB_BARRY_H_

#include <string>
#include <iterator>
#include <stdexcept>
#include "Point.h"

namespace MeshLib
{

    class Barry : public Point
    {
    public:
        Barry( Point * v, Point & p );
		Barry(){Point();};
    };


}; //namespace 
#endif

