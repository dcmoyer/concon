//TopologyException

#include "TopologyException.h"

#ifdef _DEBUG
    #undef THIS_FILE
static char THIS_FILE[]=__FILE__;
    #define new DEBUG_NEW
#endif

namespace MeshLib{

TopologyException::TopologyException()
{
    msg = (char*)"failed to satisfy the manifold constraints";
}
TopologyException::TopologyException(char *strg)
{
    msg = strg;
};

TopologyException::~TopologyException(){};


void
TopologyException::what()
{
    std::cout << msg << std::endl;
};

}
