//FException

#include "FException.h"

#ifdef _DEBUG
    #undef THIS_FILE
static char THIS_FILE[]=__FILE__;
    #define new DEBUG_NEW
#endif

namespace MeshLib{

FException::FException()
{
    msg = (char*)"failed to open file";
}
FException::FException(char *strg)
{
    msg = strg;
};

FException::~FException(){};


void
FException::what()
{
    std::cout << msg << std::endl;
};


}//end of namespace



