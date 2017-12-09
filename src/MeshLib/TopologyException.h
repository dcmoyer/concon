//TopologyException

#ifndef _TOPOLOGY_EXCEPTION_H_
#define _TOPOLOGY_EXCEPTION_H_

#include <iostream>


namespace MeshLib
{

    class TopologyException 
    {
    public:
        TopologyException();
        TopologyException(char * strg);

        ~TopologyException();

        void what();

    private:

        char * msg;
    };


}//end of namespace MeshLib
#endif //TopologyException

