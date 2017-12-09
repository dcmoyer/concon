#ifndef _VERTEXPLANETRAIT_H_
#define _VERTEXPLANETRAIT_H_

#include <string>
#include <iterator>
#include "Trait.h"
//#include "iterators.h"
#include "AVLTree.h"
//#include "Solid.h"
#include "BarryCentric.h"
//#include "int_alias.h"

namespace MeshLib
{

	class int_alias
	{
	public:
		int_alias(){};
		int_alias(int num){id = num;};
		~int_alias(){id = 0;};
	public:
		bool operator< (const int_alias & a) const
        {
            return a.id >= id; 
        };
		bool operator== (const int_alias & a) const
        {
            return a.id == id; 
        };
	public:
		int id;
	};


    class vertexPlaneTrait : public Trait
    {
    public:
        vertexPlaneTrait();

        vertexPlaneTrait(int * list);

        ~vertexPlaneTrait();

        void set(int  *list);

		void reset();

        void clear();


        void read() ;
        void write() ;

    public:
        int *facelist, theFace;
		bool active, found;
		double s,t;
		Barry b;
		Point loc;
		AVL::Tree<int_alias> face_tree;
//	int size;

    };

	

} //namespace 
#endif

