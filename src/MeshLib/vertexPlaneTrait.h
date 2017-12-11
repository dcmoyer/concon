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

