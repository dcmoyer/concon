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



#ifndef MESHLIB_HPP
#define MESHLIB_HPP

//No Dep.
#include "Point.h"
#include "FException.h"
#include "TopologyException.h"

//Depend on Point.h
#include "Trait.h"
#include "Vector4.h"
#include "BarryCentric.h"
#include "Chain.h"
//Depends on Chain.h
#include "BitMap.h"


//Circular.
//Yes. Circular.
//
//These are circular because they're really all describing the same thing;
//if file length were not an issue, they would be put in the same file.
//Technically, Solid could be broken from this cycle by breaking up the
//iterators file. Solid encapsulates all of the above classes to represent
//a mesh.
#include "Vertex.h"
#include "Edge.h"
#include "HalfEdge.h"
#include "Face.h"
//Dependent only through iterators header.
#include "Solid.h"
//Header Only Circular
#include "iterators.h"

//Header Only
#include "AVLTree.h"
#include "TopologyException.h"
#include "string_token_iterator.h"
#include "TopologyException.h"
#include "FaceNormalTrait.h"
#include "vertexPlaneTrait.h"
#include "SList.h"
#include "AVLTree.h"
#include "DList.h"
#include "Chain.h"

#endif

