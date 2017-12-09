
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

