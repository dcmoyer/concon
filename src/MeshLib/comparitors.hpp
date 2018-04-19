//Copyright 2017 Yalin Wang, Boris Gutman, and Daniel Moyer
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

#ifndef MESHLIB_COMPARITORS
#define MESHLIB_COMPARITORS

#include "MeshLib/Solid.h"
#include "MeshLib/iterators.h"

#include <vector>

namespace MeshLib {

//
//  is_equal_in_fv
//    This method checks the equivalence of two meshes.
//    It only checks:
//      -vertex position and order
//      -face sets and face order
//      -edge order (implicitly)
//
//  WARNING: This explicitly checks order. Thus, if I permute face or vertex
//  labels, and according change the structure to match, the mesh will be
//  reported as not equivalent. This is intended!
//
//  (many downstream representations rely on Face or Vertex correspondence)
//
bool is_equal_in_fv(
  MeshLib::Solid& left, MeshLib::Solid& right,
  bool verbose=false);

}//end of namespace MeshLib

#endif



