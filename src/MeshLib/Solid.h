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

#ifndef _MESHLIB_SOLID_H_
#define _MESHLIB_SOLID_H_

#define MAX_LINE 2048

#include <math.h>
#include <fstream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <float.h>
#include <string.h>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include "Edge.h"
#include "Face.h"
#include "HalfEdge.h"
#include "Vertex.h"
#include "Point.h"
#include "Trait.h"
#include "AVLTree.h"
#include "SList.h"
#include "Vector4.h"
#include "TopologyException.h"
#include "FException.h"
#include "FaceNormalTrait.h"
#include "vertexPlaneTrait.h"
#include "DList.h"

//#include "BitMap.h"
#include "Chain.h"

//#include "common_routines.h"
//#include "pca2.h"
//#include "DList.h"
//#include "BitMap.h"

//TODO:Remove this
int round_to_int( float a );
bool point_match(MeshLib::Point &p1, MeshLib::Point &p2);
double cotsum(MeshLib::Point &p1, MeshLib::Point &p2, MeshLib::Point &p3);

namespace MeshLib{

class Vertex;
class HalfEdge;
class Edge;
class Face;

class SolidHalfEdgeIterator;
class SolidVertexIterator;
class SolidEdgeIterator;
class SolidFaceIterator;
class Configure;

//!  Solid.
/*!
  This class define solid(mesh) structure.
*/

class Solid 
{

public:

  // has to be gong someday!

  // pointer to Vertices, Halfedges, Edges, Face and Solid
  //!  Vertex pointer.
  typedef Vertex   * tVertex;
  //!  Halfedge pointer.
  typedef HalfEdge * tHalfEdge;
  //!  Edge pointer.
  typedef Edge     * tEdge;
  //!  Face pointer.
  typedef Face     * tFace;

  //!  Solid pointer.
  typedef Solid    * tSolid;

  //constructor and destructor
  //!  Constructor.
  Solid(){};
  //!  Destructor.
  virtual ~Solid();


  void MetricDistortion(Solid & solid, float * corners, float * edges, float * vert_angles, float * vert_dist);

  //copy operator
  //!  Solid copy operator.
  /*!
    \param solid an Solid which is the copying source..

  */
  void copy( Solid & solid );

  //
  //!  Add more geometry to the mesh.
  /*!
    \param solid an Solid which is the copying source..
    \param uv the UV coordinates
  */
  void add( Solid & solid );

  void add_arrow(Point direction, Point origin, double scale, double r = -1.0, double g = -1.0, double b = -1.0);

  void add_point_sphere(Point origin, double radius, double r, double g, double b);

  //adds texture rgb trait based on UV map and a bitmap file provided. 
  //Must have UV map, so currently only for alias obj files (11/11/13)
  //void add_texture(char * fn_bmp);

  //reads .m files without manifold assumptions and removes non-manifold faces
  void read_NM( const char * input, Point *test = NULL, int n_test_pts = 0 );

  //reads ALIAS ascii obj format, fixes non-manifold topology    
  void readOBJ( const char * input);

  //reads ascii ply format
  void readPLY( const char * input);
  
  //reads ascii ply format, fixes non-manifold topology
  void readPLY_NM(const char * input);

  //reads off format
  void readOFF( const char * input);

  //labels faces from a "polygon soup" as non-manifold
  void label_non_manifold_faces(Point *V, Point *F, bool * use_face, int v_counter,int f_counter);

  void label_non_manifold_faces_test(Point *V, Point *F, bool * use_face, int v_counter, int f_counter, Point * test, int n_test_pts);

  //file io
  //!  Read solid from istream.
  /*!
    \param is an istream to get information from.
  */
  void read(  std::istream & is );
  //!  Write solid to isstream.
  /*!
    \param os an ostream to put information to.
  */
  void write( std::ostream & os );

  //reset id2() values for vertex-array indexing
  void reset_id2s();

  //!  Read solid from input file.
  /*!
    \param input the char pointer to indicate input file name.
    Now a generic.
  */
  void read(const char * input, std::string filetype="");

  void read_minc_obj(const std::string& filename);

  void read_old(const char * input, bool reset_ids = false, int v_id_start = 0, int f_id_start = 0);
  void read_fs_binary(const char * input, bool reset_ids = false, int v_id_start = 0, int f_id_start = 0);

  //rotates mesh and/or conformal position on the sphere by the input Euler angles
  void EulerRotate(double alpha, double beta, double gamma, int flag) ;

  //aligns main axes of the shape with the coordinates z and x, returns Euler angles used for transformation
  Point align_main_axes();

  //sets uniform rgb trait for all vertices
  void color_solid(double r, double g, double b);

  //colors vertices according to attribute array
  void color_solid(float * atr, double min, double max, bool reverse = false);

  //colors vertices according to attribute array
  void color_solid(float * atr, bool reverse = false);

  //sets new attribute array over vertices
  void set_scalar_attribute(float * atr, char* key);

  //gets attribute array over vertices based on known key. default value set where key is not found
  void get_scalar_attribute(float * atr, char* key, double default_val);

  //gets attribute 3D vector array over vertices based on known key. default value set where key is not found
  void get_vector_attribute(Point * atr, char* key, Point default_val);

  //returns the standard location of a spherical point at given index and bandwidth
  Point point_on_sphere(int index, int bandwidth);

  //tiles the sphere with standard faces
  int sphere_faces(int bandwidth, int* m_face_x, int* m_face_y, int* m_face_z);

  //warps sphere by tangent vector field, assumes geodesic length = arcsin(||v||)
  void warp_sphere_arcsin(Point *v);

  //creates a mesh of a regular sphere with given resolution
  void generate_regular_sphere(int bandwidth, int cur_verts = 0, int cur_faces = 0);

  //sets pole vertices of a regular shape to be average of the pole boundaries
  void reset_regular_poles(int bw);

  //closes pole boundaries of a regularly sampled spherical shape
  void close_std_shape(int bw, int cur_verts = 0, int cur_faces = 0);

  //closes a given boundary by adding a new vertex at the average position of boundary verts
  void close_boundary(int *list, int & max_id, int & max_id_f, bool set_fill);
  
  //closes all boundaries by adding a new vertex at the average position of boundary verts
  void close_boundaries(bool set_fill = false);

  //removes all connected components except the largest
  void keep_largest_connected_piece(bool mark_verts_only = false);

  //removes connected components that consist of only two faces & no boundaries
  void removeDegenerateComponents();


  //returns center of mass of the shape
  Point center_mass();

  Point Min_corner();
  Point Max_corner();

  void reverse_vertex_order();

  //computes exact volume, assuming surface is closed
  double volume();

  //!  Write solid to output file.
  /*!
    \param output the char pointor to indicate output file name.
  */
  void write( const char * output);

  //write selectively, only writing verts/faces with mask value above threshold
  void write( const char * output, float * mask, float thresh );
  
  //!  Update trait string with new normal vectors.
      
  void UpdateNormalTrait();


  //number of vertices, faces, edges
  //!  Get # of vertices.
  /*!
   \return an int which is the # of vertices.
  */
  int  numVertices();
  //!  Get # of edges.
  /*!
    \return an int which is the # of edges.
  */
  int  numEdges();
  //!  Get # of face.
  /*!
    \return an int which is the # of faces.
  */
  int  numFaces();

  double surfaceArea();

  //returns max(max_x-min_x,max_y-min_y,max_z-min_z)
  double max_eucledian_extent(); //returns max(max_x-min_x,max_y-min_y,max_z-min_z)

  double average_edge_length();

  double resize_by_area(double desired_surface_area);

  //is boundary
  //!  Is boundary for vertex.
  /*!      
    \param v a tVertex to be test it is boundary or not.
    \return the boolean result that the vertex is boundary or not.
  */
  bool    isBoundary( tVertex  v );

  //!  Is boundary for edge.
  /*!      
    \param e a tEdge to be test it is boundary or not.
    \return the boolean result that the edge is boundary or not.
  */
  bool    isBoundary( tEdge    e );

  //!  Is boundary for halfedge.
  /*!      
    \param he a tHalfEdge to be test it is boundary or not.
    \return the boolean result that the halfedge is boundary or not.
  */
  bool    isBoundary( tHalfEdge  he );

  //acess vertex - id
  //!  Acess vertex from id.
  /*!      
    \param id an int which is the ID of the vertex which needs to be access.
    \return the tVertex which is the vertex with its ID = given id.
  */
  tVertex idVertex( int id );

  //!  Get vertex ID from given vertex.
  /*!      
    \param v a tVertex that needs to find its id.
    \return a int which is the id of the given vertex.
  */
  int     vertexId( tVertex  v );

  //access face - id
  //!  Acess face from id.
  /*!      
    \param id an int which is the ID of the face which needs to be access.
    \return the tFace which is the face with its ID = given id.
  */
  tFace   idFace( int id );

  //!  Get face ID from given face.
  /*!      
    \param v a tFace that needs to find its id.
    \return a int which is the id of the given face.
  */
  int     faceId( tFace  f );

  //access edge - edge key, vertex
  //!  Acess edge from ids.
  /*!      
    \param id0 an int which is the source vertex ID of edge.
    \param id1 an int which is the target vertex ID of edge.
    \return the tEdge which is the edge with its IDs = given ids.
  */
  tEdge   idEdge( int id0, int id1 );
  //!  Acess edge from vertices.
  /*!      
    \param id0 a tVertex which is the source vertex of the edge.
    \param id1 a tVertex which is the target vertex of the edge.
    \return the tEdge which is the edge with its source/target vertex = given verticess.
  */
  tEdge   vertexEdge( tVertex v0, tVertex v1 );

  //access halfedge - halfedge key, vertex
  //!  Acess halfedge from ids.
  /*!      
    \param id0 an int which is the source vertex ID of the halfedge.
    \param id1 an int which is the target vertex ID of the halfedge.
    \return the tHalfEdge which is the halfedge with its IDs = given ids.
  */
  tHalfEdge   idHalfedge( int id0, int id1 );

  //!  Acess halfedge from vertices.
  /*!      
    \param id0 a tVertex which is the source vertex of the halfedge.
    \param id1 a tVertex which is the target vertex of the halfedge.
    \return the tHalfEdge which is the halfedge with its source/target vertex = given verticess.
  */
  tHalfEdge   vertexHalfedge( tVertex v0, tVertex v1 );

  //!  Acess halfedge from corner.
  /*!      
    \param v a tVertex which is on the oppsite side of the halfedge.
    \param f a tFace which is the face thich is associated with the halfedge.
    \return the tHalfEdge which is the halfedge from its face and corner information.
  */
  tHalfEdge   corner( tVertex v, tFace f);

  //halfedge->face
  //!  Acess face from halfedge.
  /*!      
    \param he a tHalfEdge.
    \return the tFace which is associated with the given halfedge.
  */
  tFace   halfedgeFace( tHalfEdge he );

  //halfedge->vertex
  //!  Acess vertex from halfedge.
  /*!      
    \param he a tHalfEdge.
    \return the tVertex which is associated with the given halfedge (target vertex).
  */
  tVertex halfedgeVertex( tHalfEdge he );

  //edge->vertex
  //!  Acess vertex from edge.
  /*!      
    \param e a tEdge.
    \return the tVertex which is associated with the given edge (smaller vertex id).
    \return For an unboundary edge, condition edgeVertex1(edge).id() < edgeVertex(edge).id()
    \return is strictly enforced.
    \return For a boundary edge, no such condition exists.
  */
  tVertex edgeVertex1( tEdge  e );

  //!  Acess vertex from edge.
  /*!      
    \param e a tEdge.
    \return the tVertex which is associated with the given edge (larger vertex id).
  */
  tVertex edgeVertex2( tEdge  e );

  //edge->face
  //!  Acess face from edge.
  /*!      
    \param e a tEdge.
    \return the tFace which is the 1st associated with halfedge which belongs to the given edge.
  */
  tFace edgeFace1( tEdge  e );

  //!  Acess face from edge.
  /*!      
    \param e a tEdge.
    \return the tFace which is the snd associated with halfedge which belongs to the given edge.
  */
  tFace edgeFace2( tEdge  e );

  //Euler operations
  //!  Access most clw out halfedge from a vertex.
  /*!      
    \param v a tVertex which is the base center for clw access.
    \return the tHalfEdge which is moust clw out halfedge from the given vertex.
  */
  tHalfEdge vertexMostClwOutHalfEdge( tVertex  v );

  //!  Access next ccw out halfedge from a halfedge.
  /*!      
    \param he a tHalfEdge which is the base center for ccw access.
    \return the tHalfEdge which is next ccw out halfedge from the given halfedge.
  */
  tHalfEdge vertexNextCcwOutHalfEdge( tHalfEdge  he );

  //!  Access most ccw out halfedge from a vertex.
  /*!      
    \param v a tVertex which is the base center for clw access.
    \return the tHalfEdge which is most ccw out halfedge from the given vertex.
  */
  tHalfEdge vertexMostCcwOutHalfEdge( tVertex  v );

  //!  Access next clw out halfedge from a halfedge.
  /*!
    \param he a tHalfEdge which is the base center for clw access.
    \return the tHalfEdge which is next clw out halfedge from the given halfedge.
  */
  tHalfEdge vertexNextClwOutHalfEdge( tHalfEdge  he );

  //!  Access most clw in halfedge from a vertex.
  /*!      
    \param v a tVertex which is the base center for clw access.
    \return the tHalfEdge which is most clw in halfedge from the given vertex.
  */
  tHalfEdge vertexMostClwInHalfEdge( tVertex  v );

  //!  Access next ccw in halfedge from a halfedge.
  /*!
    \param he a tHalfEdge which is the base center for ccw access.
    \return the tHalfEdge which is next ccw in halfedge from the given halfedge.
  */
  tHalfEdge vertexNextCcwInHalfEdge( tHalfEdge  he );

  //!  Access most ccw in halfedge from a vertex.
  /*!      
    \param v a tVertex which is the base center for ccw access.
    \return the tHalfEdge which is most ccw in halfedge from the given vertex.
  */
  tHalfEdge vertexMostCcwInHalfEdge( tVertex  v );

  //!  Access next ccw in halfedge from a halfedge.
  /*!
    \param he a tHalfEdge which is the base center for clw access.
    \return the tHalfEdge which is next ccw in halfedge from the given halfedge.
  */
  tHalfEdge vertexNextClwInHalfEdge( tHalfEdge  he );

  //!  Access most clw in halfedge from a face.
  /*!      
    \param f a tFace which is the base center for clw access.
    \return the tHalfEdge which is most clw in halfedge from the given face.
  */
  tHalfEdge faceMostClwHalfEdge( tFace  f );

  //!  Access most ccw in halfedge from a face.
  /*!      
    \param f a tFace which is the base center for ccw access.
    \return the tHalfEdge which is most ccw in halfedge from the given face.
  */
  tHalfEdge faceMostCcwHalfEdge( tFace  f );

  //!  Access next ccw in halfedge from a face.
  /*!      
    \param f a tFace which is the base center for ccw access.
    \return the tHalfEdge which is next ccw in halfedge from the given face.
  */
  tHalfEdge faceNextCcwHalfEdge( tHalfEdge  he );

  //!  Access next ccw in halfedge from a halfedge.
  /*!      
    \param f a tHalfEdge which is the base center for ccw access.
    \return the tHalfEdge which is next ccw in halfedge from the given halfedge.
  */
  tHalfEdge faceNextClwHalfEdge( tHalfEdge  he );

  //!  Get edge length.
  /*!      
    \param e a tEdge to be get length.
    \return a double which is the length of this edge.
  */
  double edgeLength( tEdge e );

  //Dynamic part

  //!  Split an edge.
  /*!      
    \param e a tEdge to be split.
    \return a tVertex which is the new vertex after edge split.
  */
  tVertex   edgeSplit( tEdge e );
  tVertex   trianglesubdivision( tFace  f , tEdge e1 , tEdge e2 , tEdge e3 );

  //!split a face to add a new point

  void add_point_face(Point & p, Face * f);

  //!split an edge to add a new point
  void add_point_edge(Point & p, Edge *edge);

  //!add point to the mesh that is closest to the space point
  //returns a 4-vector with entries as follows:
  //res(projection type, element id, coord1, coord2)    or
  //res(projection type, element id1, element id1, coord)  
  //projectiontype: 1.0 - vertex, 2.0 - edge, 3.0 - face
  //element id: one entry for face, vertex, 2 entries for edge
  //coord1, coord2 - for Face::apply_coords
  //coord - for Edge::project_segment
  Vector4 project_point(Point p, bool add_to_mesh = true);

  //!use result of project points
  Point apply_coords(Vector4 coords, int * min_id = NULL);

  void segment_regular_by_curve(int bw, Point * curve, int N, float * mask, bool dilate_by_edge_length = false);

  //! 
  Vertex * closest_vert_to_point(Point p);

  Vertex * closest_vert_to_point(Point p, int init_id, int s_rad);

  void Hausdorff(Solid * mesh, int s_rad, int fix_iters, int * ids);

  //!  Label boundary edges.
  /*!
      \param void
      \return void
  */
  void labelBoundaryEdges();


  //! removeDanglingVertices.
  /*!
      \param void
      \return void
  */
  void removeDanglingVertices();

  void UpdateNormals(void);

  void _FaceNormal(MeshLib::Solid *pMesh);
  void _VertexNormal(MeshLib::Solid *pMesh);

  void readBYU(  std::istream & is );
  void writeBYU( std::ostream & os );

  protected:

  //std::map<EdgeKey, tEdge, ltedgekey >              m_edges;

  //!  AVL tree contains edges.
  AVL::Tree<Edge>                                   m_edges;
  //!  AVL tree contains vertices.
  AVL::Tree<Vertex>                                 m_verts;
  //!  AVL tree contains faces.
  AVL::Tree<Face>                                   m_faces;

  //  bool * m_use_face;
  //  Point * m_F, *m_V; // storing initial "polygon soup"

  protected:

  friend class SolidVertexIterator;
  friend class SolidEdgeIterator;
  friend class SolidFaceIterator;
  friend class SolidHalfEdgeIterator;
  friend class DSolid;
  friend class WSolid;
  friend class TDSolid;
  friend class Configure;
  friend class SolidDelegate;

  //!  Create a new vertex.
  /*!      
    \param id a int which is the id of the new vertex (default = 0).
    \return a tVertex which is the new vertex.
  */
  tVertex   createVertex(   int id = 0 );
  //!  Destroy a vertex
  /*!
    \param vertex a vertex which is to be destroyed
  */
  void destroyVertex(Vertex * vertex);

  //! Destroy an Edge
  /*!
    \param edge an edge which is to be destroyed
  */
  void destroyEdge(Edge * edge);

  //! Destroy a HalfEdge
  /*!
    \param he a halfedge which is to be destroyed
  */
  void destroyHalfEdge(HalfEdge * he);

  //! Destroy a Face
  /*!
    \param face a face which is to be destroyed
  */
  void destroyFace(Face * face);

  //!  Create a new edge.
  /*!      
    \param start a int which is the id of the souce vertex of the edge.
    \param end a int which is the id of the target vertex of the edge.
    \return a tEdge which is the new edge.
  */
  tEdge     createEdge( int start, int end );

  //!  Create a new edge.
  /*!      
    \param id a int which is the id of the new face (default = 0).
    \return a tFace which is the new face.
  */
  tFace     createFace( int id = 0 );

  //!  Create a new face.
  /*!      
    \param v a int pointer which contains 3 vertices position information.
    \param id a int which is the id of the new face.
    \return a tFace which is the new face.
  */
  tFace     createFace( int * v, int id ); //create a triangle

  //! Collapse Edge
  /*!
    \param edge an edge to be collapsed
  */
  void    collapseEdge(Edge * edge);

  //! Collapse Edge according to vertex
  /*!
    \param edge an edge to be collapsed
    \param vertex the vertex on the edge which should be saved
  */
  void    collapseEdgeVertex(Edge * edge, Vertex * vertex);
  bool    collapsable(HalfEdge * halfedge);

};

}//end of namespace MeshLib

#endif //_MESHLIB_SOLID_H_ defined

