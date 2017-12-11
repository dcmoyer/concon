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

//#include "SparseMatrix2.h"
#include "Solid.h"

#include "iterators.h"
//#include "BitMap.h"
//#ifdef WIN32
//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif
//#endif
//
//#ifdef WIN32
//
//#include "crtdbg.h"
//
//#ifdef _DEBUG
//#ifndef DEBUG_NEW
//#define DEBUG_NEW   new( _NORMAL_BLOCK, __FILE__, __LINE__)
//#endif
//#else
//#ifndef DEBUG_NEW
//#define DEBUG_NEW
//#endif
//#endif
//
//#else
//
//#ifndef DEBUG_NEW
//#define DEBUG_NEW
//#endif
//
//#endif          // WIN32

//TODO:Remove this
int round_to_int( float a )
{
    int int_a = (int)floor(a + 0.5);
    return int_a;
}

bool point_match(MeshLib::Point &p1, MeshLib::Point &p2)
{
	int i1 = -1, i2 = -1, i3 = -1;

	for(int i = 0; i < 3; i++)
		if(p1[0] == p2[i])
			i1 = i;
	if(i1 < 0)
		return false;

	for(int i = 0; i < 3; i++)
		if(i != i1 && p1[1] == p2[i])
			i2 = i;
	if(i2 < 0)
		return false;

	for(int i = 0; i < 3; i++)
		if(i != i1 && i != i2 && p1[2] == p2[i])
			i3 = i;
	if(i3 < 0)
		return false;

	return true;
}

double cotsum(MeshLib::Point &p1, MeshLib::Point &p2, MeshLib::Point &p3)
{
	MeshLib::Point A1 = p1-p3;
	MeshLib::Point A2 = p2-p3;
	
	double epsilon = 1e-6;
	double res = fabs(A1*A2)/((A1^A2).norm() + epsilon);

	A1 = p1-p2;
	A2 = p3-p2;
	res += fabs(A1*A2)/((A1^A2).norm() + epsilon);

	A1 = p2-p1;
	A2 = p3-p1;
	res += fabs(A1*A2)/((A1^A2).norm() + epsilon);

	return res;
}

namespace MeshLib{


//!  Solid.
/*!
  This class define solid(mesh) structure.
*/
Solid::tFace Solid::createFace( int * v , int id )
{
	/*if(id == 25683)*/
	//	printf("in create face %i\n",id);

    tFace f = createFace( id );  
	/*if(id == 25683)*/
	//	printf("created\n");
    //create halfedges
    tHalfEdge hes[3];

    int i;
    for ( i = 0; i < 3; i ++ )
    {
//		printf("i = %i\n",i);
        hes[i] = new HalfEdge;

        Vertex * vert =  idVertex( v[i] );
//		printf("vert id-d\n");
//		if(vert == NULL) printf("no such vert, %i\n",v[i]);
        hes[i]->vertex() = vert;
        vert->halfedge() = hes[i];
    }

//	printf("1st loop\n");

    //linking to each other
    for ( i = 0; i < 3; i ++ )
    {
        hes[i]->he_next() = hes[(i+1)%3];
        hes[i]->he_prev() = hes[(i+2)%3];
    }

//	printf("2nd loop\n");
	

    //linking to face
    for ( i = 0; i < 3; i ++ )
    {
        hes[i]->face()   = f;
        f->halfedge()    = hes[i];
    }

	/*if(id == 25683)*/
	//	printf("done creating halfedges\n");

    //connecting with edge
    //this part is modified: fix the constraint " the target vertex id is 
    //greater than the source vertex id of halfedge(0) with interior edge "
    for ( i = 0; i < 3; i ++ )
    {
        tEdge e = createEdge( v[i], v[(i+2)%3] );
	//	printf("edge %i done\n",i);
	/*	if(id == 25683)*/
	//		printf("edge %i created\n", i);
        if ( e->halfedge(0)  == NULL )
        {
			/*if(id == 25683)*/
//				printf("halfedge(0) of %i NULL\n", i);
            e->halfedge(0) = hes[i];
            hes[i]->edge() = e;
        }
        else
        {
			/*if(id == 25683)*/
//				printf("halfedge(0) of %i exists\n", i);
            if ( e->halfedge(1) != NULL )
			{
				printf("edge %i %i already has two halfedges\n", v[i], v[(i+2)%3] );
			//	destroyHalfEdge(hes[i]);
			//	return NULL;
				throw TopologyException();
			}

            e->halfedge(1) = hes[i];
            hes[i]->edge() = e;

            if ( e->halfedge(0)->target()->id() < e->halfedge(0)->source()->id() )
            {
				/*if(id == 25683)*/
//					printf("last if at %i\n", i);

                HalfEdge * temp = e->halfedge(0);
                e->halfedge(0 ) = e->halfedge(1);
                e->halfedge(1 ) = temp;
            }

        }
    }

	/*if(id == 25683)*/
	//	printf("done creating edges, returning\n");

    return f;
};
//this function should always be the last to execute
//the vertex should always be isolated before destroy
void Solid::destroyVertex(Vertex * vertex)
{
    assert(vertex->halfedge () == NULL);
    assert(vertex->trait () == NULL);

    m_verts.remove (vertex);
    delete vertex;
}

//note this detroy will generate isolated vertex
void Solid::destroyEdge (Edge * edge)
{
    assert(edge->halfedge (0) == NULL && edge->halfedge (1) == NULL);
    assert(edge->trait () == NULL);
    m_edges.remove (edge);
    delete edge;
}

//note this function also can not manipulate the generation of isolated vertex
//this function assume no vertex has link to the halfedge to be destroyed
void Solid::destroyHalfEdge(HalfEdge * he)
{
    assert(he->target ()->halfedge () != he);
    assert(he->trait () == NULL);
    HalfEdge * hes = he->he_sym ();
    if ( hes != NULL )
    {
        hes->he_sym () = NULL;
        if ( hes == hes->edge ()->halfedge (1) )
        {
            Edge * e = hes->edge ();
            e->halfedge (0) = hes;
            e->halfedge (1) = NULL;
        }
    }
    else
    {
        Edge * e = he->edge ();
        e->halfedge (0) = NULL;
        e->halfedge (1) = NULL;
        destroyEdge(e);
    }
    HalfEdge * hep = he->he_prev ();
    HalfEdge * hen = he->he_next ();
    if ( hep!= NULL ) hep->he_next () = NULL;
    if ( hen!= NULL ) hen->he_prev () = NULL;
    if ( he->face ()->halfedge () == he )
        he->face ()->halfedge () = he->he_next ();
    delete he;
}

//this will happen only after all the halfedge on the face is deleted
void Solid::destroyFace (Face * face)
{
    assert(face->trait () == NULL);
    HalfEdge * he = face->halfedge ();
    for ( int i=0; i<3; i++ )
    {
        HalfEdge * hen = he->he_next ();
        assert(he->target ()->halfedge () == NULL || he->target ()->halfedge ()->face () != face);
        destroyHalfEdge(he);
        he = hen;
    }
//	if(face != NULL)
		assert(face->halfedge () == NULL);
    m_faces.remove (face);
    delete face;
}


void Solid::collapseEdge (Edge * edge)
{
    collapseEdgeVertex(edge,edge->halfedge (0)->source () );
}

bool Solid::collapsable (HalfEdge * halfedge)
{
    Edge * edge = halfedge->edge ();
    Vertex * source = halfedge->source ();
    Vertex * target = halfedge->target ();
    Face * fL = halfedge->face ();
    Face * fR = halfedge->he_sym ()!= NULL?halfedge->he_sym()->face ():NULL;
    Vertex * vL = halfedge->he_next ()->target ();
    Vertex * vR = halfedge->he_sym ()!= NULL? halfedge->he_sym ()->he_next ()->target ():NULL;
    std::vector <HalfEdge *> thlist;
    for ( VertexInHalfedgeIterator vhiter( this, target );!vhiter.end(); ++ vhiter )
    {
        HalfEdge * he = *vhiter;

        if ( he->face() == fL ) continue;
        if ( he->face() == fR ) continue;

        thlist.push_back ( he );
    }

    for ( VertexInHalfedgeIterator svhiter( this, source );!svhiter.end(); ++ svhiter )
    {
        HalfEdge * he = *svhiter;

        if ( he->face() == fL ) continue;
        if ( he->face() == fR ) continue;

        thlist.push_back ( he );
    }
    if ( thlist.size () == 0 ) return false;

    thlist.clear ();
    for ( VertexInHalfedgeIterator lvhiter( this, vL );!lvhiter.end(); ++ lvhiter )
    {
        HalfEdge * he = *lvhiter;

        if ( he->face() == fL ) continue;
        if ( he->face() == fR ) continue;

        thlist.push_back ( he );
    }
    if ( thlist.size () == 0 ) return false;
    if ( vR != NULL )
    {
        thlist.clear ();
        for ( VertexInHalfedgeIterator rvhiter( this, vR );!rvhiter.end(); ++ rvhiter )
        {
            HalfEdge * he = *rvhiter;

            if ( he->face() == fL ) continue;
            if ( he->face() == fR ) continue;

            thlist.push_back ( he );
        }
        if ( thlist.size () == 0 ) return false;
    }

    if ( source->boundary() && target->boundary() && !edge->boundary() )
        return false;

    List<Vertex> verts;

    for ( VertexVertexIterator vviter( source ); !vviter.end(); ++ vviter )
    {
        Vertex * v = *vviter;
        if ( v == vL ) continue;
        if ( v == vR ) continue;
        verts.Append( v );
    }

    for ( VertexVertexIterator viter(  target ); !viter.end(); ++ viter )
    {
        Vertex * v = *viter;
        if ( verts.contains( v ) ) return false;
    }

    return true;
}

void Solid::collapseEdgeVertex (Edge * edge, Vertex * vertex)
{
    List<Edge> changekey;
    Vertex * vt = edge->other_vertex (vertex);
    for ( VertexEdgeIterator eiter(vt); !eiter.end (); ++eiter )
    {
        Edge * e = *eiter;
        if ( edge->halfedge (0)->face ()->include_edge (e) )
            continue;
        if ( edge->halfedge (1)!= NULL && edge->halfedge (1)->face ()->include_edge (e) )
            continue;
        changekey.Append (e);
    }
    List<Face> del_f;
    del_f.Append (edge->halfedge (0)->face ());
    if ( edge->halfedge (1) != NULL )
        del_f.Append (edge->halfedge (1)->face ());
    List<HalfEdge> retain_he;
    List<HalfEdge> sym_he;
    for ( int i=0; i<2; i++ )
    {
        HalfEdge * he = edge->halfedge (i);
        if ( he!= NULL )
        {
            if ( he->target () == vertex )
            {
                retain_he.Append (he->he_next ()->he_sym ());
                sym_he.Append (he->he_prev ()->he_sym ());
            }
            else
            {
                retain_he.Append (he->he_prev ()->he_sym ());
                sym_he.Append (he->he_next ()->he_sym ());
            }
        }
    }
    std::vector <Vertex * > iso_v;
    for ( ListIterator<Face> fiter(del_f); !fiter.end (); ++fiter )
    {
        Face * face = *fiter;
        HalfEdge * he = face->halfedge ();
        for ( int i=0; i<3; i++ )
        {
            if ( he->target ()->halfedge () == he )
            {
                Vertex * v = he->target ();
                HalfEdge * ohe = NULL;
                for ( VertexInHalfedgeIterator hiter(this, v);!hiter.end (); ++hiter )
                {
                    HalfEdge * heo = *hiter;
                    if ( !del_f.contains (heo->face ()) )
                    {
                        ohe = heo;
                        break;
                    }
                }
                if ( ohe == NULL )
                {
                    v->halfedge () = NULL;
                    iso_v.push_back (v);
                }
                else
                {
                    v->halfedge () = ohe;
                }
            }
            he = he->he_next ();
        }
    }
    for ( ListIterator<Face> fiter(del_f); !fiter.end (); ++fiter )
    {
        Face * f = *fiter;
        destroyFace(f);
    }
    ListIterator<HalfEdge> siter(sym_he);
    for ( ListIterator<HalfEdge> hiter(retain_he); !hiter.end () && !siter.end (); ++hiter, ++siter )
    {
        HalfEdge * he = *hiter;
        HalfEdge * se = *siter;
        if ( he == NULL && se == NULL ) continue;
        if ( he != NULL && se != NULL )
        {
            he->he_sym () = se;
            Edge * e = se->edge ();
            for ( int i=0; i<2; i++ )
            {
                HalfEdge * hes = e->halfedge (i);
                if ( hes!=NULL )
                {
                    if ( hes->target () == vt )
                    {
                        hes->target () = vertex;
                    }
                    else if ( hes->source () == vt )
                        hes->source () = vertex;
                    if ( hes->target () == vertex && vertex->halfedge () == NULL ) vertex->halfedge () = hes;
                    hes->edge () = he->edge ();
                }
            }
            e->halfedge (0) = NULL;
            e->halfedge (1) = NULL;
            destroyEdge(e);
        }
        else if ( he!= NULL && se == NULL )
        {
            he->he_sym () = NULL;
        }
        else
        {
            assert(he == NULL && se != NULL);
            changekey.Append (se->edge ());
        }
    }

    for ( ListIterator<Edge> eiter(changekey); !eiter.end (); ++eiter )
    {
        Edge * e = *eiter;
        m_edges.remove (e);
        if ( e->ekey ().s () == vt->id () )
        {
            e->ekey () = EdgeKey(vertex->id (),e->ekey ().t ());
        }
        else if ( e->ekey ().t () == vt->id () )
            e->ekey () = EdgeKey(e->ekey ().s (),vertex->id ());
        assert(m_edges.find (e) == NULL);
        m_edges.insert (e);
        for ( int i=0; i<2; i++ )
        {
            HalfEdge * he = e->halfedge (i);
            if ( he != NULL )
            {
                if ( he->target () == vt )
                    he->target () = vertex;
                else if ( he->source () == vt )
                    he->source () = vertex;
                if ( he->target () == vertex && vertex->halfedge () == NULL )
                    vertex->halfedge () = he;
            }
        }
    }
    std::vector<Vertex *>::iterator iter = std::find(iso_v.begin (),iso_v.end (),vt);
    if ( iter == iso_v.end () )
    {
        vt->halfedge () = NULL;
        destroyVertex(vt);
    }

    for ( int i=0; i<(int)iso_v.size (); i++ )
    {
        Vertex * v = iso_v[i];
        if ( v!= vertex )
            destroyVertex(v);
    }
}

void Solid::color_solid(double r, double g, double b)
{
//	printf("in color solid\n");
	AVL::TreeIterator<Vertex> viter( m_verts );
    for ( ; !viter.end() ; ++ viter )
    {
        tVertex v = *viter;
		v->color_vertex(r,g,b);
	}
//	printf("solid colored\n");
}

void Solid::EulerRotate(double alpha, double beta, double gamma, int flag) 
{
    float w1,w2,w3;
    char *conformal = new char[128];

    AVL::TreeIterator<Vertex> viter( m_verts );
    for ( ; !viter.end() ; ++ viter )
    {
        tVertex v = *viter;

        std::string key( "conformal" );
        std::string s = Trait::getTraitValue( v->string(), key );
        if ( flag != 2 )
            if ( s.length() > 0 )
            {

                sscanf( s.c_str(), "%g %g %g", &w1, &w2, &w3 );
                Point pc(w1, w2, w3);
                pc.EulerRotate(alpha,beta,gamma);

                sprintf(conformal,"%g %g %g", pc[0],pc[1],pc[2]);
                std::string rotconformal = std::string(conformal);
                Trait::updateTraitString(v->string(),key, rotconformal);
            }

        if ( flag != 1 )
            v->point().EulerRotate(alpha,beta,gamma);
    }

	delete[] conformal;
}

//TODO:Review this
//Point Solid::align_main_axes()
//{
//	pca2 PCA;
//
//	/*for (j = 1; j <= m; j++)
//    {
//    mean[j] = 0.0;
//    for (i = 1; i <= n; i++)
//        {
//        mean[j] += data[i][j];
//        }
//    mean[j] /= (double)n;
//    }*/
//
//	Point *P = new Point[3];
//
//	const int nf = numFaces();
//	double ** data = new double*[numFaces()+1], ** vecs = new double*[4];
//	
//	for(int i = 0; i < numFaces()+1; i++)
//		data[i] = new double[4];
//
//	for(int i = 0; i < 4; i++)
//		vecs[i] = new double[4];
//
//	data[0][0] = data[0][1] = data[0][2] = data[0][3] = 0.0;
//
//	Point cm = center_mass();
//
//	int i = 1;
//	for(SolidFaceIterator fiter(this); !fiter.end(); ++fiter)
//	{
//		Face *f = *fiter;
//		data[i][0] = 0.0;
//
//		Point p = (f->barycenter()-cm)*f->area();
//		data[i][1] = p(0);
//		data[i][2] = p(1);
//		data[i][3] = p(2);
//		i++;
//	}
//
//
//
//	PCA.pca_eigvecs(data, numFaces(), 3, vecs, 'S');
//
//	for (int j = 1; j <= 3; j++) {
//       for (i = 1; i <= 3; i++)  {
//          printf("%12.4f", vecs[j][4-i]); 
//			P[3-i][j-1] = vecs[j][4-i];
//	
//	   
//			}
//          printf("\n");  }
//
//
//	Point x(1.0,0.0,0.0);
//	Point z(0.0,0.0,1.0);
//
//	/*Point p1_rot = P[0].PT(z,P[2]);
//	Point p2_rot = P[1].PT(z,P[2]);
//	Point p3_rot = P[2].PT(z,P[2]);*/
//
//	
//
//	//Solid dirs;
//
//	/*dirs.add_arrow(P[0],cm,5,1.0,0.0,0.0);
//	dirs.add_arrow(P[1],cm,5,0.0,1.0,0.0);
//	dirs.add_arrow(P[2],cm,5,0.0,0.0,1.0);*/
//
//	double gamma1 = 2*M_PI - atan2(P[2][1],P[2][0]);
//	double beta1 = M_PI/2.0 - atan2(P[2][2],sqrt(P[2][1]*P[2][1] + P[2][0]*P[2][0]));
//
//	Point p1_rot = P[0];//.EulerRotate(0.0,beta1,gamma1);
//	Point p2_rot = P[1];//.EulerRotate(0.0,beta1,gamma1);
//	Point p3_rot = P[2];//.EulerRotate(0.0,beta1,gamma1);
//
//	p1_rot.EulerRotate(0.0,beta1,gamma1);
//	p2_rot.EulerRotate(0.0,beta1,gamma1);
//	p3_rot.EulerRotate(0.0,beta1,gamma1);
//
//	/*p2_rot.print();
//	p3_rot.print();*/
//
//	double alpha1 = 2*M_PI - atan2(p2_rot[1],p2_rot[0]);
//
//	EulerRotate(alpha1,beta1,gamma1,0);
//
//	p1_rot = P[0];//.EulerRotate(0.0,beta1,gamma1);
//	p2_rot = P[1];//.EulerRotate(0.0,beta1,gamma1);
//	p3_rot = P[2];//.EulerRotate(0.0,beta1,gamma1);
//
//	p1_rot.EulerRotate(alpha1,beta1,gamma1);
//	p2_rot.EulerRotate(alpha1,beta1,gamma1);
//	p3_rot.EulerRotate(alpha1,beta1,gamma1);
//
//
//	/*for(SolidVertexIterator viter(this); !viter.end(); ++viter)
//	{
//		Vertex *v = *viter;
//		Point p = v->point();
//		p = p.PT(z,P[2]);
//		p = p.PT(x,p2_rot);
//		v->point() = p;
//	}*/
//
//	/*p1_rot = p1_rot.PT(x,P[1]);
//	p2_rot = p2_rot.PT(x,P[1]);*/
//
//	//cm = center_mass();
//	//
//	//dirs.add_arrow(Point(0.0,1.0,0.0),cm,5,1.0,0.0,0.0);
//	//dirs.add_arrow(x,cm,5,0.0,1.0,0.0);
//	//dirs.add_arrow(z,cm,5,0.0,0.0,1.0);
//
//	//dirs.add_arrow(p1_rot,cm,7.5,1.0,0.2,0.0);	
//	//dirs.add_arrow(p2_rot,cm,7.5,0.0,1.0,0.2);
//	//dirs.add_arrow(p3_rot,cm,7.5,0.0,0.2,1.0);
//	//
//	//dirs.write("dirs.m");
//
//
//	for (int j = 0; j <= numFaces(); j++)
//		delete [] data[j];
//	for (int j = 0; j <= 3; j++)
//		delete [] vecs[j];
//	
//
//	delete [] vecs;
//	delete [] data;
//	delete [] P;
//
//	return Point(alpha1,beta1,gamma1);
//}


void Solid::UpdateNormalTrait() 
{
    char *normal = new char[128];

    AVL::TreeIterator<Vertex> viter( m_verts );
    for ( ; !viter.end() ; ++ viter )
    {
        Vertex * v = *viter;

        std::string key( "normal" );
        std::string s = Trait::getTraitValue( v->string(), key );
      
        if ( s.length() > 0 )
        {          
			Point n = v->normal();          

            sprintf(normal,"%g %g %g", n[0],n[1],n[2]);
            std::string newnormal = std::string(normal);
            Trait::updateTraitString(v->string(),key, newnormal);
        }      
    }
}





int Solid::numVertices() 
{
    /*
    
    int sum = 0;

    AVL::TreeIterator<Vertex> viter( m_verts );
    for( ; !viter.end(); ++viter  )
    {
        sum ++;
    }
    return sum;
    
    */
    return m_verts.getSize();
};

int Solid::numEdges() 
{
/*
    int sum = 0;

    AVL::TreeIterator<Edge> eiter( m_edges );
    for( ; !eiter.end(); ++eiter )
    {
        sum ++;
    }
    return sum;
*/

    return m_edges.getSize();
};

int Solid::numFaces() 
{
    /*
    int sum = 0;

    AVL::TreeIterator<Face> fiter( m_faces );
    for( ; !fiter.end(); ++fiter )
    {
        sum ++;
    }
    return sum;
    */

    return m_faces.getSize();
};

double Solid::surfaceArea()
{
    double area = 0.0;

    AVL::TreeIterator<Face> fiter( m_faces );
    for ( ; !fiter.end(); ++fiter )
    {
        tFace f = *fiter;
        area += f->area();
    }
    return area;
}

double Solid::resize_by_area(double desired_surface_area)
{
	double old_area = surfaceArea();
	double ratio = sqrt(desired_surface_area / old_area);
	AVL::TreeIterator<Vertex> viter(m_verts);
	for (; !viter.end(); ++viter)
	{
		Vertex * v = *viter;
		Point p = v->point();
		p *= ratio;
		v->point() = p;
	}

	return ratio;
}


double Solid::average_edge_length()
{
	double length = 0.0;

    AVL::TreeIterator<Edge> eiter( m_edges );
    for ( ; !eiter.end(); ++eiter )
    {
        Edge * e = *eiter;
		length += e->length();
    }
	
    return length/(double)numEdges();
}

double Solid::max_eucledian_extent() //returns max(max_x-min_x,max_y-min_y,max_z-min_z)
{
	Point Pmin, Pmax;
	int i = 0;

	printf("about to do extent loop\n");

	printf("num verts = %i\n",numVertices());

	AVL::TreeIterator<Vertex> viter( m_verts );

	for(; !viter.end(); ++viter)
	{
		tVertex  v = * viter;

		if(i == 0)
		{
			Pmin = v->point();
			Pmax = v->point();
		}

		for(int j = 0; j < 3; j++)
		{
			if(Pmax[j] < v->point()[j])
				Pmax[j] = v->point()[j];
			if(Pmin[j] > v->point()[j])
				Pmin[j] = v->point()[j];
		}

		i++;
	}

	Point dif = Pmax - Pmin;
	printf("thru extent loop\n");
	dif.print();


	return std::max(dif[0],std::max(dif[1],dif[2]));
}


Point Solid::point_on_sphere(int index, int bandwidth)
{
	double x,y,z, dphi, theta, dheight;
	int dsize;

	if(index == 4*bandwidth*bandwidth)
		return Point(0.0,0.0,-1.0);
	if(index == 4*bandwidth*bandwidth+1)
		return Point(0.0,0.0,1.0);

	dsize = 2*bandwidth;
	dphi = (double)(index%dsize);

	dheight = (double)(index - (index%dsize))/dsize;
	theta = (2*dheight+1.0)/(4*bandwidth);

	z=(-1.0)*cos(M_PI*(theta));
	x=sqrt(1-pow(z,2))*cos(2*M_PI*(dphi/dsize));
	y=sqrt(1-pow(z,2))*sin(2*M_PI*(dphi/dsize));

	return Point(x,y,z);
}

void Solid::warp_sphere_arcsin(Point *p)
{
//	printf("in warp sphere arcsin\n");
	//freopen("test_dem.txt","w",stdout);
	for(SolidVertexIterator viter(this); !viter.end(); ++viter)
	{
		Vertex * v = *viter;
		Point loc = v->point();

		//printf("in solid %i %lf\n",v->id(),p[v->id()-1].norm());

		double CS = sqrt(1.0-p[v->id2()-1].norm2());
		v->point() = loc*CS + p[v->id2()-1];	
	}	
}

void Solid::generate_regular_sphere(int bandwidth, int cur_verts, int cur_faces)
{
//	int cur_verts = 0;//m_verts.getSize();
	int counter = 0, id, verts = 4*bandwidth*bandwidth;
	while (counter < verts)
    {
		id = counter + 1;
        Point p;
		p = point_on_sphere(counter,bandwidth);
        tVertex v  = createVertex( cur_verts + id );
        v->point() = p;
        v->id()    = cur_verts + counter + 1;
		v->id2()   = cur_verts + counter + 1;
		counter++;
   }
	
	int size=2*bandwidth;
	int* m_face_x = (int *) malloc(sizeof(int) * (size-1)*size*2);
	int* m_face_y = (int *) malloc(sizeof(int) * (size-1)*size*2);
	int* m_face_z = (int *) malloc(sizeof(int) * (size-1)*size*2);
	int faces = sphere_faces(bandwidth, m_face_x, m_face_y,  m_face_z);

	counter = 0;
	while(counter < faces)
    {	
		id = counter + 1;
		int v[3];
		v[0] = cur_verts + m_face_x[counter];
		v[1] = cur_verts + m_face_y[counter];
		v[2] = cur_verts + m_face_z[counter];
		createFace( v, id + cur_faces);
		counter++;
   }

	labelBoundaryEdges();

	free(m_face_x);
	free(m_face_y);
	free(m_face_z);
}


int Solid::sphere_faces(int bandwidth, int* m_face_x, int* m_face_y, int* m_face_z)

{

	int i,j,k, size, local[6];

	size=2*bandwidth;
	
	k=0;
	i=0;
	while(i<size-1){

		j=0;
		while(j<size/2)
		
		{
			local[0]=size*i + j*2 + 1;
			local[1]=size*i + j*2 + 2;
			local[2]=size*i + j*2 + 3;
			local[3]=local[0] + size ;
			local[4]=local[1] + size ;
			local[5]=local[2] + size ;

			if(j==size/2-1)
			{
				local[2]=local[2]-size;
				local[5]=local[5]-size;
			}
			if(i%2==0)
			{
				m_face_x[4*k]  =  local[0];
				m_face_y[4*k]  =  local[1];
				m_face_z[4*k]  =  local[4];

				m_face_x[4*k+1]=  local[0];
				m_face_y[4*k+1]=  local[4];
				m_face_z[4*k+1]=  local[3];

				m_face_x[4*k+2]=  local[1];
				m_face_y[4*k+2]=  local[2];
				m_face_z[4*k+2]=  local[4];

				m_face_x[4*k+3]=  local[2];
				m_face_y[4*k+3]=  local[5];
				m_face_z[4*k+3]=  local[4];

			}
	
			else

			{
				m_face_x[4*k]  =  local[0];
				m_face_y[4*k]  =  local[1];
				m_face_z[4*k]  =  local[3];

				m_face_x[4*k+1]=  local[1];
				m_face_y[4*k+1]=  local[4];
				m_face_z[4*k+1]=  local[3];

				m_face_x[4*k+2]=  local[1];
				m_face_y[4*k+2]=  local[2];
				m_face_z[4*k+2]=  local[5];

				m_face_x[4*k+3]=  local[1];
				m_face_y[4*k+3]=  local[5];
				m_face_z[4*k+3]=  local[4];
	
			}
			j++;
			k++;
		}
		i++;
	//	printf("on i = %i\n",i);
	}
	return 4*k;
}

Point Solid::center_mass()

{
	int j;
	Point ap[3], n, c;
	double area, total_area = 0;
	AVL::TreeIterator<Face> fiter( m_faces );
	Point mass_center( 0.0 , 0.0 , 0.0 );

	for( ; !fiter.end(); ++ fiter )
	{
		Solid::tFace f = *fiter;
		j = 0;
		for(FaceVertexIterator fviter(f); !fviter.end(); ++fviter)
		{
			Solid::tVertex v = *fviter;
			ap[j] = v->point();
			j++;
		}
		n = (ap[1] - ap[0])^(ap[2]-ap[0]);
		area = n.norm();
		c = (ap[0]+ap[1]+ap[2])/3.0;
		mass_center += c * area;
		total_area += area;
	}
	mass_center /= total_area;

	return mass_center;
}

Point Solid::Min_corner()
{
	SolidVertexIterator viter(this);
	Point p_min;
	int i = 0;
	for(;!viter.end(); ++ viter)
	{
		Vertex * v = * viter;
		if(i == 0) p_min = v->point();
		else
		{
			if(v->point()[0] < p_min[0]) p_min[0] = v->point()[0];
			if(v->point()[1] < p_min[1]) p_min[1] = v->point()[1];
			if(v->point()[2] < p_min[2]) p_min[2] = v->point()[2];
		}
		i++;
	}

	return p_min;
}

Point Solid::Max_corner()
{
	SolidVertexIterator viter(this);
	Point p_max;
	int i = 0;
	for(;!viter.end(); ++ viter)
	{
		Vertex * v = * viter;
		if(i == 0) p_max = v->point();
		else
		{
			if(v->point()[0] > p_max[0]) p_max[0] = v->point()[0];
			if(v->point()[1] > p_max[1]) p_max[1] = v->point()[1];
			if(v->point()[2] > p_max[2]) p_max[2] = v->point()[2];
		}
		i++;
	}

	return p_max;
}




double Solid::volume()
{
	SolidVertexIterator viter(this);
	SolidFaceIterator fiter(this);
	Point p_min, p_orig, p[3];
	int i = 0;
	double pos_sum = 0.0, neg_sum = 0.0, vol;

	/*for(;!viter.end(); ++ viter)
	{
		Vertex * v = * viter;
		if(i == 0) p_min = v->point();
		else
		{
			if(v->point()[0] < p_min[0]) p_min[0] = v->point()[0];
			if(v->point()[1] < p_min[1]) p_min[1] = v->point()[1];
			if(v->point()[2] < p_min[2]) p_min[2] = v->point()[2];
		}
		i++;
	}

	p_orig = p_min - Point(1.0,1.0,1.0);*/
	p_orig = center_mass();

	for(; !fiter.end(); ++ fiter)
	{
		Face * f = * fiter;		
		i = 0;
		for(FaceVertexIterator fviter(f); ! fviter.end(); ++ fviter)
		{
			Vertex * v = * fviter;
			p[i] = v->point();
			i++;
		}
		//tetrahedral volume
		vol = fabs((p[0]-p_orig)*((p[1]-p_orig)^(p[2]-p_orig)))/6.0;
		if((f->norm())*(p_orig - p[0]) > 0.0) pos_sum += vol;
		else neg_sum += vol;
	}

	return fabs(pos_sum - neg_sum);
}



void Solid::reverse_vertex_order()
{
	SolidHalfEdgeIterator hviter(this);
		for(; !hviter.end(); ++ hviter)
		{
			HalfEdge * he1 = *hviter;
			Vertex * temp = he1->source();
			he1->source() = he1->target();
			he1->target() = temp;
		}
}

void Solid::close_boundary(int *list, int & max_id, int & max_id_f,  bool set_fill)
{
	int N = list[0];
//	int max_id = 0, max_id_f = 0;

	Point p;

	if(max_id == 0)
	{
		AVL::TreeIterator<Vertex> viter( m_verts );
		for ( ; !viter.end(); ++viter )
		{
			Vertex *v = *viter;
			max_id = std::max(max_id,v->id());
		}
	}
	max_id++;

	if(max_id_f == 0)
	{
		AVL::TreeIterator<Face> fiter( m_faces );
		for ( ; !fiter.end(); ++fiter )
		{
	        Face *f = *fiter;
			max_id_f = std::max(max_id_f,f->id());
		}
	}	
	//max_id_f++;

	for(int i = 1; i <= N; i++)
		p += idVertex(list[i])->point();

	p /= (double)N;

	tVertex v_new  = createVertex( max_id);
    v_new->point() = p;
    v_new->id()    = max_id;
	v_new->id2()   = numVertices();

	int vv[3];

	//printf("new vertex created\n");
	//p.print();
	//printf("id = %i, id2 = %i\n", v_new->id(),v_new->id2());

	//printf("N = %i\n",N);

	for(int j = 0; j < N; j++)
	{
		max_id_f++;

        vv[0] = list[j + 1];
        vv[1] = list[(j+1)%N + 1];
		vv[2] = max_id;
	//	printf("about to create face\n");
        tFace f = createFace(vv, max_id_f );   
		idVertex(list[j+1])->boundary() = false;

		std::string key("fill");
		std::string val(" ");

		if (set_fill)
			Trait::updateTraitString(f->string(), key, val);	 //setting fill parameter

	//	printf("on j = %i\n",j);
	//	printf("Face %i  %i %i %i created\n",f->id(),vv[0],vv[1],vv[2]);
		
//		break;
    }
//	labelBoundaryEdges();

	
}

void Solid::close_boundaries(bool set_fill)
{
	int *list = new int[numVertices()];
	int max_id = 0, max_id_f = 0;
	int bnd = 0;

	AVL::TreeIterator<Vertex> viter( m_verts );
    for ( ; !viter.end(); ++viter )
	{
		Vertex *v = *viter;

		int N = 0;

		if(v->boundary())
		{
			printf("boundary found\n");
			N = 1;
			list[N] = v->id();
			list[0] = N;
			
			for(VertexVertexIterator vviter(v); !vviter.end(); ++vviter)
			{
				Vertex *vv = *vviter;
				Edge *e = idEdge(v->id(),vv->id());


				if(e->boundary())
				{
					N++;
					list[0] = N;
					
					HalfEdge *he = e->halfedge(0);

					if(he->source()->id() == vv->id())
						list[N] = vv->id();
					else
					{
						list[N] = v->id();
						list[N-1] = vv->id();
					}
					break;
				}
			}

			printf("%i\n%i\n", list[1], list[2]);

			while(list[1] != list[N])
			{
				Vertex *vv = idVertex(list[N]);
		/*		printf("on vertex %i\n",vv->id());

				if(vv->boundary())
					printf("boundary vert\n");
				else
					printf("NOT boundary vert\n");*/

				if(0)//vv->id() == 3864)
				{
					//printf("testing 3864\n");
					HalfEdge *he = vv->halfedge();


					for(SolidEdgeIterator eiter(this); !eiter.end(); ++eiter)
					{
						Edge * e = *eiter;

						if(e->include_vertex(vv))
						{
							printf("vert %i included\n",e->other_vertex(vv)->id());
							if(e->boundary())printf("boundary\n");
						}
						
					}

					/*isBoundary(vv) = false;

					for(VertexOutHalfedgeIterator vihiter(this,vv); !vihiter.end(); ++vihiter)
					{
						HalfEdge *he = *vihiter;
						printf("vert by halfedge %i\n",he->target()->id());
					}*/

					exit(0);


//					he->set_verbose(true);

					if(he == NULL)
						printf("no halfedge\n");
					
					printf("past he 1\n");

					he = vv->most_clw_out_halfedge();
					if(he == NULL)
						printf("no most ccw halfedge\n");

					printf("past he 2\n");




				}
			//	else
			//		exit(0);

				for(VertexVertexIterator vviter(vv); !vviter.end(); ++vviter)
			//	for(VertexInHalfedgeIterator vihiter(this,vv); !vihiter.end(); ++vihiter)
				{
					Vertex *v_cur = *vviter;
					Edge *e = idEdge(v_cur->id(),vv->id());
			//		HalfEdge *he = *vihiter;
			//		printf("on Vvertex %i,\n",vn++);
			//		Vertex *v_cur = he->target();
			//		printf("on Vvertex %i,\n",v_cur->id());
					//if(he->edge()->boundary() && v_cur->id() != list[N-1])
					if(e->boundary() && v_cur->id() != list[N-1])
					{
						N++;
						list[N] = v_cur->id();
				//		printf("%i\n",list[N]);
						break;
					}
				}
				if(N == 2)
					break;
			}
			

			printf("boundary found, N = %i\n",N);
			printf(" first id = %i, last id  = %i\n",list[1],list[N]);

			if(N > 2)
				N--;
			else
			{
				printf("2-vertex boundary!\n");
			}

	//		printf("Now, N = %i\n",N);
	//		printf(" first id = %i, last id  = %i\n",list[1],list[N]);
			bnd++;

			printf("boundary %i has %i edges\n", bnd,N);
			list[0] = N;
			
			
			close_boundary(list, max_id, max_id_f, set_fill);
			//break;

		}
		
	}

	delete [] list;

}


void Solid::keep_largest_connected_piece(bool mark_verts_only)
{
	int *cur_list_he = new int[numFaces() + 1];
	int *component_index = new int[numFaces()];
	int *component_size = new int[numFaces()];
	bool *visited = new bool[numFaces()];
	int *id_f = new int[numFaces()];

	int num_visited = 0, N_comps = 0, max_id = 0;
	for(int i = 0; i < numFaces(); i++ )		
	{
		visited[i] = false;
		component_size[i] = 0;
	}
		
	for(SolidFaceIterator fiter(this); !fiter.end(); ++fiter)		
	{
		Face * f = *fiter;
		if(max_id < f->id())
			max_id = f->id();
	}

	int count = 0;
	int *ID2 = new int[max_id+1];

	for(SolidFaceIterator fiter(this); !fiter.end(); ++fiter)		
	{
		Face * f = *fiter;
		ID2[f->id()] = count;
		count++;
	}

	count = 0;
	for(SolidFaceIterator fiter(this); !fiter.end(); ++fiter)		
	{
		Face * f = *fiter;
		int I = ID2[f->id()];
	//	int I = f->id()-1;
	//	int I = f->id2()-1;

		id_f[count] = f->id();
		count++;
		
		if(!visited[I])
		{	
			//printf("first I in new = %i\n",I);

			N_comps++;
			component_index[I] = N_comps - 1;
			component_size[N_comps-1]++;

			visited[I] = true;
			num_visited++;
			cur_list_he[0] = 1;
			cur_list_he[cur_list_he[0]] = f->id();

			for(int n = 1; n <= cur_list_he[0]; n++)
			{						
				Face * cur_F = idFace(cur_list_he[n]);
				visited[ID2[cur_F->id()]] = true;//needed here?
			//	visited[cur_F->id()-1] = true;//needed here?
			//	visited[cur_F->id2()-1] = true;//needed here?

			//	printf("cur_F id2 - 1 = %i\n",cur_F->id2()-1);

				
				for( FaceHalfedgeIterator fhiter(cur_F); !fhiter.end(); ++ fhiter )//populating new cluster
				{
					HalfEdge * he = *fhiter;//face belogning to edge
					HalfEdge * he_sym = he->he_sym();
		
					if(he_sym == NULL)
						continue;

					Face * ff = he_sym->face();

					int II = ID2[ff->id()];
				//	int II = ff->id()-1;
				//	int II = ff->id2()-1;
						
					if(!visited[II])
					{
						cur_list_he[0]++;
						cur_list_he[cur_list_he[0]] = ff->id();
						visited[II] = true;
						component_index[II] = N_comps - 1;
						component_size[N_comps-1]++;
						num_visited++;	
				//		printf("Face id2 %i added\n",ff->id2()-1);
				//		printf("total n = %i \n",n);

					}
				}
			}
		}
	}

	//printf("past labeling loop\n");


	int max_ind = 0, max_comps = 0;

	for(int i = 0; i < N_comps; i++)
	{
		printf("comp %i size = %i\n", i, component_size[i]);
		if(component_size[i] > max_comps)
		{
			max_comps = component_size[i];
			max_ind = i;
		}
	}

//	printf("max comp found\n");


	int n_face_orig = numFaces();
//	Edge ** ee = new Edge*[3];

	DList<Edge>  Elist;


	for (int i = 0; i < n_face_orig; i++)
	{
	//	printf("on face %i out of %i\n",i,n_face_orig);//i == 25604
		//printf("component id %i\n",component_index[i]);
		if(component_index[i] != max_ind)
		{

			tFace f = idFace(id_f[i]);
			tHalfEdge he = f->halfedge();      

			for(FaceEdgeIterator feiter(f); !feiter.end(); ++feiter)
			{
				Edge *e = *feiter;
		//		ee[j] = e;
		//		j++;
				if(!Elist.contains(e))
					Elist.Append(e);
			}

	//		if(i == 25604)
	//			printf("set ee\n");

			
			for(FaceVertexIterator fviter(f); !fviter.end(); ++fviter)
			{
				Vertex * v = *fviter;
				 if ( v->trait() != NULL )
					   v->trait()->clear( v->trait() );
				
				 if(v->halfedge() != NULL)
					v->halfedge() = NULL;	
			}

//			if(i == 25604)
//				printf("set verts to null\n");

			List<HalfEdge> hes;

			do
			{
				he = he->he_next();
				hes.Append(he);
			}while ( he != f->halfedge() );

	//		if(i == 25604)
	//			printf("halfedge list filled\n");

			ListNode<HalfEdge> * node = hes.head();
			do
			{
				if ( node->data()->trait() != NULL )
				{
					node->data()->trait()->clear( node->data()->trait() );
				}
				delete node->data();
				node = node->next();
			}while ( node != NULL );

	//		if(i == 25604)
	//			printf("halfedges in list nulled\n");

			if ( f->trait() != NULL )
			{
				f->trait()->clear( f->trait() );
			}

		//	if(i == 25604)
		//		printf("face traits removed\n");

			m_faces.remove(f);

			if(f != NULL)
				delete f;

	//		if(i == 25604)
	//			printf("faces removed\n");

			//if(0)
			//for(j=0;j<3;j++)
			//{
			//	tEdge e = ee[j];
			//	if(e != NULL)
			//	{
			//		if(e->ekey().s() != 0 && e->ekey().t() != 0)//not allowed 0 vertex id's in CCBBM
			//		{
			//			if ( e->trait() != NULL )
			//				e->trait()->clear( e->trait() );
			//			m_edges.remove(e);
			//		}
			//		delete e;
			//	}
			//}

	//		if(i == 25604)
	//			printf("edges removed\n");
		}
	}

	printf("thru main loop\n");
	
	for( DListIterator<Edge> Eiter( Elist ); !Eiter.end(); ++ Eiter )
	{
		Edge *e = * Eiter;

		if(e != NULL)
		{
			if(e->ekey().s() != 0 && e->ekey().t() != 0)//not allowed 0 vertex id's in ccbbm
			{
				if ( e->trait() != NULL )
					e->trait()->clear( e->trait() );
				m_edges.remove(e);
			}
			delete e;
		}
	}


	printf("thru edge loop\n");

	removeDanglingVertices();
	reset_id2s();




	delete [] cur_list_he ;
	delete [] component_index ;
	delete [] component_size ;
	delete [] visited ;
	delete [] id_f ;
//	delete [] ee ;
	delete [] ID2 ;


//	close_boundaries();


}

//removes connected components that consist of only two faces & no boundaries
void Solid::removeDegenerateComponents()
{
	bool * degenerate = new bool[numFaces()];
	int * id_f = new int[numFaces()];
	int * id_f_sym = new int[numFaces()];
	Edge ** ee = new Edge*[3];
	int i = 0;

	AVL::TreeIterator<Face> fiter( m_faces );
    for ( ; !fiter.end(); ++fiter )
    {
        tFace f = *fiter;
        tHalfEdge he = f->halfedge();
        tHalfEdge he_n = he->he_next();
        tHalfEdge he_p = he->he_prev();

		int sym_id = he->he_sym()->face()->id();

		id_f[i] = f->id();
		id_f_sym[i] = sym_id;

		if(!he->edge()->boundary() && !he_p->edge()->boundary() && !he_n->edge()->boundary())
		{
			if( sym_id == he_n->he_sym()->face()->id() && sym_id == he_p->he_sym()->face()->id())
			{	degenerate[i] = true; printf("degenerate %i\n",f->id());}
			else
				degenerate[i] = false;
		}
		else
			degenerate[i] = false;
		i++;
	}

	int n_face_orig = numFaces();

	for (i = 0; i < n_face_orig; i++)
	if(degenerate[i] && id_f_sym[i] > id_f[i])
	{
		int j = 0;

		tFace f = idFace(id_f[i]);
		tFace f_sym = f->halfedge()->he_sym()->face();

        tHalfEdge he = f->halfedge();
        tHalfEdge he_sym = f_sym->halfedge();


		for(FaceEdgeIterator feiter(f); !feiter.end(); ++feiter)
		{
			Edge *e = *feiter;
			ee[j] = e;
			j++;
		}

		
		for(FaceVertexIterator fviter(f); !fviter.end(); ++fviter)
		{
			Vertex * v = *fviter;
			 if ( v->trait() != NULL )
                   v->trait()->clear( v->trait() );
			 m_verts.remove(v);
			 delete v;
        }

		List<HalfEdge> hes;
		List<HalfEdge> hes_sym;
        do
        {
            he = he->he_next();
            hes.Append(he);
        }while ( he != f->halfedge() );

		 do
        {
            he_sym = he_sym->he_next();
            hes_sym.Append(he_sym);
        }while ( he_sym != f_sym->halfedge() );

        ListNode<HalfEdge> * node = hes.head();
        do
        {
            if ( node->data()->trait() != NULL )
            {
                node->data()->trait()->clear( node->data()->trait() );
            }
            delete node->data();
            node = node->next();
        }while ( node != NULL );

		ListNode<HalfEdge> * node_sym = hes_sym.head();
        do
        {
            if ( node_sym->data()->trait() != NULL )
            {
                node_sym->data()->trait()->clear( node_sym->data()->trait() );
            }		
            delete node_sym->data();
            node_sym = node_sym->next();
        }while ( node_sym != NULL );


        if ( f->trait() != NULL )
        {
            f->trait()->clear( f->trait() );
        }
		if ( f_sym->trait() != NULL )
        {
            f_sym->trait()->clear( f_sym->trait() );
        }

		m_faces.remove(f);
		m_faces.remove(f_sym);


        delete f;
        delete f_sym;


		for(j=0;j<3;j++)
		{
			tEdge e = ee[j];
		    if ( e->trait() != NULL )
			{
				e->trait()->clear( e->trait() );
			}
			m_edges.remove(e);
			delete e;
		}
	}
	
	delete [] degenerate;
	delete [] ee;
	delete [] id_f;
	delete [] id_f_sym;
}




		




						
						





	



Solid::~Solid()
{
	bool verb = false;

	if(0)//numVertices() > 65000)
		verb = false;
	if(verb)
	printf("in destructor SOLID\n");
    //remove vertices
    AVL::TreeIterator<Vertex> viter( m_verts );
    for ( ; !viter.end(); ++viter )
    {
        tVertex v = *viter;
		if(verb)printf("v->id = %i\n",v->id());
        if ( v->trait() != NULL )
        {
            v->trait()->clear( v->trait() );
        }
        delete v;
    }

	if(verb)
	printf("thru verts\n");

    //remove faces
    AVL::TreeIterator<Face> fiter( m_faces );
    for ( ; !fiter.end(); ++fiter )
    {
        tFace f = *fiter;

        tHalfEdge he = f->halfedge();

        List<HalfEdge> hes;

        do
        {
            he = he->he_next();
            hes.Append(he);
        }while ( he != f->halfedge() );

        ListNode<HalfEdge> * node = hes.head();
        do
        {
            if ( node->data()->trait() != NULL )
            {
                node->data()->trait()->clear( node->data()->trait() );
            }
            delete node->data();
            node = node->next();
        }while ( node != NULL );

        if ( f->trait() != NULL )
        {
            f->trait()->clear( f->trait() );
        }

        delete f;
    }

	if(verb)
	printf("thru faces\n");


    //remove edges
    AVL::TreeIterator<Edge> eiter( m_edges );
    for ( ; !eiter.end(); ++ eiter )
    {
        tEdge e = *eiter;
        if ( e->trait() != NULL )
        {
            e->trait()->clear( e->trait() );
        }
        delete e;
    }

	if(verb)
	printf("thru edges\n");
	
	//if(m_F != NULL)
	//	delete [] m_F;
	//if(m_V != NULL)
	//	delete [] m_V;
	//if(m_use_face != NULL)
	//	delete [] m_use_face;

	if(verb)
		printf("done destructor SOLID\n");


};

void Solid::write( std::ostream & os )
{

	char buffer_array[65536];
    os.rdbuf()->pubsetbuf(buffer_array, sizeof(buffer_array));

    //remove vertices
    AVL::TreeIterator<Vertex> viter( m_verts );
    for ( ; !viter.end() ; ++ viter )
    {
        tVertex v = *viter;

        os << "Vertex " << v->id();

        for ( int i = 0; i < 3; i ++ )
        {
            os << " " << v->point()[i];  
        }

        if ( v->string().size() > 0 )
        {
            os << " {" << v->string() << "}";
        }

        os << std::endl;
    }

    for ( AVL::TreeIterator<Face> fiter(m_faces); !fiter.end(); ++ fiter )
    {
        tFace f = *fiter;
        os << "Face " << f->id();
        tHalfEdge he = f->halfedge();
        do
        {
            os << " " << he->target()->id();
            he = he->he_next();
        }while ( he != f->halfedge() );

        if ( f->string().size() > 0 )
        {
            os << " {" << f->string() << "}";
        }

        os << std::endl;
    }


    for ( AVL::TreeIterator<Edge> eiter(m_edges); !eiter.end(); ++ eiter )
    {
        tEdge e = *eiter;
        if ( e->string().size() > 0 )
        {
            os << "Edge " << edgeVertex1(e)->id()<< " " << edgeVertex2(e)->id() << " ";
            os << "{" << e->string() << "}" << std::endl;
        }
    }

    for ( AVL::TreeIterator<Face> pfiter(m_faces); !pfiter.end(); ++ pfiter )
    {
        tFace f = *pfiter;
        tHalfEdge he = f->halfedge();
        do
        {
            if ( he->string().size() > 0 )
            {
                os << "Corner " << he->vertex()->id() << " " << he->face()->id() << " ";
                os <<"{" << he->string() << "}" << std::endl;
            }
            he = he->he_next();
        }while ( he != f->halfedge() );
    }


};

void Solid::read( std::istream & is )
{

    char line[MAX_LINE];
    int id, count = 0;

    while ( is && !is.eof() && is.getline(line, MAX_LINE) )
    {
        if ( strlen( line ) == 0 ) continue;

        std::string s(line);

        string_token_iterator iter(s, " \n"); 

        std::string str = *iter;

        if ( str == "Vertex" )
        {
			count++;

            str = *(++iter);

            id = atoi( str.c_str() );


            Point p;
            for ( int i = 0 ; i < 3; i ++ )
            {
                str = *(++iter);
                p[i] = atof( str.c_str() );
            }

            tVertex v  = createVertex( id );

            v->point() = p;
            v->id()    = id;
			v->id2()   = count;

            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
            {
                v->string() = s.substr( sp+1, ep-sp-1 );
            }

            continue;

        }

        if ( str == "Face" )
        {
            str = *(++iter);

            id = atoi( str.c_str() );

            int v[3];
            for ( int i = 0; i < 3; i ++ )
            {
                str = *(++iter);
                v[i] = atoi( str.c_str() );
            }

            tFace f = createFace( v, id );

            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
            {
                f->string() = s.substr( sp+1, ep-sp-1 );
            }
            continue;
        }

        //read in edge attributes
        if ( str == "Edge" )
        {
            str = *(++iter);
            int id0 = atoi( str.c_str() );

            str = *(++iter);
            int id1 = atoi( str.c_str() );

            tEdge edge = idEdge( id0, id1 );

            str = *(++iter);

            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
            {
                edge->string() = s.substr( sp+1, ep-sp-1 );
            }

            continue;

        }

        //read in edge attributes
        if ( str == "Corner" )
        {
            str = *(++iter);
            int id0 = atoi( str.c_str() );

            str = *(++iter);
            int id1 = atoi( str.c_str() );


            Vertex * v = idVertex( id0 );
            Face   * f = idFace( id1 );
            tHalfEdge he = corner( v, f );

            str = *(++iter);

            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
            {
                he->string() = s.substr( sp+1, ep-sp-1 );
            }
            continue;
        }


    }

    labelBoundaryEdges();

    removeDanglingVertices();

	reset_id2s();

};

//void Solid::add_texture(char * fn_bmp)//based on bitmap files, only for alias obj files for now (11/11/13)
//{
//	BitMap BM;
//
//	BM.read(fn_bmp);
//
//	std::string u_key = std::string("u_text");
//	std::string v_key = std::string("v_text");
//
//	for(SolidVertexIterator viter(this); !viter.end(); ++viter)
//	{
//		Vertex *v = *viter;
//		double U, V;
//		std::string  U_str = MeshLib::Trait::getTraitValue(v->string(), u_key);
//		std::string  V_str = MeshLib::Trait::getTraitValue(v->string(), v_key);
//		sscanf(U_str.c_str(),"%lf",&U);	
//		sscanf(V_str.c_str(),"%lf",&V);	
//
//		Point RGB = BM.interpolate(U,1.0 - V);
//
//		v->color_vertex(RGB[2]/255.0,RGB[1]/255.0,RGB[0]/255.0);
//	}
//}


void Solid::read_NM( const char * input, Point *test, int n_test_pts )
{
	std::string fname = input;
	std::string format = fname.substr(fname.length()-3,3);



	if(!strcmp(format.c_str(),"ply"))
	{
		readPLY_NM(input);
		return;
	}
	else if(!strcmp(format.c_str(),"off"))
	{
		readOFF(input);
		return;
	}
	else if(!strcmp(format.c_str(),"obj"))
	{
		readOBJ(input);
		return;
	}

	char * cline = (char*)calloc(1024,sizeof(char));
	std::string sline;


    FILE * is = fopen( input, "r" );
    if ( is == NULL )
    {
        fprintf(stderr,"Error is opening file %s\n", input );
        return;
    }

    char line[MAX_LINE];
    int count = 0, fcount = 0;

	int v_s = -1, v_counter = 0;
	int f_s = -1, f_counter = 0;

	//find the number of verts, faces, attributes
	while(!feof(is))
	{
		fgets(cline,1024,is);
		sline = cline;
//		printf("%s",cline);
		v_s  = sline.rfind("Vertex ",0);
		f_s  = sline.rfind("Face ",0);

		if(v_s == 0)	v_counter++;
		if(f_s == 0)	f_counter++;

	}

	if(v_s == 0)	v_counter--;
	if(f_s == 0)	f_counter--;


	bool * use_face = new bool[f_counter];
	for(int i = 0; i < f_counter; i++) use_face[i] = true;

	Point *V = new Point[v_counter];
	int * v_id = new int[v_counter];

	int * F_id = new int[f_counter];
	Point * F = new Point[f_counter];

	std::string * vstring = new std::string[v_counter];
	std::string * fstring = new std::string[f_counter];

	int max_vid = 0;

	rewind(is);
	
	while ( !feof(is) )
    {
        fgets( line, MAX_LINE, is );
		char * str = strtok( line, " \r\n\t");
        if ( feof(is) || strlen( line ) == 0 ) break;

		if ( strcmp(str, "Vertex" ) == 0 )
        {
			count++;
            str = strtok(NULL," \r\n\t");

            v_id[count-1] = atoi( str );

			if(max_vid <  v_id[count-1])
				max_vid =  v_id[count-1];
		
            for ( int i = 0 ; i < 3; i ++ )
            {
                str = strtok(NULL," \r\n\t");
                V[count-1][i] = atof( str );
            }

            str = strtok( NULL, "\r\n\t");

            if ( str == NULL ) continue;

            std::string s(str);
            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
            {
                vstring[count-1] = s.substr( sp+1, ep-sp-1 );
            }

   //         continue;

        }
	}

//	printf("past vertex reading\n");
	
	int *orig_vids = new int[max_vid+1];
	for(int kk = 0; kk < max_vid+1; kk++)
		orig_vids[kk] = 0;
	for(int kk = 0; kk < v_counter; kk++)
		orig_vids[ v_id[kk]] = kk+1;
	
	rewind(is);
	
	while ( !feof(is) )
    {
        fgets( line, MAX_LINE, is );
		char * str = strtok( line, " \r\n\t");
        if ( feof(is) || strlen( line ) == 0 ) break;

        if ( strcmp(str,"Face") == 0 )
		{
	
			fcount++;	
            str = strtok(NULL, " \r\n\t");

            F_id[fcount-1] = atoi( str );
            for ( int i = 0; i < 3; i ++ )
            {
                str = strtok(NULL," \r\n\t");
           
				int ind = atoi( str );		
				F[fcount-1][i] = (float) orig_vids[ind]; 
				
				if(orig_vids[ind] == 0)
					printf("orig_vids[%i] is 0, face %i\n",ind,fcount-1);

            }			
		
            str = strtok( NULL, "\r\n\t");
            if ( str== NULL || strlen( str ) == 0 ) continue;

            std::string s(str);

            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
				fstring[fcount-1] = s.substr( sp+1, ep-sp-1 );          

		}
	}

//	printf("past face reading\n");

	if(n_test_pts > 0)
		label_non_manifold_faces_test(V, F,  use_face, v_counter,f_counter, test,  n_test_pts);
		
	else
		label_non_manifold_faces(V, F,  use_face, v_counter,f_counter);

//	printf("past nm labeling\n");


	for(int i = 0; i < v_counter; i++)
	{
		tVertex v  = createVertex( v_id[i] );
        v->point() = V[i];
        v->id()    = v_id[i];
		v->string() = vstring[i];
	}

	for(int i = 0; i < f_counter; i++)
		if(use_face[i])
		{
		//	printf("using face %i\n",i);
			int v[3];
			for(int j = 0; j < 3; j++)
				v[j] = v_id[(int)F[i][j]-1];//Floc[i][j];
			
			tFace f  = createFace(v, F_id[i] );
			f->string() = fstring[i];
		}

	

	//	printf("mesh created\n");

	rewind(is);
	
	while ( !feof(is) )
    {
		fgets( line, MAX_LINE, is );
		char * str = strtok( line, " \r\n\t");
        if ( feof(is) || strlen( line ) == 0 ) break;

        //read in edge attributes
        if ( strcmp(str,"Edge")==0 )
        {
            str = strtok(NULL, " \r\n\t");
            int id0 = atoi( str );

            str = strtok(NULL, " \r\n\t");
            int id1 = atoi( str );

            tEdge edge = idEdge( id0, id1 );


            str = strtok(NULL, "\r\n\t");

            std::string s(str);
            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
            {
                edge->string() = s.substr( sp+1, ep-sp-1 );
            }

            continue;

        }


        //read in edge attributes
        if ( strcmp(str,"Corner")==0 )
        {
            str = strtok(NULL," \r\n\t");
            int id0 = atoi( str );

            str = strtok(NULL," \r\n\t");
            int id1 = atoi( str );


            Vertex * v = idVertex( id0 );
            Face   * f = idFace( id1 );
            tHalfEdge he = corner( v, f );

            str = strtok(NULL,"\r\n\t");
            std::string s(str);

            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
            {
                he->string() = s.substr( sp+1, ep-sp-1 );
            }
            continue;
        }


    }


    fclose(is);


    labelBoundaryEdges();
    removeDanglingVertices();
	reset_id2s();


	delete [] F;
	delete [] V;
	delete [] use_face;
	delete [] F_id;
	delete [] v_id;
	
	delete [] vstring;
	delete [] fstring;
	delete [] orig_vids;

	//printf("exiting read_NM\n");

};


void Solid::read(const char * input, bool reset_ids, int v_id_start, int f_id_start)
{
//	read_NM(input);
//	return;

	//testing grounds for anything

	//testing sparse storage structure SparseMatrix2

	//SparseMatrix2<double> SP_double(4,5,0.0);


//	printf("v start, f start %i %i\n", v_id_start, f_id_start);

//	v_id_start = 5000;
//	f_id_start = 5000;


	
	std::string fname = input;
	std::string format = fname.substr(fname.length()-3,3);

	if(!strcmp(format.c_str(),"ply"))
	{
		readPLY_NM(input);
		return;
	}
	else if(!strcmp(format.c_str(),"off"))
	{
		readOFF(input);
		return;
	}
	else if(!strcmp(format.c_str(),"obj"))//already corrects for NM topology						   
	{
		readOBJ(input);
		return;
	}

	


    FILE * is = fopen( input, "r" );
    if ( is == NULL )
    {
        fprintf(stderr,"Error is opening file %s\n", input );
        return;
    }

    char line[MAX_LINE];
    int id, id_f, count = 0, count_f = 0;


	int * old_ids = NULL;
	int * old_ids_f = NULL;

	if(reset_ids)
	{
		int max_id = 1;
		int max_id_f = 1;
		while ( !feof(is) )
		{
			fgets( line, MAX_LINE, is );
			char * str = strtok( line, " \r\n\t");
			if ( feof(is) || strlen( line ) == 0 ) break;
			if ( strcmp(str, "Vertex" ) == 0 )
			{			
				str = strtok(NULL," \r\n\t");
				id = atoi( str );
				if(max_id < id)
					max_id = id;
			}
			if ( strcmp(str, "Face" ) == 0 )
			{			
				str = strtok(NULL," \r\n\t");
				id_f = atoi( str );
				if(max_id_f < id_f)
					max_id_f = id_f;
			}
		}
		old_ids = new int[max_id+1];	
		old_ids_f = new int[max_id_f+1];	
	}

	rewind(is);

	while ( !feof(is) )
    {
        fgets( line, MAX_LINE, is );
		std::string s(line);

		char * str = strtok( line, " \r\n\t");

        if ( feof(is) || strlen( line ) == 0 ) break;

		 if ( strcmp(str, "Vertex" ) == 0 )
        {
			count++;
            str = strtok(NULL," \r\n\t");

            id = atoi( str );

			if(reset_ids)
			{
				old_ids[id] = count;
				id = count;			
			}


            Point p;
            for ( int i = 0 ; i < 3; i ++ )
            {
                str = strtok(NULL," \r\n\t");
                p[i] = atof( str );
            }

			id += v_id_start;
            tVertex v  = createVertex( id );

            v->point() = p;
            v->id()    = id;
			v->id2()   = count;

            str = strtok( NULL, "\r\n\t");

            if ( str == NULL ) continue;

            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
            {
                v->string() = s.substr( sp+1, ep-sp-1 );
            }

            continue;

        }


        if ( strcmp(str,"Face") == 0 )
		{
	//	if(tru_face[face_iter])
			count_f++;
      //  {
		

            str = strtok(NULL, " \r\n\t");

            id = atoi( str );

			if(reset_ids)
			{
				old_ids_f[id] = count_f;
				id = count_f;			
			}


            int v[3];
            for ( int i = 0; i < 3; i ++ )
            {
                str = strtok(NULL," \r\n\t");
                v[i] = atoi( str );
				if(reset_ids)
					v[i] = old_ids[v[i]];

				v[i] += v_id_start;
            }


			tFace f = createFace(v, id + f_id_start);


            str = strtok( NULL, "\r\n\t");
            if ( str== NULL || strlen( str ) == 0 ) continue;

			//std::string s(str);
            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
            {
                f->string() = s.substr( sp+1, ep-sp-1 );
            }

			if (0)//count_f == 3)
			{
				printf("s = %s\n", s.c_str());
				printf("in read, verts = %i\n string: %s\n", numVertices(), f->string().c_str());
			}
            
	//	}
	//	face_iter++;
		continue;
		}

        //read in edge attributes
        if ( strcmp(str,"Edge")==0 )
        {
            str = strtok(NULL, " \r\n\t");
            int id0 = atoi( str );

            str = strtok(NULL, " \r\n\t");
            int id1 = atoi( str );

            tEdge edge = idEdge( id0, id1 );


            str = strtok(NULL, "\r\n\t");

     //       std::string s(str);
            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
            {
                edge->string() = s.substr( sp+1, ep-sp-1 );
            }

            continue;

        }

        //read in edge attributes
        if ( strcmp(str,"Corner")==0 )
        {
            str = strtok(NULL," \r\n\t");
            int id0 = atoi( str );

            str = strtok(NULL," \r\n\t");
            int id1 = atoi( str );

			if(reset_ids)
				id0 = old_ids[id0];


            Vertex * v = idVertex( id0 );
            Face   * f = idFace( id1 );
            tHalfEdge he = corner( v, f );

            str = strtok(NULL,"\r\n\t");
   //         std::string s(str);

            int sp = (int)s.find("{");
            int ep = (int)s.find("}");

            if ( sp >= 0 && ep >= 0 )
            {
                he->string() = s.substr( sp+1, ep-sp-1 );
            }
            continue;
        }
    }

    fclose(is);

    labelBoundaryEdges();
    removeDanglingVertices();
	reset_id2s();

	if(reset_ids)
	{
		delete [] old_ids;
		delete [] old_ids_f;
	}


	//printf("string face 1: %s\n",idFace(1)->string().c_str());

};


void Solid::readOBJ( const char * input)
{

  //TODO:...what type of OBJ is this!?
  // Fix or rewrite.
  printf("The MeshLib::Solid::readOBJ function is suspicious enough that "
    "I have disabled it. DMoyer 171129");
  exit(1);

/*
	printf("reading ascii mtl obj file...\n");fflush(stdout);
    FILE * is = fopen( input, "r" );
    if ( is == NULL )
    {
        fprintf(stderr,"Error is opening file %s\n", input );
        return;
    }

//    char line[MAX_LINE];
    int id;//, verts, faces;

	rewind(is);

	char * cline = (char*)calloc(1024,sizeof(char));
	std::string sline;
	std::string substr;

	int v_s, v_counter = 0;
	int vn_s,vn_counter = 0;
	int vt_s, vt_counter = 0;
	int f_s, f_counter = 0;

	//find the number of verts, faces, attributes
	while(!feof(is))
	{
		fgets(cline,1024,is);
		sline = cline;

		v_s  = sline.rfind("v ",0);
		vn_s  = sline.rfind("vn ",0);
		vt_s  = sline.rfind("vt ",0);
		f_s  = sline.rfind("f ",0);

		if(v_s == 0)	v_counter++;
		if(vn_s == 0)	vn_counter++;
		if(vt_s == 0)	vt_counter++;
		if(f_s == 0)	f_counter++;
	}

	if(v_s == 0)	v_counter--;
	if(vn_s == 0)	vn_counter--;
	if(vt_s == 0)	vt_counter--;
	if(f_s == 0)	f_counter--;
	printf("Num v, vn, vt, f = %i %i %i %i\n",v_counter,vn_counter,vt_counter,f_counter);	
	
	Point *V = new Point[v_counter];
	Point *VN = new Point[vn_counter];
	Point *VT = new Point[vt_counter];
	Point *F = new Point[f_counter];
	Point *FT = new Point[f_counter];
	int * VT_v = new int[v_counter];

	bool * use_face = new bool[f_counter];
	for(int i = 0; i < f_counter; i++) {use_face[i] = true;}
	for(int i = 0; i < v_counter; i++) {VT_v[i] = 0;}
	
	
	rewind(is);
	
	//set position of first data entry
	do
	{
		fgets(cline,1024,is);
		sline = cline;
		v_s  = sline.rfind("v ",0);
		vn_s  = sline.rfind("vn ",0);
		vt_s  = sline.rfind("vt ",0);
		f_s  = sline.rfind("f ",0);
	}while(!feof(is)&& f_s != 0 && v_s != 0 && vn_s != 0 && f_s != 0);

       
	std::string sub;
	int tmp_v_counter = 0;
	int tmp_vn_counter = 0;
	int tmp_vt_counter = 0;
	int tmp_f_counter = 0;

	int div = 1;

	if(vt_counter > 0)
		div++;
	if(vn_counter > 0)
		div++;


	
	while ( !feof(is) )
    {
        Point p;

		if(sline.length() > 2)
			sub = sline.substr(2,sline.length());

		v_s  = sline.rfind("v ",0);
		vn_s  = sline.rfind("vn ",0);
		vt_s  = sline.rfind("vt ",0);
		f_s  = sline.rfind("f ",0);

		if(v_s == 0)
		{
			sscanf(sub.c_str(),"%lf %lf %lf",&V[tmp_v_counter][0],&V[tmp_v_counter][1],&V[tmp_v_counter][2]);	
			tmp_v_counter++;
		}
		if(vn_s == 0)
		{
			sscanf(sub.c_str(),"%lf %lf %lf",&VN[tmp_vn_counter][0],&VN[tmp_vn_counter][1],&VN[tmp_vn_counter][2]);	
			tmp_vn_counter++;
		}
		if(vt_s == 0)
		{			
			sscanf(sub.c_str(),"%lf %lf %lf",&VT[tmp_vt_counter][0],&VT[tmp_vt_counter][1]);	
			tmp_vt_counter++;
		}
		if(f_s == 0)
		{	
			if(0)
			{
			int sp = (int)sline.find("f ");
		    int ep = (int)sline.find("/",sp);
         
			substr = sline.substr( sp+1, ep-sp-1 );
			sscanf(substr.c_str(),"%lf",&F[tmp_f_counter][0]);

			sp = (int)sline.find(" ",ep);
			ep = (int)sline.find("/",sp);

			substr = sline.substr( sp+1, ep-sp-1 );
			sscanf(substr.c_str(),"%lf",&F[tmp_f_counter][1]);

			sp = (int)sline.find(" ",ep);
			ep = (int)sline.find("/",sp);

			substr = sline.substr( sp+1, ep-sp-1 );
			sscanf(substr.c_str(),"%lf",&F[tmp_f_counter][2]);
			}
			else
			{
				int jj = 0;
				int k_v = 0;
				char * cur = strtok(cline,"f /");
				int II;
			
				while(cur != NULL)
				{
					sscanf(cur,"%i",&II);
				
					if(jj%div == 0)//assuming geometry, texture and normal info for all faces
					{
						F[tmp_f_counter][k_v] = (float)II;
						k_v++;
					}
					if(jj%div == 1)
					{
						VT_v[(int)F[tmp_f_counter][k_v-1]-1] = II;
						FT[tmp_f_counter][k_v-1] = (float)II;
					}
								
					cur = strtok(NULL,"f /");	
					jj++;
				}	
			}

			tmp_f_counter++;
		}
	
		fgets(cline,1024,is);
		sline = cline;		
   }

	free(cline);

	printf("done reading\n");

	label_non_manifold_faces(V, F, use_face, v_counter,f_counter);


	for(int i = 0; i < v_counter; i++)
	{
		tVertex v  = createVertex( i+1 );
		v->point() = V[i];
		v->id()    = i + 1;
		v->id2()   = i + 1;

		if(vn_counter > i)
			v->normal() = VN[i];
		if(VT_v[i] != 0 && VT_v[i] <= vt_counter)
		{
			v->set_scalar_attribute(VT[VT_v[i]-1][0],"u_text");
			v->set_scalar_attribute(VT[VT_v[i]-1][1],"v_text");
		}
	}

	printf("verts set\n");


	for(int i = 0; i < f_counter; i++)
	{
		int v[3], vt[3];
		v[0] = (int)F[i][0];
		v[1] = (int)F[i][1];
		v[2] = (int)F[i][2];

		vt[0] = (int)FT[i][0];
		vt[1] = (int)FT[i][1];
		vt[2] = (int)FT[i][2];

		if(use_face[i])
		{
			tFace f = createFace( v, i+1 );

			if(vt[0] <= vt_counter && vt[1] <= vt_counter && vt[2] <= vt_counter &&
				vt[0] > 0 && vt[1] > 0  && vt[2] > 0)
				f->set_UV_attribute(VT[vt[0]-1],VT[vt[1]-1],VT[vt[2]-1],"UV");
		}
	}

	printf("faces set\n");

	fclose(is);

//	printf("file closed\n");


    labelBoundaryEdges();
	//printf("boundary labeled\n");

    removeDanglingVertices();
//	printf("dangling edges removed\n");

	reset_id2s();

	delete [] VN;
	delete [] VT;
	delete [] VT_v;
	delete [] V;
	delete [] F;
	delete [] FT;
	delete [] use_face;
*/
}





void Solid::readPLY( const char * input)
{
	printf("reading ascii ply file...\n");fflush(stdout);
    FILE * is = fopen( input, "r" );
    if ( is == NULL )
    {
        fprintf(stderr,"Error is opening file %s\n", input );
        return;
    }

//    char line[MAX_LINE];
	int id, verts, faces, d = 0, face_property = 0;

	rewind(is);

	char * cline = (char*)calloc(1024,sizeof(char));
	std::string sline;

	//find the number of verts
	do
	{
		fgets(cline,1024,is);
		sline = cline;
//		printf("%s",cline);
		d++;
	}while(!feof(is)&& sline.rfind("format ascii",0) == std::string::npos && d < 4);

	if(sline.rfind("format ascii",0) == std::string::npos)
	{
		printf("not an ascii ply file: this reader can only read ascii\n");
		exit(1);
	}

	//find the number of verts
	do
	{
		fgets(cline,1024,is);
		sline = cline;
//		printf("%s",cline);
	}while(!feof(is)&& sline.rfind("element vertex",0) == std::string::npos);

	if(sline.rfind("element vertex",0) == std::string::npos)
	{
		printf("invalid ply file: no 'element vertex' key found\n");
		exit(1);
	}

	std::string substr = sline.substr(15,sline.length()-15);
	sscanf(substr.c_str(),"%i",&verts);

	do
	{
		fgets(cline,1024,is);
		sline = cline;
//		printf("%s",cline);
	}while(!feof(is)&& sline.rfind("element face",0) == std::string::npos);

	if(sline.rfind("element face",0) == std::string::npos)
	{
		printf("invalid ply file: no 'element face' key found\n");
		exit(1);
	}

	substr = sline.substr(13,sline.length()-13);
	sscanf(substr.c_str(),"%i",&faces);

	do
	{
		fgets(cline,1024,is);
		sline = cline;
		if (sline.rfind("property", 0) != std::string::npos)
			face_property++;
//		printf("%s",cline);
	}while(!feof(is)&& sline.rfind("end_header",0) == std::string::npos);

	if(sline.rfind("end_header",0) == std::string::npos)
	{
		printf("invalid ply file: no 'end_header' key found\n");
		exit(1);
	}

	printf("%i verts, %i faces\n", verts, faces); fflush(stdout);
	printf("%i face properties\n", face_property); fflush(stdout);

	int counter = 0;

    while ( !feof(is) && counter < verts)
    {
		fgets(cline,1024,is);
   		id = counter + 1;
        Point p;
		sscanf(cline,"%lf %lf %lf",&p[0],&p[1],&p[2]);
        tVertex v  = createVertex( id );
        v->point() = p;
        v->id()    = counter + 1;
		v->id2()   = counter + 1;
	//	if(counter%1000 == 0)
	//	printf("counter = %i\n",counter);
		counter++;
   }
	printf("verts read\n");

	if(counter < verts)
	{
		printf("unexpected end of file in ply file\n");
		exit(1);
	}


	counter = 0;

    while ( !feof(is) && counter < faces)
    {
		fgets(cline,1024,is);
   		id = counter + 1;
		int v[3];
    //int dummy;
		char * str = strtok(cline, " \r\n\t");

		for (int k = 0; k < face_property + 2; k++)
		{
			str = strtok(NULL, " \r\n\t");
			//	puts(str);

			if ( k < face_property - 1 ){
				//dummy = atoi(str);
				continue;
			} else {
				v[k - face_property + 1] = atoi(str);
      }

		}
		v[0] += 1;
		v[1] += 1;
		v[2] += 1;

		//printf("%i %i %i %i\n",dummy, v[0]-1,v[1]-1,v[2]-1);

        createFace( v, id );
	//	if(counter%1000 == 0)

	//	printf("face counter = %i\n",counter);

		counter++;
   }

	printf("faces read\n");

	if(counter < faces)
	{
		printf("unexpected end of file in ply file\n");
		exit(1);
	}

	//printf("out of loop\n");
	
	fclose(is);

//	printf("file closed\n");


    labelBoundaryEdges();
	//printf("boundary labeled\n");

    removeDanglingVertices();
//	printf("dangling edges removed\n");

}

void Solid::readPLY_NM(const char * input)
{
	printf("reading ascii ply file...\n"); fflush(stdout);
	FILE * is = fopen(input, "r");
	if (is == NULL)
	{
		fprintf(stderr, "Error is opening file %s\n", input);
		return;
	}



	//    char line[MAX_LINE];
	int verts, faces, d = 0, face_property = 0;

	rewind(is);

	char * cline = (char*)calloc(1024, sizeof(char));
	std::string sline;

	//find the number of verts
	do
	{
		fgets(cline, 1024, is);
		sline = cline;
		//		printf("%s",cline);
		d++;
	} while (!feof(is) && sline.rfind("format ascii", 0) == std::string::npos && d < 4);

	if (sline.rfind("format ascii", 0) == std::string::npos)
	{
		printf("not an ascii ply file: this reader can only read ascii\n");
		exit(1);
	}

	//find the number of verts
	do
	{
		fgets(cline, 1024, is);
		sline = cline;
		//		printf("%s",cline);
	} while (!feof(is) && sline.rfind("element vertex", 0) == std::string::npos);

	if (sline.rfind("element vertex", 0) == std::string::npos)
	{
		printf("invalid ply file: no 'element vertex' key found\n");
		exit(1);
	}

	std::string substr = sline.substr(15, sline.length() - 15);
	sscanf(substr.c_str(), "%i", &verts);

	do
	{
		fgets(cline, 1024, is);
		sline = cline;
		//		printf("%s",cline);
	} while (!feof(is) && sline.rfind("element face", 0) == std::string::npos);

	if (sline.rfind("element face", 0) == std::string::npos)
	{
		printf("invalid ply file: no 'element face' key found\n");
		exit(1);
	}

	substr = sline.substr(13, sline.length() - 13);
	sscanf(substr.c_str(), "%i", &faces);

	do
	{
		fgets(cline, 1024, is);
		sline = cline;
		if (sline.rfind("property", 0) != std::string::npos)
			face_property++;
		//		printf("%s",cline);
	} while (!feof(is) && sline.rfind("end_header", 0) == std::string::npos);

	if (sline.rfind("end_header", 0) == std::string::npos)
	{
		printf("invalid ply file: no 'end_header' key found\n");
		exit(1);
	}

	printf("%i verts, %i faces\n", verts, faces); fflush(stdout);
	printf("%i face properties\n", face_property);

	int v_counter = verts, f_counter = faces;

	Point *V = new Point[v_counter];
	Point *F = new Point[f_counter]; 
	bool * use_face = new bool[f_counter];
	for (int i = 0; i < f_counter; i++) { use_face[i] = true; }



	int counter = 0;

	while (!feof(is) && counter < verts)
	{
		fgets(cline, 1024, is);
		//id = counter + 1;
		Point p;
		sscanf(cline, "%lf %lf %lf", &p[0], &p[1], &p[2]);

		V[counter] = p;
		//tVertex v = createVertex(id);
		//v->point() = p;
		//v->id() = counter + 1;
		//v->id2() = counter + 1;
		//	if(counter%1000 == 0)
		//		printf("counter = %i\n",counter);
		counter++;
	}
	printf("verts read\n");

	if (counter < verts)
	{
		printf("unexpected end of file in ply file\n");
		exit(1);
	}


	counter = 0;

	while (!feof(is) && counter < faces)
	{
		fgets(cline, 1024, is);
		//id = counter + 1;
		int v[3];
		//int dummy;
		Point P;

		char * str = strtok(cline, " \r\n\t");			

		for (int k = 0; k < face_property + 2; k++)
		{
			str = strtok(NULL, " \r\n\t");
  		// puts(str);
			if (k < face_property-1) {
        continue;
        //dummy = atoi(str);
			} else {
        v[k - face_property+1] = atoi(str);
      }

		}
		

	//	sscanf(cline, "%i %i %i %i", &dummy, &v[0], &v[1], &v[2]);
		for (int k = 0; k < 3; k++)	P[k] = (float)v[k] + 1.0;
		F[counter] = P;

	/*	if(counter%1000 == 0)
			printf("Face %i %f %f %f\n", counter, P[0] - 1.0, P[1] - 1.0, P[2] - 1.0);*/ 
		//	printf("face counter = %i\n",counter);

		counter++;
	}

	printf("faces read\n");

	if (counter < faces)
	{
		printf("unexpected end of file in ply file\n");
		exit(1);
	}

	label_non_manifold_faces(V, F, use_face, v_counter, f_counter);


	for (int i = 0; i < v_counter; i++)
	{
		tVertex v = createVertex(i + 1);
		v->point() = V[i];
		v->id() = i + 1;
		v->id2() = i + 1;
	}

	printf("verts set\n");


	for (int i = 0; i < f_counter; i++)
	{
		int v[3];
		v[0] = (int)F[i][0];// +1;
		v[1] = (int)F[i][1];// +1;
		v[2] = (int)F[i][2];// +1;
  
		if (use_face[i])
		{
			createFace(v, i + 1);
	  	}
	}

	printf("faces set\n");

	fclose(is);
	free(cline);

	//	printf("file closed\n");


	labelBoundaryEdges();
	//printf("boundary labeled\n");

	removeDanglingVertices();
	//	printf("dangling edges removed\n");

	reset_id2s();

	
	delete[] V;
	delete[] F;	 
	delete[] use_face;

}

void Solid::readOFF( const char * input)
{
	printf("reading off file...\n");
    FILE * is = fopen( input, "r" );
    if ( is == NULL )
    {
        fprintf(stderr,"Error is opening file %s\n", input );
        return;
    }

//    char line[MAX_LINE];
    int id, verts, faces;

	rewind(is);

	char * cline = (char*)calloc(1024,sizeof(char));
//	std::string sline;

	fgets(cline,1024,is);
	fgets(cline,1024,is);
	sscanf(cline,"%i %i",&verts,&faces);

	printf("%i verts, %i faces\n", verts, faces);

	int counter = 0;

    while ( !feof(is) && counter < verts)
    {
		fgets(cline,1024,is);
   		id = counter + 1;
        Point p;
		sscanf(cline,"%lf %lf %lf",&p[0],&p[1],&p[2]);
        tVertex v  = createVertex( id );
        v->point() = p;
        v->id()    = counter + 1;
		v->id2()    = counter + 1;
		counter++;
//		printf("on vertex %i\n",counter);
   }


	counter = 0;

    while ( !feof(is) && counter < faces)
    {
		fgets(cline,1024,is);
   		id = counter + 1;
		int v[3], dummy;
		//if(counter == 25682)
		//	printf("%s\n",cline);
		sscanf(cline,"%i %i %i %i",&dummy, &v[0],&v[1],&v[2]);
		v[0] += 1;
		v[1] += 1;
		v[2] += 1;

	//	printf("about to create face\n");
        createFace( v, id );
		counter++;
//		printf("on face %i\n",counter);
   }
	
	fclose(is);

    labelBoundaryEdges();
    removeDanglingVertices();
	reset_id2s();
}

void Solid::label_non_manifold_faces_test(Point *V, Point *F, bool * use_face, int v_counter, int f_counter, Point * test, int n_test_pts)
{

  printf("label_non_manifold_faces_test: What the hell is this?"
    "mem leaks everywhere, removed for now");
  exit(1);
  /*
	printf("v_counter, f_counter %i %i\n", v_counter, f_counter);
	int *record = new int[3*f_counter];
	for(int i = 0; i < 3*f_counter; i++) record[i] = i;
	
//	Point * V = m_V;
//	Point * F = m_F; 
//	bool * use_face = m_use_face;
	
	//Point *E = new Point[3*f_counter];//edges before halfedge, allows non-manifold edges, index v1, index v2, 0
	int *S = new int[3*f_counter];
	int *T = new int[3*f_counter];

	//bool * use_face = new bool[f_counter];
	//bool * boundary_face = new bool[f_counter];
	int * boundary_edges = new int[f_counter];
	for(int i = 0; i < f_counter; i++) {	use_face[i] = true; boundary_edges[i] = 0;}
		
	DList<int> * E_F = new DList<int>[3*f_counter]; //face id's for every edge
	DList<int> * F_E = new DList<int>[f_counter]; //edge id's for every face
	DList<int> * V_F = new DList<int>[v_counter]; //face id's for every vertex
//	DList<int> * V_E = new DList<int>[v_counter]; //edge id's for every vertex (not currently used)
	//int * E_F = new int*[3*f_counter]; //faces for every edge, allows non-manifold edges (stores face id)
	
	int * E_FN = new int[3*f_counter];//number of faces per edge (should be 2 if manifold)
	int * E_FNu = new int[3*f_counter];//number of remaining faces labeled as "use" per edge (should be 2 if output manifold)
	int * V_FN = new int[v_counter];//number of faces  per vertex 
//	int * V_EN = new int[v_counter];//number of edges  per vertex (not currently used)
	int * V_Fc = new int[v_counter];//number of face clusters per vertex (should be 1 if manifold)
	int * max_clusters = new int[f_counter];
	int e_counter = 0;

	for(int i = 0; i < v_counter; i++) {V_FN[i] = 0; V_Fc[i] = 0;}// V_EN[i] = 0;}
	for(int i = 0; i < f_counter; i++) {max_clusters[i] = 0;}

//	printf("past all allocs\n");



	AVL::Tree<EdgeKey> loc_edges;

	for(int i = 0; i < f_counter; i++)
	{
		int v[3];
		
		v[0] = (int)F[i][0];
		v[1] = (int)F[i][1];
		v[2] = (int)F[i][2];

		

		for(int j = 0; j < 3; j++)
		{
			if(v[j] == 0)
				printf("face %i, vert %i 0\n",i,j);

			V_F[v[j]-1].Append(&record[i]);
			V_FN[v[j]-1]++;
			for(int k = j+1; k < 3; k++)
			{
			
				EdgeKey loc_key(v[j],v[k]);
				EdgeKey * e = loc_edges.find( & loc_key );

				if(e == NULL)
				{					
					S[e_counter] = (int)std::min(v[j],v[k]);
					T[e_counter] = (int)std::max(v[j],v[k]);

				//	V_EN[v[j]-1]++;
				//	V_EN[v[k]-1]++;
				//	V_E[v[j]-1].Append(&record[e_counter]);
				//	V_E[v[k]-1].Append(&record[e_counter]);

					E_FN[e_counter] = 1;
					E_FNu[e_counter] = 1;
					E_F[e_counter].Append(&record[i]);
					e = new EdgeKey(v[j],v[k],e_counter);
					loc_edges.insert(e);
					F_E[i].Append(&record[e_counter]);
					e_counter++;
				}
				else
				{					
					int id_e = e->id();
					E_FN[id_e]++;
					E_FNu[id_e]++;
					E_F[id_e].Append(&record[i]);
					F_E[i].Append(&record[id_e]);
				}
			}
		}
	}

	printf("edges = %i\n",e_counter);


	for(int i = 0; i < f_counter; i++)
		use_face[i] = false;


	int * closest = new int[n_test_pts];
	double * min_dist   = new double[n_test_pts];


	for(int i = 0; i < n_test_pts; i++)
	{
		closest[i] = 0;
		min_dist[i] = 1e10;
	}
	printf("past closest alloc\n");
	
	for(int i = 0; i < v_counter; i++)
		for(int j = 0; j < n_test_pts; j++)
			if(min_dist[j] > (V[i] - test[j]).norm())
			{
				min_dist[j] = (V[i] - test[j]).norm();
				closest[j] = i;
			}

	printf("past closest set\n");

	//closest[0] = 0;

	int * cur_list_t = new int[f_counter+1];
	int * tmp_list_t = new int[f_counter+1];

	int radius = 3;
	n_test_pts = 1;

	bool * visited = new bool[f_counter];

	for(int i = 0; i < n_test_pts; i++)
	{
		int ind = closest[i];

		printf("ind = %i\n",ind);

		cur_list_t[0] = 0;
		printf("set cur list\n");
				 
		for( DListIterator<int> VFiter( V_F[ind] ); !VFiter.end(); ++ VFiter )
		{
			int * f_rec = *VFiter;
			cur_list_t[0]++;
			cur_list_t[cur_list_t[0]] = f_rec[0];
			printf("face %i added to init curlist\n",f_rec[0]);
		}

//		printf("first loop fileld\n");


		for(int j = 1; j <= cur_list_t[0]; j++)
			use_face[cur_list_t[j]] = true;

	//	printf("use face fileld\n");

		for(int j = 0; j < radius; j++)
		{
			for(int k = 0; k < f_counter; k++)
				visited[k] = false;

		//	printf("visited fileld\n");

			for(int k = 1; k <= cur_list_t[0]; k++)
				visited[cur_list_t[k]] = true;

		//	printf("visited fileld 2\n");


			tmp_list_t[0] = 0;

		//	printf("tmp_list_t 0 set\n");


			for(int k = 1; k <= cur_list_t[0]; k++)
				for( DListIterator<int> FEiter( F_E[cur_list_t[k]] ); !FEiter.end(); ++ FEiter )
				{
					int * e_rec = *FEiter;
					for( DListIterator<int> EFiter( E_F[e_rec[0]] ); !EFiter.end(); ++ EFiter )
					{
						int * f_rec2 = *EFiter;
						printf("on face %i\n",f_rec2[0]);
						
						if(!visited[f_rec2[0]])
						{
							tmp_list_t[0]++;
							tmp_list_t[tmp_list_t[0]] = f_rec2[0];
							visited[f_rec2[0]] = true;
						}

						printf("new face added\n");

					}
				}

						printf("about to set use_face\n");

			for(int k = 0; k <= tmp_list_t[0]; k++)
			{
				cur_list_t[k] = tmp_list_t[k];
				use_face[tmp_list_t[k]] = true;
			}	
			printf("set use_face\n");

		}

	}

	printf("thru neighborhood construction\n");
	
*/	
}

void Solid::label_non_manifold_faces(Point *V, Point *F, bool * use_face, int v_counter, int f_counter)
{

//	printf("v_counter, f_counter %i %i\n", v_counter, f_counter);
	int *record = new int[3*f_counter];
	for(int i = 0; i < 3*f_counter; i++) record[i] = i;
	
	Point *E = new Point[3*f_counter];//edges before halfedge, allows non-manifold edges, index v1, index v2, 0
	int *S = new int[3*f_counter];
	int *T = new int[3*f_counter];

	int * boundary_edges = new int[f_counter];
	for(int i = 0; i < f_counter; i++) {	use_face[i] = true; boundary_edges[i] = 0;}
		
	DList<int> * E_F = new DList<int>[3*f_counter]; //face id's for every edge
	DList<int> * F_E = new DList<int>[f_counter]; //edge id's for every face
	DList<int> * V_F = new DList<int>[v_counter]; //face id's for every vertex
//	DList<int> * V_E = new DList<int>[v_counter]; //edge id's for every vertex (not currently used)
	//int * E_F = new int*[3*f_counter]; //faces for every edge, allows non-manifold edges (stores face id)
	
	int * E_FN = new int[3*f_counter];//number of faces per edge (should be 2 if manifold)
	int * E_FNu = new int[3*f_counter];//number of remaining faces labeled as "use" per edge (should be 2 if output manifold)
	int * V_FN = new int[v_counter];//number of faces  per vertex 
//	int * V_EN = new int[v_counter];//number of edges  per vertex (not currently used)
	int * V_Fc = new int[v_counter];//number of face clusters per vertex (should be 1 if manifold)
	int * max_clusters = new int[f_counter];
	int e_counter = 0;

	for(int i = 0; i < v_counter; i++) {V_FN[i] = 0; V_Fc[i] = 0;}// V_EN[i] = 0;}
	for(int i = 0; i < f_counter; i++) {max_clusters[i] = 0;}

	AVL::Tree<EdgeKey> loc_edges;

	for(int i = 0; i < f_counter; i++)
	{
		int v[3];
		
		v[0] = (int)F[i][0];
		v[1] = (int)F[i][1];
		v[2] = (int)F[i][2];
		

		for(int j = 0; j < 3; j++)
		{
			V_F[v[j]-1].Append(&record[i]);
			V_FN[v[j]-1]++;
			for(int k = j+1; k < 3; k++)
			{
			
				EdgeKey loc_key(v[j],v[k]);
				EdgeKey * e = loc_edges.find( & loc_key );

				if(e == NULL)
				{					
					S[e_counter] = (int)std::min(v[j],v[k]);
					T[e_counter] = (int)std::max(v[j],v[k]);

				//	V_EN[v[j]-1]++;
				//	V_EN[v[k]-1]++;
				//	V_E[v[j]-1].Append(&record[e_counter]);
				//	V_E[v[k]-1].Append(&record[e_counter]);

					E_FN[e_counter] = 1;
					E_FNu[e_counter] = 1;
					E_F[e_counter].Append(&record[i]);
					e = new EdgeKey(v[j],v[k],e_counter);
					//EdgeKey new_e = EdgeKey(v[j],v[k],e_counter);
					//loc_edges.insert(&new_e);
					loc_edges.insert(e);
					F_E[i].Append(&record[e_counter]);
					e_counter++;
				}
				else
				{					
					int id_e = e->id();
					E_FN[id_e]++;
					E_FNu[id_e]++;
					E_F[id_e].Append(&record[i]);
					F_E[i].Append(&record[id_e]);
				}
			}
		}
	}

	//need to delete temp edges before loc_edges is destroyed
	for ( AVL::TreeIterator<EdgeKey> leiter(loc_edges); !leiter.end(); ++leiter )
    {
        EdgeKey * e = *leiter;
		if(e != NULL)
			delete  e;
	}

//	printf("edges = %i\n",e_counter);

	for(int i = 0; i < e_counter; i++)
	{
		if(E_FNu[i] > 2)
		{
		//	printf("edge %i has more than two faces\n",i);
		
			int * tmp_faces = new int[E_FN[i]];
			int loc_faces = 0;
			for( DListIterator<int> iter( E_F[i] ); !iter.end(); ++ iter )
			{
				int * rec = *iter;

				if(use_face[rec[0]])
				{
				//	printf("id face = %i\n",rec[0]);
					tmp_faces[loc_faces] = rec[0];
					loc_faces++;
				}
			}
			E_FNu[i] = loc_faces;


			for(int j = 0; j < loc_faces; j++)
				for(int k = j+1; k < loc_faces; k++)
					if(point_match(F[tmp_faces[j]],F[tmp_faces[k]]) && use_face[tmp_faces[j]])
					{
						use_face[tmp_faces[j]] = false;
						for( DListIterator<int> f_iter( F_E[tmp_faces[j]] ); !f_iter.end(); ++ f_iter )
						{
							int * f_rec = *f_iter;
							E_FNu[f_rec[0]]--;
						}
						if(j%10000 == 0 || i%10000 == 0)
							printf("duplicate %i %i\n",i,j);
					}

			delete [] tmp_faces;
		}
	}

	//find number of boundary edges for each face
	for(int i = 0; i < f_counter; i++)
	{
		//int boundary_edges = 0;
		if(use_face[i])//save some time
		for( DListIterator<int> iter( F_E[i] ); !iter.end(); ++ iter )
		{
			int * rec = *iter;

			if(E_FNu[rec[0]] <= 1)
				boundary_edges[i]++;
		}
		if(boundary_edges[i] == 3)//erase isolated faces
		{
		//	printf("isolated face %i\n",i);
			for( DListIterator<int> iter( F_E[i] ); !iter.end(); ++ iter )
			{
				int * rec = *iter;
				E_FNu[rec[0]]--;
			}
			use_face[i] = false;
		}
	}

	//get rid of face pairs which share more than one edge
	for(int i = 0; i < f_counter; i++)
	{
		bool degenerate = false;
		if(use_face[i])
			for( DListIterator<int> FEiter1( F_E[i] ); !FEiter1.end(); ++ FEiter1 )
			{
				int * e_rec1 = *FEiter1;
				for( DListIterator<int> EFiter( E_F[e_rec1[0]] ); !EFiter.end(); ++ EFiter )
				{
					int *f_rec = *EFiter;
					if(f_rec[0] != i && use_face[f_rec[0]])//different face
					{
						for( DListIterator<int> FEiter2( F_E[f_rec[0]] ); !FEiter2.end(); ++ FEiter2 )
						{
							int * e_rec2 = *FEiter2;
							if(e_rec2[0] != e_rec1[0] /*different edge*/ && F_E[i].contains(&record[e_rec2[0]]))
							{
								degenerate = true;
								use_face[i] = false;
								use_face[f_rec[0]] = false;
			//					printf("degenerate faces %i %i\n",i,f_rec[0]);

								break;
							}
						}
					}
					if(degenerate)
					{
						for( DListIterator<int> iter2( F_E[f_rec[0]] ); !iter2.end(); ++ iter2 )
						{
							int * rec2 = *iter2;
							E_FNu[rec2[0]]--;
						}
					}
				}
			}
	}
	//return;

	int ** V_Fclabel = new int*[v_counter];
	for(int i = 0; i < v_counter; i++)//find connected clusters of faces about a vertex
	{
		V_Fc[i] = 0;
		V_Fclabel[i] = new int[V_FN[i]];

		for(int j = 0; j < V_FN[i]; j++)
			V_Fclabel[i][j] = 0;

		int *cur_list = new int[V_FN[i] + 1];
				
		int j = 0;
		
		for(DListIterator<int> VFiter( V_F[i] ); !VFiter.end(); ++ VFiter )		
		{
			int *f_rec = *VFiter;//face belonging to vertex

			if(V_Fclabel[i][j] == 0 && use_face[f_rec[0]] )
			{
				cur_list[0] = 1;
				cur_list[cur_list[0]] = f_rec[0];

				for(int n = 1; n <= cur_list[0]; n++)
				{
						
					int cur_F = cur_list[n];

					for( DListIterator<int> FEiter( F_E[cur_F] ); !FEiter.end(); ++ FEiter )//populating new cluster
					{
						int * e_rec = *FEiter;//face belogning to edge
						if(S[e_rec[0]] == i+1 || T[e_rec[0]] == i+1)//edge contains current vertex	
						for( DListIterator<int> EFiter( E_F[e_rec[0]] ); !EFiter.end(); ++ EFiter )
						{
							int *f_rec2 = *EFiter;
							bool contains  = false;

							if(use_face[f_rec2[0]])
							for(int k = cur_list[0]; k >= 1; k--)//checking if current face is in current cluster already
								if(f_rec2[0] == cur_list[k])
								{
									contains = true;
									break;
								}

							if(!contains && use_face[f_rec2[0]])
							{
								cur_list[0]++;
								cur_list[cur_list[0]] = f_rec2[0];	
							}						
						}
					}
				}

				V_Fc[i]++;			
				
				for(int n = 1; n <= cur_list[0]; n++)
				{
					int cur_F = cur_list[n];
					int k = 0;

					for(DListIterator<int> VFiter2( V_F[i] ); !VFiter2.end(); ++ VFiter2 )//assigning cluster index to face-vert pair
					{
						int *f_rec3 = *VFiter2;
						if(cur_F == f_rec3[0])
						{
							V_Fclabel[i][k] = V_Fc[i];
							break;
						}
						k++;
					}
				}
			}
			j++;
		}

		delete []  cur_list;

//		if(V_Fc[i] > 1)
//			printf("vert %i has %i clusters\n",i, V_Fc[i]);


//		if(V_Fc[i] == 0) 
//			printf("vert %i has %i clusters\n",i, V_Fc[i]);
		

		//now reset max clusters for each face
		for(DListIterator<int> VFiter( V_F[i] ); !VFiter.end(); ++ VFiter )
		{
			int *f_rec = *VFiter;//face belonging to vertex

			if(max_clusters[f_rec[0]] < V_Fc[i])
				max_clusters[f_rec[0]] = V_Fc[i];
		}

	}//END VERTEX-FACE CLUSTER LOOP


	for(int i = 0; i < e_counter; i++)
	{
		if(E_FNu[i] > 2)
		{
	//		printf("edge %i still has more than two faces\n",i);
		
			int * tmp_faces = new int[E_FNu[i]];
			int loc_faces = 0;
			for( DListIterator<int> iter( E_F[i] ); !iter.end(); ++ iter )
			{
				int * rec = *iter;

				if(use_face[rec[0]])
				{
		//			printf("id face = %i\n",rec[0]);
					tmp_faces[loc_faces] = rec[0];
					loc_faces++;
				}
			}
			E_FNu[i] = loc_faces;


			
			for(int j = 0; j < loc_faces; j++)
				if(max_clusters[tmp_faces[j]] > 1 && use_face[tmp_faces[j]] )
				{
					bool cluster_face = false;
					for(int k = 0; k < 3; k++)
						if((int)F[tmp_faces[j]][k] != S[i] 
							&& (int)F[tmp_faces[j]][k] != T[i] 
							&& V_Fc[(int)F[tmp_faces[j]][k]-1] > 1)//make sure it's the vertex not in this edge
						{cluster_face = true; /*printf("cluster vert %i\n",(int)F[tmp_faces[j]][k]);*/}

					if(cluster_face)
					{
						use_face[tmp_faces[j]] = false;
						for( DListIterator<int> iter2( F_E[tmp_faces[j]] ); !iter2.end(); ++ iter2 )
						{
							int * rec2 = *iter2;
							E_FNu[rec2[0]]--;
						}
					//	printf("cluster face\n");
					}
					
				}

			delete [] tmp_faces;
		}
	}


	//A round of preferential face deletion, only taking out faces with 3 NM edges

	int * num_NM_edges = new int[f_counter];

	for(int i = 0; i < f_counter; i++)
		num_NM_edges[i] = 0;


//	printf("about to do preferential face deletion\n");

	for(int i = 0; i < e_counter; i++)
	{
		if(E_FNu[i] > 2)
		{		
			int * tmp_faces = new int[E_FNu[i]];
			int loc_faces = 0;
			for( DListIterator<int> iter( E_F[i] ); !iter.end(); ++ iter )
			{
				int * rec = *iter;

				if(use_face[rec[0]])
				{
					tmp_faces[loc_faces] = rec[0];
					loc_faces++;
				}
			}
			E_FNu[i] = loc_faces;

			if(E_FNu[i] > 2)
				for(int j = 0; j < E_FNu[i]; j++)
					num_NM_edges[tmp_faces[j]]++;

			delete [] tmp_faces;
		}
	}


//	printf("done setting NM_edge number\n");


	if(1)
	for(int i = 0; i < e_counter; i++)
	{
		if(E_FNu[i] > 2)
		{
		//	printf("edge %i still has more than two faces\n",i);
			
			for( DListIterator<int> EFiter( E_F[i] ); !EFiter.end(); ++ EFiter )
			{
				int * f_rec = *EFiter;

				if(use_face[f_rec[0]] && num_NM_edges[f_rec[0]] > 2)
				{
					use_face[f_rec[0]] = false;
					for( DListIterator<int> EFiter( F_E[f_rec[0]] ); !EFiter.end(); ++ EFiter )
					{
						int * e_rec = *EFiter;
						E_FNu[e_rec[0]]--;
						for( DListIterator<int> EFiter2( E_F[e_rec[0]] ); !EFiter2.end(); ++ EFiter2 )
						{
							int * f_rec2 = * EFiter2;
							if(num_NM_edges[f_rec2[0]] > 0)
								num_NM_edges[f_rec2[0]]--;
						}

					}
				//	printf("face %i has 3 NM edges\n", f_rec[0]);
				}
			}			
		}
	}


	//One more round of preferential face deletion, taking out faces with at least 2 NM edges
	for(int i = 0; i < f_counter; i++)
		num_NM_edges[i] = 0;

	for(int i = 0; i < e_counter; i++)
	{
		if(E_FNu[i] > 2)
		{
			int * tmp_faces = new int[E_FNu[i]];
			int loc_faces = 0;
			for( DListIterator<int> iter( E_F[i] ); !iter.end(); ++ iter )
			{
				int * rec = *iter;

				if(use_face[rec[0]])
				{
					tmp_faces[loc_faces] = rec[0];
					loc_faces++;
				}
			}
			E_FNu[i] = loc_faces;

			if(E_FNu[i] > 2)
				for(int j = 0; j < E_FNu[i]; j++)
					num_NM_edges[tmp_faces[j]]++;

			delete [] tmp_faces;
		}
	}

	if(1)
	for(int i = 0; i < e_counter; i++)
	{
		if(E_FNu[i] > 2)
		{
		//	printf("edge %i still has more than two faces\n",i);
			
			for( DListIterator<int> EFiter( E_F[i] ); !EFiter.end(); ++ EFiter )
			{
				int * f_rec = *EFiter;

				if(use_face[f_rec[0]] && num_NM_edges[f_rec[0]] > 1)
				{
					use_face[f_rec[0]] = false;
					for( DListIterator<int> EFiter( F_E[f_rec[0]] ); !EFiter.end(); ++ EFiter )
					{
						int * e_rec = *EFiter;
						E_FNu[e_rec[0]]--;
						for( DListIterator<int> EFiter2( E_F[e_rec[0]] ); !EFiter2.end(); ++ EFiter2 )
						{
							int * f_rec2 = * EFiter2;
							if(num_NM_edges[f_rec2[0]] > 0)
								num_NM_edges[f_rec2[0]]--;
						}

					}
				//	printf("face %i has 2 NM edges\n", f_rec[0]);
				}
			}			
		}
	}

	delete [] num_NM_edges;




	for(int i = 0; i < e_counter; i++)
	{
		if(E_FNu[i] > 2)
		{
		//	printf("edge %i still has more than two faces\n",i);
		
			int * tmp_faces = new int[E_FNu[i]];
			int loc_faces = 0;
			for( DListIterator<int> iter( E_F[i] ); !iter.end(); ++ iter )
			{
				int * rec = *iter;

				if(use_face[rec[0]])
				{
			//		printf("id face = %i\n",rec[0]);
					tmp_faces[loc_faces] = rec[0];
					loc_faces++;
				}
			}
			E_FNu[i] = loc_faces;


			if(E_FNu[i] > 2)
			for(int j = 0; j < loc_faces; j++)
				if(boundary_edges[tmp_faces[j]] > 1 && use_face[tmp_faces[j]] && E_FNu[i] > 2)
				{
					use_face[tmp_faces[j]] = false;
					E_FNu[i]--;
			//		printf("2 boundary face\n");
				}
			delete [] tmp_faces;
		}
	}

	for(int i = 0; i < e_counter; i++)
	{
		if(E_FNu[i] > 2)
		//if(E_FN[i] > 2)
		{
		//	printf("edge %i still has more than two faces\n",i);
		
			int * tmp_faces = new int[E_FNu[i]];
			//int * tmp_faces = new int[E_FN[i]];
			int loc_faces = 0;
			for( DListIterator<int> iter( E_F[i] ); !iter.end(); ++ iter )
			{
				int * rec = *iter;

				if(use_face[rec[0]])
				{
					use_face[rec[0]] = false; // get rid of ALL bad edges
				//	printf("id face = %i\n",rec[0]);
					tmp_faces[loc_faces] = rec[0];
					loc_faces++;
				}
			}
			E_FNu[i] = loc_faces;


			if(E_FNu[i] > 2)
			for(int j = 0; j < loc_faces; j++)
				if(boundary_edges[tmp_faces[j]] > 0 && use_face[tmp_faces[j]] && E_FNu[i] > 2)
				{
					use_face[tmp_faces[j]] = false;
					E_FNu[i]--;
				//	printf("1 boundary face\n");
				}
			delete [] tmp_faces;
		}
	}

	//return;

	for(int i = 0; i < e_counter; i++)
	{
		if(E_FNu[i] > 2)
		{
		//	printf("edge %i still has more than two faces\n",i);
		
			int * tmp_faces = new int[E_FNu[i]];
			int loc_faces = 0;
			for( DListIterator<int> iter( E_F[i] ); !iter.end(); ++ iter )
			{
				int * rec = *iter;

				if(use_face[rec[0]])
				{
				//	printf("id face = %i\n",rec[0]);
					tmp_faces[loc_faces] = rec[0];
					loc_faces++;
				}
			}
			E_FNu[i] = loc_faces;


			double * cotan_sum = new double[loc_faces];
			int first_ind = -1, second_ind = -1; 

			double first_cotsum = -1;
			double second_cotsum = -1;

			if(E_FNu[i] > 2)
			{
				int k = 0;
				for(int j = 0; j < loc_faces; j++)//keep 2 faces with 2 lowest |cot| sums, i.e. best-formed triangles
					if(use_face[tmp_faces[j]])
					{
						cotan_sum[j] = cotsum(V[(int)F[tmp_faces[j]][0]],V[(int)F[tmp_faces[j]][1]],V[(int)F[tmp_faces[j]][2]]);
							//cotsum to be defined in common routines


				//		printf("cotsum of face %i %lf\n",tmp_faces[j], cotan_sum[j]);

						if(k == 0)
						{
							first_cotsum = cotan_sum[j];
							first_ind = j;
						}
						if(k == 1)
						{
							if(first_cotsum > cotan_sum[j])
							{								
								second_cotsum = first_cotsum;
								first_cotsum = cotan_sum[j];

								second_ind = first_ind;
								first_ind = j;
							}
							else
							{
								second_cotsum = cotan_sum[j];
								second_ind = j;
							}
						}
						else
						{
							if(first_cotsum > cotan_sum[j])
							{								
								second_cotsum = first_cotsum;
								first_cotsum = cotan_sum[j];

								second_ind = first_ind;
								first_ind = j;
							}
							else if(second_cotsum > cotan_sum[j])
							{
								second_cotsum = cotan_sum[j];
								second_ind = j;
							}						
						}
						k++;
					
					}
			}
			else 
			{
	//			cotan_sum[0] = cotsum(V[(int)F[tmp_faces[0]][0]],V[(int)F[tmp_faces[0]][1]],V[(int)F[tmp_faces[0]][2]]);
	//			cotan_sum[1] = cotsum(V[(int)F[tmp_faces[1]][0]],V[(int)F[tmp_faces[1]][1]],V[(int)F[tmp_faces[1]][2]]);
				first_ind = 0;
				second_ind = 1;
			}

			delete [] cotan_sum;//new 8-27-14

	//		printf("using faces %i %i\n", tmp_faces[first_ind], tmp_faces[second_ind]);
	//		printf("cotsums %lf %lf\n", cotan_sum[first_ind], cotan_sum[second_ind]);

			for(int j = 0; j < loc_faces; j++)
				if(j != first_ind && j != second_ind)
					use_face[tmp_faces[j]] = false;

			
			delete [] tmp_faces;

		}		
	}


	int max_total_clusters = 2;
	int declust_round = 0;
	
	//if(0)
	while(max_total_clusters > 1)
	{
		max_total_clusters = 1;
		declust_round++;
	//	printf("declustering round %i\n",declust_round);
	for(int i = 0; i < v_counter; i++)//find connected clusters of faces about a vertex
	{
		V_Fc[i] = 0;
		//V_Fclabel[i] = new int[V_FN[i]];
		int max_clust_count = 0;
		int max_clust_id = 1;

		for(int j = 0; j < V_FN[i]; j++)
			V_Fclabel[i][j] = 0;

		int *cur_list = new int[V_FN[i] + 1];
				
		int j = 0;

		for(DListIterator<int> VFiter( V_F[i] ); !VFiter.end(); ++ VFiter )		
		{
			int *f_rec = *VFiter;//face belonging to vertex


			if(V_Fclabel[i][j] == 0 && use_face[f_rec[0]] )
			{
				cur_list[0] = 1;
				cur_list[cur_list[0]] = f_rec[0];

				for(int n = 1; n <= cur_list[0]; n++)
				{
					int cur_F = cur_list[n];

					for( DListIterator<int> FEiter( F_E[cur_F] ); !FEiter.end(); ++ FEiter )//populating new cluster
					{
						int * e_rec = *FEiter;//face belogning to edge
							
		
						if(S[e_rec[0]] == i+1 || T[e_rec[0]] == i+1)//edge contains current vertex
						for( DListIterator<int> EFiter( E_F[e_rec[0]] ); !EFiter.end(); ++ EFiter )
						{
							int *f_rec2 = *EFiter;
							bool contains  = false;

							if(use_face[f_rec2[0]])
							for(int k = cur_list[0]; k >= 1; k--)//checking if current face is in current cluster already
								if(f_rec2[0] == cur_list[k])
								{
									contains = true;
									break;
								}

							if(!contains && use_face[f_rec2[0]])
							{
								cur_list[0]++;
								cur_list[cur_list[0]] = f_rec2[0];							
							}
						}
					}
				}

				V_Fc[i]++;

				if(max_total_clusters < V_Fc[i])
					max_total_clusters = V_Fc[i];

				if(V_Fc[i] == 1)
					max_clust_count = cur_list[0];
				else if(max_clust_count < cur_list[0])
				{
					max_clust_count = cur_list[0];
					max_clust_id = V_Fc[i];
				}
				
				for(int n = 1; n <= cur_list[0]; n++)
				{
					int cur_F = cur_list[n];
					int k = 0;

					for(DListIterator<int> VFiter2( V_F[i] ); !VFiter2.end(); ++ VFiter2 )//assigning cluster index to face-vert pair
					{
						int *f_rec3 = *VFiter2;
						if(cur_F == f_rec3[0])
						{
							V_Fclabel[i][k] = V_Fc[i];
							break;
						}
						k++;
					}
				}
			}
			j++;
		} //END VFITER1

	
	//	if(V_FN[i] > 0)
		delete []  cur_list;
	//	else
	//		cur_list = NULL;

//		if(V_Fc[i] > 1)
//			printf("vert %i has %i clusters\n",i, V_Fc[i]);

		//	if(V_Fc[i] == 0) 
	//		printf("vert %i has %i clusters\n",i, V_Fc[i]);		

		//now set all but biggest cluster to !use
		int ii = 0;
		if(V_Fc[i] > 0)
		for(DListIterator<int> VFiter( V_F[i] ); !VFiter.end(); ++ VFiter )
		{
			int *f_rec = *VFiter;//face belonging to vertex

			if(V_Fclabel[i][ii] != max_clust_id)
				use_face[f_rec[0]] = false;
			ii++;
		}
	
	}
	}//END WHILE FOR VERTEX-FACE CLUSTER



	int *cur_list_he = new int[f_counter + 1];
	bool *visited = new bool[f_counter];

//	if(1)//for(int kk = 0; kk < 3; kk++)
//	{

	int to_use = 0, num_visited = 0;
	for(int i = 0; i < f_counter; i++ )		
	{
		visited[i] = false;
		if(use_face[i]) to_use++;
	}	

//	printf("fixing halfedge direction\n");

	//if(0)
	//while(num_visited < to_use)
	
	for(int I = 0; I < f_counter; I++)		
		if(use_face[I] && !visited[I])
		{
			if(I != 0)
			//	printf("in next component\n");

			visited[I] = true;
			num_visited++;
			cur_list_he[0] = 1;
			cur_list_he[cur_list_he[0]] = I;

			for(int n = 1; n <= cur_list_he[0]; n++)
			{						
				int cur_F = cur_list_he[n];
				visited[cur_F] = true;
				
				for( DListIterator<int> FEiter( F_E[cur_F] ); !FEiter.end(); ++ FEiter )//populating new cluster
				{
					int * e_rec = *FEiter;//face belogning to edge
						
					for( DListIterator<int> EFiter( E_F[e_rec[0]] ); !EFiter.end(); ++ EFiter )
					{
						int *f_rec2 = *EFiter;
						
						if(use_face[f_rec2[0]] && !visited[f_rec2[0]])
						{
							cur_list_he[0]++;
							cur_list_he[cur_list_he[0]] = f_rec2[0];

						//	printf("cur_F, neighbor face %i %i\n", cur_F, f_rec2[0]);

							//this is where you check and, if needed, fix halfedge direction
							int u1 = -1919, u2 = -1919, v1 = -1919, v2 = -1919;

							for(int k = 0; k < 3; k++)
							{
								if((int)F[cur_F][k] == S[e_rec[0]]) u1 = k;
								if((int)F[cur_F][k] == T[e_rec[0]]) u2 = k;
								if((int)F[f_rec2[0]][k] == S[e_rec[0]]) v1 = k;
								if((int)F[f_rec2[0]][k] == T[e_rec[0]]) v2 = k;
							}

              //u1, u2, v1, v2 need a init check here
              if(u1 == -1919 || u2 == -1919 || v1 == -1919 || v2 == -1919){
                printf("[solid/label_non_manifold_faces]:"
                  "initialization error.");
                exit(1);
              }

							if( (3+u1-u2)%3 == (3+v1-v2)%3 )//reverse order if needed
							{
							//	printf("reversing face %i\n", f_rec2[0]);
							//	F[f_rec2[0]].print();
								
								F[f_rec2[0]][v1] = (float)T[e_rec[0]];
								F[f_rec2[0]][v2] = (float)S[e_rec[0]];
							//	printf("now,\n");
							//	F[f_rec2[0]].print();
							//	use_face[f_rec2[0]] = false;
							}

							visited[f_rec2[0]] = true;
							num_visited++;
						}						
					}
				}
			}
		}

	
		//printf("done halfedge direction\n\n");

		bool moebius_loops = false;

		//if(1)
		for(int n = 0; n < f_counter; n++)
		{						
				int cur_F = n;
				//visited[cur_F] = true;
				
				if(use_face[cur_F])
				for( DListIterator<int> FEiter( F_E[cur_F] ); !FEiter.end(); ++ FEiter )//populating new cluster
				{
					int * e_rec = *FEiter;//face belogning to edge
						
					for( DListIterator<int> EFiter( E_F[e_rec[0]] ); !EFiter.end(); ++ EFiter )
					{
						int *f_rec2 = *EFiter;
						
						if(use_face[f_rec2[0]] && f_rec2[0] != cur_F)
						{
							//this is where you check and, if needed, fix halfedge direction

							int u1 = -1919, u2 = -1919, v1 = -1919, v2 = -1919;

							for(int k = 0; k < 3; k++)
							{
								if((int)F[cur_F][k] == S[e_rec[0]]) u1 = k;
								if((int)F[cur_F][k] == T[e_rec[0]]) u2 = k;
								if((int)F[f_rec2[0]][k] == S[e_rec[0]]) v1 = k;
								if((int)F[f_rec2[0]][k] == T[e_rec[0]]) v2 = k;
							}

              //u1, u2, v1, v2 need a init check here
              if(u1 == -1919 || u2 == -1919 || v1 == -1919 || v2 == -1919){
                printf("[solid/label_non_manifold_faces]:"
                  "initialization error.");
                exit(1);
              }

							if( (3+u1-u2)%3 == (3+v1-v2)%3 )//reverse order if needed
							{
							//	printf("conflicting faces %i %i\n",cur_F,f_rec2[0]);
								F[cur_F].print();
								F[f_rec2[0]].print();
								use_face[f_rec2[0]] = false;
								use_face[cur_F] = false;
								moebius_loops = true;
							}
							
						}						
					}
				}
		}



	if(moebius_loops)//need to fix non-manifold vertices that may have been created by face deletion above
	{
		max_total_clusters = 2;
		declust_round = 0;
	
		while(max_total_clusters > 1)
		{
			max_total_clusters = 1;
			declust_round++;
			//printf("declustering round %i\n",declust_round);
			for(int i = 0; i < v_counter; i++)//find connected clusters of faces about a vertex
			{
				V_Fc[i] = 0;
				//V_Fclabel[i] = new int[V_FN[i]];
				int max_clust_count = 0;
				int max_clust_id = 1;

				for(int j = 0; j < V_FN[i]; j++)
					V_Fclabel[i][j] = 0;

				int *cur_list = new int[V_FN[i] + 1];
				
				int j = 0;

				for(DListIterator<int> VFiter( V_F[i] ); !VFiter.end(); ++ VFiter )		
				{
					int *f_rec = *VFiter;//face belonging to vertex
			
					if(V_Fclabel[i][j] == 0 && use_face[f_rec[0]] )
					{
						cur_list[0] = 1;
						cur_list[cur_list[0]] = f_rec[0];

						for(int n = 1; n <= cur_list[0]; n++)
						{
							int cur_F = cur_list[n];

							for( DListIterator<int> FEiter( F_E[cur_F] ); !FEiter.end(); ++ FEiter )//populating new cluster
							{
								int * e_rec = *FEiter;//face belogning to edge
							
		
								if(S[e_rec[0]] == i+1 || T[e_rec[0]] == i+1)//edge contains current vertex
								for( DListIterator<int> EFiter( E_F[e_rec[0]] ); !EFiter.end(); ++ EFiter )
								{
									int *f_rec2 = *EFiter;
									bool contains  = false;

									if(use_face[f_rec2[0]])
									for(int k = cur_list[0]; k >= 1; k--)//checking if current face is in current cluster already
										if(f_rec2[0] == cur_list[k])
										{
											contains = true;
											break;
										}

									if(!contains && use_face[f_rec2[0]])
									{
										cur_list[0]++;
										cur_list[cur_list[0]] = f_rec2[0];							
									}
								}
							}
						}

						V_Fc[i]++;

						if(max_total_clusters < V_Fc[i])
							max_total_clusters = V_Fc[i];

						if(V_Fc[i] == 1)
							max_clust_count = cur_list[0];
						else if(max_clust_count < cur_list[0])
						{
							max_clust_count = cur_list[0];
							max_clust_id = V_Fc[i];
						}
				
						for(int n = 1; n <= cur_list[0]; n++)
						{
							int cur_F = cur_list[n];
							int k = 0;

							for(DListIterator<int> VFiter2( V_F[i] ); !VFiter2.end(); ++ VFiter2 )//assigning cluster index to face-vert pair
							{
								int *f_rec3 = *VFiter2;
								if(cur_F == f_rec3[0])
								{
									V_Fclabel[i][k] = V_Fc[i];
									break;
								}
								k++;
							}
						}
					}
					j++;
				} //END VFITER1

	
				if(V_FN[i] > 1)//don't know why, but otherwise it's a memory error
					delete []  cur_list;
				else
					cur_list = NULL;

			//	if(V_Fc[i] > 1)
			//		printf("vert %i has %i clusters\n",i, V_Fc[i]);

		//	if(V_Fc[i] == 0) 
	//		printf("vert %i has %i clusters\n",i, V_Fc[i]);		

		//now set all but biggest cluster to !use
				int ii = 0;
				if(V_Fc[i] > 0)
				for(DListIterator<int> VFiter( V_F[i] ); !VFiter.end(); ++ VFiter )
				{
					int *f_rec = *VFiter;//face belonging to vertex

					if(V_Fclabel[i][ii] != max_clust_id)
						use_face[f_rec[0]] = false;
					ii++;
				}
	
			}
		}//END WHILE FOR VERTEX-FACE CLUSTER 3
	}
	
	delete [] cur_list_he;
	delete [] visited;	


	delete [] record;
	delete [] E ;
	delete [] S;
	delete [] T;

	delete [] boundary_edges;
		
	delete [] E_F;
	delete [] F_E;
	delete [] V_F;
//	delete [] V_E;
	
	delete [] E_FN;
	delete [] E_FNu;
	delete [] V_FN;
//	delete [] V_EN;
	delete [] V_Fc;
	delete [] max_clusters;

	for(int i = 0; i < v_counter; i++)
		delete [] V_Fclabel[i];

	delete [] V_Fclabel;

}

void Solid::close_std_shape(int bw, int cur_verts, int cur_faces)
{
//	printf("bw = %i\n", bw);

	Point south;
	Point north;
	int i;

	for(i = 1; i <= 2*bw; i++)
		south += idVertex(i + cur_verts)->point();

	south /= 2*bw;


//	printf("south done\n");
//	south.print();

	for(i = 2*bw*(2*bw-1)+1; i <= 4*bw*bw; i++)
		north += idVertex(i + cur_verts)->point();

	north /= 2*bw;
//	north.print();
//	exit(0);


//	printf("north done\n");


	tVertex v_south  = createVertex( 4*bw*bw + 1 + cur_verts);
    v_south->point() = south;
    v_south->id()    = 4*bw*bw + 1 + cur_verts;
	v_south->id2()   = 4*bw*bw + 1 + cur_verts;

//	printf("south created\n");


	tVertex v_north  = createVertex( 4*bw*bw + 2 + cur_verts);
    v_north->point() = north;
    v_north->id()    = 4*bw*bw + 2 + cur_verts;
	v_north->id2()   = 4*bw*bw + 2 + cur_verts;

//	printf("north created\n");


	int n_faces_open = numFaces();
	int v_i[3];

	for(i = 0; i < 2*bw; i++) //south
	{
		//printf("on i = %i\n",i);
		
		v_i[0] = v_south->id();
		v_i[2] = i%(2*bw) + 1 + cur_verts;
		v_i[1] = (i+1)%(2*bw) + 1 + cur_verts;

        createFace( v_i, n_faces_open + i + 1 + cur_faces );
//		printf("f->id() = %i\n",f->id());
		
	}

//	printf("south faces\n");



	for(i = 0; i < 2*bw; i++) //north
	{

		int ii = 2*bw*(2*bw - 1);
	//	printf("on i = %i, ii = %i\n",i,ii);

		v_i[0] = v_north->id();
		v_i[1] = ii + i%(2*bw) + 1 + cur_verts;
		v_i[2] = ii + (i+1)%(2*bw) + 1 + cur_verts;

        createFace( v_i, n_faces_open + i + 1 + 2*bw + cur_faces );
	}

//	printf("north faces\n");



}

void Solid::reset_regular_poles(int bw)
{
	Point Sval;
	Point Nval;	
	for(int i = 0; i < 2*bw; i++)
		Sval += idVertex(i+1)->point();
	Sval = Sval/(double)(2*bw);

	for(int i = 4*bw*bw - 2*bw; i < 4*bw*bw; i++)
		Nval += idVertex(i+1)->point();
	Nval = Nval/(double)(2*bw);

	idVertex(4*bw*bw)->point() = Sval;
	idVertex(4*bw*bw+1)->point() = Nval;
}




//Euler operation

Solid::tHalfEdge Solid::vertexMostClwOutHalfEdge( Solid::tVertex  v )
{
    return v->most_clw_out_halfedge();
};

Solid::tHalfEdge Solid::vertexMostCcwOutHalfEdge( Solid::tVertex  v )
{
    return v->most_ccw_out_halfedge();
};

Solid::tHalfEdge  Solid::corner( tVertex v, tFace f)
{
    Solid::tHalfEdge he = f->halfedge();
    do
    {
        if ( he->vertex() == v )
            return he;
        he = he->he_next();
    }while ( he != f->halfedge() );
    return NULL;
};

Solid::tHalfEdge Solid::vertexNextCcwOutHalfEdge( Solid::tHalfEdge  he )
{
    return he->ccw_rotate_about_source();
};

Solid::tHalfEdge Solid::vertexNextClwOutHalfEdge( Solid::tHalfEdge  he )
{
    if ( he->he_sym() == NULL ) throw TopologyException();
    return he->clw_rotate_about_source();
};

Solid::tHalfEdge Solid::vertexMostClwInHalfEdge( Solid::tVertex  v )
{
    return v->most_clw_in_halfedge();
};

Solid::tHalfEdge Solid::vertexMostCcwInHalfEdge( Solid::tVertex  v )
{
    return v->most_ccw_in_halfedge();
};

Solid::tHalfEdge Solid::vertexNextCcwInHalfEdge( Solid::tHalfEdge  he )
{
    if ( he->he_sym() == NULL ) throw TopologyException();;
    return he->ccw_rotate_about_target();
};

Solid::tHalfEdge Solid::vertexNextClwInHalfEdge( Solid::tHalfEdge  he )
{
    return he->clw_rotate_about_target();
};

Solid::tHalfEdge Solid::faceNextClwHalfEdge( Solid::tHalfEdge  he )
{
    return he->he_prev();
};

Solid::tHalfEdge Solid::faceNextCcwHalfEdge( Solid::tHalfEdge  he )
{
    return he->he_next();
};

Solid::tHalfEdge Solid::faceMostCcwHalfEdge( Solid::tFace  face )
{
    return face->halfedge();
};

Solid::tHalfEdge Solid::faceMostClwHalfEdge( Solid::tFace  face )
{
    return face->halfedge()->he_next();
};


//access id->v
Solid::tVertex Solid::idVertex( int id ) 
{
    Vertex v;
    v.id() = id;
    return m_verts.find( &v );
};

//access v->id
int Solid::vertexId( Solid::tVertex  v )
{
    return v->id();
};

//access id->f
Solid::tFace Solid::idFace( int id )
{
    Face f;
    f.id() = id;

    return m_faces.find( &f );
};

//acess f->id
int Solid::faceId( Solid::tFace  f )
{
    return f->id();
};


//access id->edge
Solid::tEdge Solid::idEdge( int id0, int id1 )
{
    Edge e(id0,id1);
    return  m_edges.find( &e );//insert a new vertex, with id as the key
};

//access vertex->edge
Solid::tEdge Solid::vertexEdge( tVertex v0, tVertex v1 )
{
    int id0 = v0->id();
    int id1 = v1->id();
    Edge e(id0,id1);
    return m_edges.find( &e );//insert a new vertex, with id as the key
};


//access id->edge
//specified by source->target order(id0->id1)
Solid::tHalfEdge Solid::idHalfedge( int id0, int id1 )
{
    tEdge e = idEdge( id0, id1 );
    if ( ! e ) throw TopologyException();;
    tHalfEdge he = e->halfedge(0);
    if ( he->source()->id() == id0 && he->target()->id() == id1 ) return he;
    he = e->halfedge(1);
    if ( he->source()->id() != id0 || he->target()->id() != id1 ) throw TopologyException();
    return he;
};

//access vertex->edge
//specified by source->target order(v0->v1)
Solid::tHalfEdge Solid::vertexHalfedge( tVertex v0, tVertex v1 )
{
    tEdge e = vertexEdge( v0, v1 );
    tHalfEdge he = e->halfedge(0);
    if ( he->vertex() == v1 && he->he_prev()->vertex() == v0 ) return he;
    he = e->halfedge(1);
    if ( he->vertex() != v1 || he->he_prev()->vertex() != v0 ) throw TopologyException();
    return he;
};


//access e->v
Solid::tVertex Solid::edgeVertex1( Solid::tEdge  e )
{
    if ( e->halfedge(0 ) == NULL ) throw TopologyException();
    return e->halfedge(0)->source();
};

//access e->v
Solid::tVertex Solid::edgeVertex2( Solid::tEdge  e )
{
    if ( e->halfedge(0 ) == NULL ) throw TopologyException();
    return e->halfedge(0)->target();
};

//access e->f
Solid::tFace Solid::edgeFace1( Solid::tEdge  e )
{
    if ( e->halfedge(0) == NULL ) throw TopologyException();
    return e->halfedge(0)->face();
};

//access e->f
Solid::tFace Solid::edgeFace2( Solid::tEdge  e )
{
    if ( e->halfedge(1) == NULL ) throw TopologyException();
    return e->halfedge(1)->face();
};

//access he->f
Solid::tFace Solid::halfedgeFace( Solid::tHalfEdge  he )
{
    return he->face();
};


//access he->v
Solid::tVertex Solid::halfedgeVertex( Solid::tHalfEdge  he )
{
    return he->vertex();
};

bool Solid::isBoundary( Solid::tVertex  v )
{
    return v->boundary();
};

bool Solid::isBoundary( Solid::tEdge  e )
{
    if ( e->halfedge(0) == NULL || e->halfedge(1) == NULL ) return true;
    return false;
};

bool Solid::isBoundary( Solid::tHalfEdge  he )
{
    if ( he->he_sym() == NULL ) return true;
    return false;
};



//create new gemetric simplexes

Solid::tVertex Solid::createVertex( int id )
{
    tVertex v = new Vertex();

    v->id() = id;
    m_verts.insert( v );
    return v;//insert a new vertex, with id as the key
};


Solid::tEdge Solid::createEdge( int start, int end )
{
    Edge edge(start,end);

    tEdge e = m_edges.find( & edge );

    if ( e != NULL ) return e;

    e = new Edge( start, end );

    m_edges.insert( e );
    return e;

};

Solid::tFace Solid::createFace( int id )
{

    tFace f = new Face();

    f->id() = id;
    m_faces.insert( f );

    return f;//insert a new face, with id as the key
};

void Solid::copy( Solid & solid )
{   
    for ( AVL::TreeIterator<Vertex> viter(m_verts); !viter.end(); ++viter )
    {
        Vertex * v = *viter;
        Vertex * w = solid.createVertex( v->id() );

        w->point() = v->point();
        w->string()= v->string();
        w->boundary() = v->boundary();

		w->id2() = v->id2();
    }

    for ( AVL::TreeIterator<Face> fiter(m_faces); !fiter.end(); ++fiter )
    {
        Face * f = *fiter;

        HalfEdge * he = f->halfedge();
        int v[3];

        for ( int i = 0; i < 3; i ++ )
        {
            v[i] = he->vertex()->id();
            he = he->he_next();
        }

        Face * F = solid.createFace(v, f->id() );
        F->string() = f->string();
    }

    solid.labelBoundaryEdges();

};

void Solid::add_point_sphere(Point origin, double radius, double r, double g, double b)
{
	int bw = 2;
	int cur_verts = m_verts.getSize();
	int cur_faces = m_faces.getSize();

	Solid tmp;

	tmp.generate_regular_sphere(bw,cur_verts,cur_faces);
	tmp.close_std_shape(bw,cur_verts,cur_faces);

	for(SolidVertexIterator viter(&tmp); !viter.end(); ++viter)
	{

		Vertex * v = *viter;
		v->point() = v->point()*radius + origin;
	}

	if(r >= 0.0 && g >= 0.0 && b >= 0.0)
		tmp.color_solid(r,g,b);

	add(tmp);
}
	

void Solid::add_arrow(Point direction, Point origin, double scale, double r, double g, double b)
{
	int bw = 3;
	double radius = 0.1;
	int cur_verts = m_verts.getSize();
	int cur_faces = m_faces.getSize();

	double theta = direction.theta();
	double phi = direction.phi();

//	printf("theta, phi = %lf M_PI %lf M_PI\n", theta/M_PI, phi/M_PI);

	double length = direction.norm();

	Solid tmp;

//	printf("cur_verts, cur_faces = %i %i\n", cur_verts, cur_faces);
//	printf("about to generate_sphere\n");

	tmp.generate_regular_sphere(bw,cur_verts,cur_faces);
//	printf("I get thru generate_sphere\n");

	tmp.close_std_shape(bw,cur_verts,cur_faces);

//	printf("I get thru closed sphere\n");

//	tmp.write("tmp_sp.m");

	Point Np = tmp.idVertex(4*bw*bw + 2 + cur_verts)->point();
	Np = Point(0.0,0.0,1.0);
	tmp.idVertex(4*bw*bw + 2 + cur_verts)->point() = Np;

	Point Sp = tmp.idVertex(4*bw*bw + 1 + cur_verts)->point();
	Sp = Point(0.0,0.0,-1.0);
	tmp.idVertex(4*bw*bw + 1 + cur_verts)->point() = Sp;

	for(int i = 0; i < 2*bw - 1; i++)
		for(int j = 0; j < 2*bw; j++)
		{
			int ind = 2*bw*i + j + 1 + cur_verts;
			Point p = tmp.idVertex(ind)->point();
			double cur_rad = sqrt(p[0]*p[0] + p[1]*p[1]);
			cur_rad /= radius;

			p[0] /= cur_rad;
			p[1] /= cur_rad;

			if(i == 2*bw - 2)
				p[2] -= 0.2;

			if(i == 0)
				p[2] = -1.0;

			 tmp.idVertex(ind)->point() = p;

		}

	int ii = 2*bw - 1;
	for(int j = 0; j < 2*bw; j++)
	{
		int ind = 2*bw*ii + j + 1 + cur_verts;
		int ind2 = 2*bw*(ii-1) + j + 1 + cur_verts;
		Point p = tmp.idVertex(ind)->point();
		p[2] = tmp.idVertex(ind2)->point()[2];
		tmp.idVertex(ind)->point() = p;
	}

	for(SolidVertexIterator viter(&tmp); !viter.end(); ++viter)
	{
		Vertex *v = *viter;
		Point p = v->point();

		p[2] += 1.0;
		p[2] *= length;
		p /= (2.0/scale);
	//	p.print();
		p.EulerRotate(phi,-theta,0.0);
	//	printf("rotated\n");
	//	p.print();
		p += origin;
		v->point() = p;
	}

	if(r >= 0.0 && g >= 0.0 && b >= 0.0)
		tmp.color_solid(r,g,b);


//	printf("about to add\n");

	add(tmp);
//	printf("added\n"); fflush(stdout);
}


void Solid::MetricDistortion(Solid & solid, float * corners, float * edges, float * vert_angles, float * vert_dist)
{
	//float * corners = new float[mesh.numFaces() * 3];
	//float * edges = new float[mesh.numEdges()];
	//float * vert_angle = new float[mesh.numVertices()];
	//float * vert_dist = new float[mesh.numVertices()];
	//above allocation must be done prior to this funtion

	for (SolidVertexIterator viter(this); !viter.end(); ++viter)
	{
		Vertex * v = *viter;

		Point base = v->point();
		Point base_tar = solid.idVertex(v->id())->point();

		Point p0;
		Point p0_tar;
		Point p_prev;
		Point p_prev_tar;
		Point p_cur;
		Point p_cur_tar;

		int i = 0;
		double dist = 0.0;
		double angle_dif = 0.0;

		for (VertexVertexIterator vviter(v); !vviter.end(); ++vviter)
		{ 
			Vertex 	* vv = *vviter;
			if (i == 0)
			{
				p0 = vv->point() - base;
				p0_tar = solid.idVertex(vv->id())->point() - base_tar;
				p_prev = vv->point() - base;
				p_prev_tar = solid.idVertex(vv->id())->point() - base_tar;

				dist += fabs(log(p_prev.norm() / p_prev_tar.norm()));
			}
			else
			{
				p_cur = vv->point() - base;
				p_cur_tar = solid.idVertex(vv->id())->point() - base_tar;

				double angle1 = acos(p_prev*p_cur / (p_prev.norm()*p_cur.norm()));
				double angle2 = acos(p_prev_tar*p_cur_tar / (p_prev_tar.norm()*p_cur_tar.norm()));

				angle_dif += fabs(angle1 - angle2);

				p_prev = p_cur;
				p_prev_tar = p_cur_tar;
				
			}

			i++;

		}

		double angle1 = acos(p_prev*p0 / (p_prev.norm()*p0.norm()));
		double angle2 = acos(p_prev_tar*p0_tar / (p_prev_tar.norm()*p0_tar.norm()));	
		angle_dif += fabs(angle1 - angle2);

		vert_dist[v->id2() - 1] = dist / (double)i;
		vert_angles[v->id2() - 1] = angle_dif / (double)i;

	}

	int j = 0;

	for(SolidFaceIterator fiter(this); !fiter.end(); ++fiter)
	{
		Face * f = *fiter;
		for (FaceHalfedgeIterator fheiter(f); !fheiter.end(); ++ fheiter)
		{
			HalfEdge * he = *fheiter;

			Point p1 = he->source()->point() - he->target()->point();
			Point p2 = he->he_next()->target()->point() - he->target()->point();

			Point p1_tar = solid.idVertex(he->source()->id())->point() - solid.idVertex(he->target()->id())->point();
			Point p2_tar = solid.idVertex(he->he_next()->target()->id())->point() - solid.idVertex(he->target()->id())->point();

			double angle1 = acos(p1*p2 / (p1.norm()*p2.norm()));
			double angle2 = acos(p1_tar*p2_tar / (p1_tar.norm()*p2_tar.norm()));

			corners[j] = angle1 - angle2;

			j++;
		}  
	}

	j = 0;

	for (SolidEdgeIterator eiter(this); !eiter.end(); ++eiter)
	{
		Edge * e = *eiter;
		Edge * e_tar = solid.idEdge(e->ekey().s(), e->ekey().t());

		edges[j] = log(e->length() / e_tar->length());

		j++;
		

	}
}








void Solid::add( Solid & solid )
{   
//	int j = 0;
    for ( AVL::TreeIterator<Vertex> viter(solid.m_verts); !viter.end(); ++viter )
    {
        Vertex * v = *viter;
        Vertex * w = createVertex( v->id() );

        w->point() = v->point();
        w->string()= v->string();
        w->boundary() = v->boundary();
	//	j++;

	//	if (j == 1)
	//		printf("first vid %i\n", w->id());
    }

//	j = 0;
    for ( AVL::TreeIterator<Face> fiter(solid.m_faces); !fiter.end(); ++fiter )
    {
        Face * f = *fiter;

        HalfEdge * he = f->halfedge();
        int v[3];

        for ( int i = 0; i < 3; i ++ )
        {
            v[i] = he->vertex()->id();
            he = he->he_next();
        }

        Face * F = createFace(v, f->id() );
        F->string() = f->string();

//		j++;

//if (j == 1)
//			printf("first fid %i\n", f->id());
    }

    labelBoundaryEdges();

};

Vertex * Solid::edgeSplit( Edge * e )
{
    //find the max_id for vertices;
//	printf("in Solid::edgeSplit\n");

    AVL::TreeIterator<Vertex> viter(m_verts);

    int max_vid = -1;

    for ( ; !viter.end(); ++viter )
    {
        max_vid = ( max_vid > (*viter)->id() )? max_vid: (*viter)->id();
    }

//	printf("max vid found\n");

    AVL::TreeIterator<Face> fiter(m_faces);

    int max_fid = -1;

    for ( ; !fiter.end(); ++fiter )
    {
        max_fid = ( max_fid > (*fiter)->id() )? max_fid: (*fiter)->id();
    }

//	printf("max fid found\n");

    //create a new vertex

    tVertex nv = createVertex( ++max_vid );

//	printf("new vert created\n");


    tHalfEdge hes[2];
    hes[0] = e->halfedge(0);
    hes[1] = e->halfedge(1);

    nv->boundary() = (hes[1]==NULL);

    m_edges.remove( e );


    List<Face> new_faces;

    for ( int j = 0; j < 2; j ++ )
    {
        tVertex v[3];
        tHalfEdge he = hes[j];

        if ( he == NULL ) continue;

        tFace  f = he->face();
        m_faces.remove( f );
        delete f;

        for ( int i = 0; i < 3; i ++ )
        {
            v[i] = he->target();

            tEdge te = he->edge();
            if ( te->halfedge(0) == he )
            {
                te->halfedge(0) = te->halfedge(1);
            }
            te->halfedge(1) = NULL;

            he = he->he_next();
        }

        int k;
        for ( k = 0; k < 3; k ++ )
        {
            tHalfEdge ne = he->he_next();
            delete he;
            he = ne;
        }

        int vid[3];
        tVertex w[3];
        w[0] = nv; w[1] = v[0]; w[2] = v[1];
        for ( k = 0; k < 3; k ++ )
        {
            vid[k] = w[k]->id(); 
        }       

        Face * nf = createFace( vid, ++max_fid );
        new_faces.Append( nf );

        w[0] = nv; w[1] = v[1]; w[2] = v[2];
        for ( k = 0; k < 3; k ++ )
        {
            vid[k] = w[k]->id(); 
        }   
        nf = createFace( vid, ++ max_fid );
        new_faces.Append( nf );
    }   


    for ( ListIterator<Face> iter( new_faces ); !iter.end(); ++ iter )
    {
        Face * f = *iter;
        HalfEdge * he = f->halfedge();

        do
        {
            Edge     * e = he->edge();
            if ( e->halfedge(1) != NULL )
            {
                HalfEdge * h = e->halfedge(0);
                if ( h->target()->id() < h->source()->id() )
                {
                    e->halfedge(0) = e->halfedge(1);
                    e->halfedge(1) = h;
                }
            }
            he = he->he_next();
        }while ( he != f->halfedge() );
    }

    delete e;
    return nv;  
};


Vertex * Solid::trianglesubdivision( Face * f , Edge * e1 , Edge * e2 , Edge * e3 )
{
    //find the max_id for vertices;

    AVL::TreeIterator<Vertex> viter(m_verts);

    int max_vid = -1;

    for ( ; !viter.end(); ++viter )
    {
        max_vid = ( max_vid > (*viter)->id() )? max_vid: (*viter)->id();
    }

    AVL::TreeIterator<Face> fiter(m_faces);

    int max_fid = -1;

    for ( ; !fiter.end(); ++fiter )
    {
        max_fid = ( max_fid > (*fiter)->id() )? max_fid: (*fiter)->id();
    }

    //create a new vertex

    tVertex nv = createVertex( ++max_vid );

    List<Face> new_faces;
    List<Face> new_faces2;
    int k;
    int vid[3];
    tVertex fv[4];
    fv[0] = f->halfedge()->target();
    fv[1] = f->halfedge()->he_next()->target();
    fv[2] = f->halfedge()->he_next()->he_next()->target();
    fv[3] = f->halfedge()->he_next()->he_next()->he_next()->target();

    m_edges.remove( e1);
    m_edges.remove( e2);
    m_edges.remove( e3);
    m_faces.remove( f );


    tVertex w[3];
    w[0] = nv; w[1] = fv[0]; w[2] = fv[1];
    for ( k = 0; k < 3; k ++ )
    {
        vid[k] = w[k]->id(); 
    }       

    Face * nf = createFace( vid, ++max_fid );
    new_faces.Append( nf );

    w[0] = nv; w[1] = fv[1]; w[2] = fv[2];

    for ( k = 0; k < 3; k ++ )
    {
        vid[k] = w[k]->id(); 
    }       

    nf = createFace( vid, ++max_fid );
    new_faces.Append( nf);



    w[0] = nv; w[1] = fv[2]; w[2] = fv[3];
    for ( k = 0; k < 3; k ++ )
    {
        vid[k] = w[k]->id(); 
    }       

    nf = createFace( vid, ++max_fid );
    new_faces2.Append( nf );

	///commented from here on
/*
    w[0] = nv; w[1] = fv[2]; w[2] = fv[0];
    for( k = 0; k < 3; k ++ )
    {
        vid[k] = w[k]->id(); 
    }		

    nf = createFace( vid, ++max_fid );
    new_faces.Append( nf );

    m_faces.remove( f );
    


    for( ListIterator<Face> iter( new_faces ); !iter.end(); ++ iter )
    {
        Face * face = *iter;
        HalfEdge * he = face->halfedge();

        do{
            Edge     * e = he->edge();
            if( e->halfedge(1) != NULL )
            {
                HalfEdge * h = e->halfedge(0);
                if( h->target()->id() < h->source()->id() )
                {
                    e->halfedge(0) = e->halfedge(1);
                    e->halfedge(1) = h;
                }
            }
            he = he->he_next();
        }while( he != face->halfedge() );
    }

*/
    return nv;  
};

void Solid::reset_id2s()
{
	int count = 1;

	AVL::TreeIterator<Vertex> viter( m_verts );
    for ( ; !viter.end() ; ++ viter )
	{
		Vertex *v = *viter;
		v->id2() = count;
		count++;
	}

	/*count = 1;

	AVL::TreeIterator<Face> fiter( m_faces );
    for ( ; !fiter.end() ; ++ fiter )
	{
		Face *f = *fiter;
		f->id2() = count;
		count++;
	}*/
}



void Solid::add_point_face(Point & p, Face * f)
{
//	printf("in add_point_face\n");
	int i1 = f->halfedge()->target()->id();
    int i2 = f->halfedge()->he_next()->target()->id(); 
    int i3 = f->halfedge()->he_next()->he_next()->target()->id();
    tEdge e1 = idEdge(i1, i2);
    tEdge e2 = idEdge(i2, i3);
    tEdge e3 = idEdge(i3, i1);

	//double cds[2];
	//f->segment_cross(p - f->norm(),p+f->norm,&cds[0]);
	//Point np = f->apply_coords(&cds[0]);
	Point np = p;

    Vertex * nv = trianglesubdivision(f , e1 , e2 , e3);
    nv->point()=np;

}

void Solid::add_point_edge(Point & p, Edge *edge)
{
	//printf("in add_point_edge\n");
	Vertex * nv = edgeSplit( edge );
    nv->point() = p;
}


//returns a 4-vector with entries as follows:
//res(projection type, element id, coord1, coord2)		or
//res(projection type, element id1, element id1, coord)	
//projectiontype: 1.0 - vertex, 2.0 - edge, 3.0 - face
//element id: one entry for face, vertex, 2 entries for edge
//coord1, coord2 - for Face::apply_coords
//coord - for Edge::project_segment
Vector4 Solid::project_point(Point p, bool add_to_mesh)
{
//	printf("in project point\n");
	int i = 0;
	double min_dist = 0;
	double min_dist_face = 0;
	double min_dist_edge;
	int min_id = 0;
	int min_id_face = 0; 
	Edge *e_cur;
	double cds[2];

	double min_cds[2] = {0,0};
	double min_t;
	
	Edge * last_edge;
	double epsilon = 0.1;
	Vector4 res;
	//res(projection type, element id, coord1, coord2)

	Point proj;

	AVL::TreeIterator<Vertex> viter( m_verts );
    for ( ; !viter.end() ; ++ viter )
	{
		Vertex *v = *viter;
		double dist = (v->point() - p).norm();

		if(i == 0)
		{
			min_dist = dist;
			min_id = v->id();
			proj = v->point();
		}

		if(min_dist > dist)
		{
			min_dist = dist;
			min_id = v->id();
			proj = v->point();
		}
		i++;
	}
	Vertex * cur = idVertex(min_id);
	last_edge = cur->halfedge()->edge();

//	printf("closest vertex found, %i\n",min_id);
//	proj.print();

	i = 0;

	for(VertexFaceIterator vfiter(cur); !vfiter.end(); ++vfiter)
	{
		Face * f = *vfiter;
		Point n1 = p - f->norm()*min_dist;
		Point n2 = p + f->norm()*min_dist;
		if(f->segment_cross(n1,n2,&cds[0]))
		{
			Point loc_proj = f->apply_coords(&cds[0]);
			double	dist_face = (loc_proj - p).norm();

			if(i == 0)
			{
				proj = loc_proj;
				min_dist_face = dist_face;
				min_id_face = f->id();
				min_cds[0] = cds[0];
				min_cds[1] = cds[1];
			}
			if(min_dist_face > dist_face)
			{
				proj = loc_proj;
				min_dist_face = dist_face;
				min_id_face = f->id();
				min_cds[0] = cds[0];
				min_cds[1] = cds[1];
			}
			i++;
		}

	}

//	printf("closest face done, candidates = %i\n",i);
//	proj.print();

	if(i == 0)
	{
		for(VertexEdgeIterator veiter(idVertex(min_id)); !veiter.end(); ++veiter)
		{
			Edge * e = * veiter;
			Vertex * vv1 = idVertex(e->vertex(0));
			Vertex * vv2 = idVertex(e->vertex(1));
			Point p1 = vv1->point();
			Point p2 = vv2->point();

			Point dif = p-p1;
			Point curve_dir = p2 - p1;

			double t = curve_dir*dif/curve_dir.norm2();
			double line_dist = (dif - curve_dir*t).norm();
		
			if(t > 0 && t < 1) 
			{
				if(i == 0)
				{
					min_dist_edge = line_dist;
//					min_id_edge = e->other_vertex(cur)->id();//actually never used
					proj = p1 + curve_dir*t;
					e_cur = e;
					min_t = t;
				}

				if(min_dist_edge > line_dist)
				{
					min_dist_edge = line_dist;
				//	min_id_edge = e->other_vertex(cur)->id();
					proj = p1 + curve_dir*t;
					e_cur = e;
					min_t = t;
				}


				i++;
			}
		}
	//	printf("closest edge done, candidates = %i\n",i);
	//	proj.print();


		if(i == 0)//worst case, vertex is actually the projection to the space point
		{
		//	printf("in worst case\n");
			Point pp2 = last_edge->other_vertex(cur)->point();
		//	last_edge->
			
		//	pp2.print();
			proj += (pp2-proj)*epsilon;
		//	printf("last edge vertices: %i %i\n",last_edge->vertex(0), last_edge->vertex(1));
		//	exit(0);
			if(add_to_mesh)
				add_point_edge(proj,last_edge);
			res = Vector4(1.0,(double)min_id,0.0,0.0); //res(projection type, element id, coord1, coord2)
			return res;
		}
		else
		{
			if(add_to_mesh)
				add_point_edge(proj,e_cur);
			res = Vector4(2.0,(double)e_cur->ekey().s(),(double)e_cur->ekey().t(),min_t); //res(projection type, element id, coord1, coord2)
			return res;
		}
	//	printf("add point done\n");

	}
	else
	{
		//check if face projection is actually on edge
		Face *f_cur = idFace(min_id_face);


		if((proj - cur->point()).norm() < last_edge->length()*epsilon)//even worse, face normal aligned with displacement, but vertex is still the projection
		{
			Point pp2 = last_edge->other_vertex(cur)->point();
			proj += (pp2-proj)*epsilon;
			if(add_to_mesh)
				add_point_edge(proj,last_edge);
			res = Vector4(1.0,(double)min_id,0.0,0.0); //res(projection type, element id, coord1, coord2)
			return res;
		}
		else
		{
			Edge *e1 = f_cur->halfedge()->edge();
			Edge *e2 = f_cur->halfedge()->he_prev()->edge();
			Edge *e3 = f_cur->halfedge()->he_next()->edge();

			if(e1->include_point(proj)) 
			{
				Point pp2 = e1->other_vertex(cur)->point();
				proj += (pp2-proj)*epsilon;			
				if(add_to_mesh)
					add_point_edge(proj,e1);
				res = Vector4(1.0,(double)min_id,0.0,0.0); //res(projection type, element id, coord1, coord2)
				return res;
			}
			else if(e2->include_point(proj))
			{
				Point pp2 = e2->other_vertex(cur)->point();
				proj += (pp2-proj)*epsilon;
				if(add_to_mesh)
					add_point_edge(proj,e2);
				res = Vector4(1.0,(double)min_id,0.0,0.0); //res(projection type, element id, coord1, coord2)
				return res;
			}
			else if(e3->include_point(proj))
			{
				Point pp2 = e3->other_vertex(cur)->point();
				proj += (pp2-proj)*epsilon;
				if(add_to_mesh)
					add_point_edge(proj,e2);
				res = Vector4(1.0,(double)min_id,0.0,0.0); //res(projection type, element id, coord1, coord2)
				return res;
			}
			else //projection is actual internal face point
			{
				if(add_to_mesh)
					add_point_face(proj,f_cur);
				res = Vector4(3.0,(double)f_cur->id(),min_cds[0],min_cds[1]); //res(projection type, element id, coord1, coord2)
				return res;
			}

		}
//		printf("add point face done\n");

	}

	return res;
	
}

//uses a 4-vector with entries as follows:
//res(projection type, element id, coord1, coord2)		or
//res(projection type, element id1, element id1, coord)	
//projectiontype: 1.0 - vertex, 2.0 - edge, 3.0 - face
//element id: one entry for face, vertex, 2 entries for edge
//coord1, coord2 - for Face::apply_coords
//coord - for Edge::project_segment
Point Solid::apply_coords(Vector4 coords, int * min_id)
{
	Point proj;
	int type = round_to_int(coords[0]);
	int id1 = 0, id2 = 0;

//	printf("type = %i\n",type);

	if(type == 1)
	{
		id1 = round_to_int(coords[1]);
		if(id1 != 0)
		{
			if(idVertex(id1) == NULL)
			{
				printf("bad vertex %i\n", id1);
				exit(1);
			}

			proj = idVertex(id1)->point();
			if(min_id != NULL)
				min_id[0] = id1;
		}
		else
		{
			printf("in Solid:: apply_coords vertex id = 0, impossible\n");
			exit(1);
		}
	}
	else if(type == 2)
	{
		id1 = round_to_int(coords[1]);
		id2 = round_to_int(coords[2]);

		if(id1 != 0 && id2 != 0)
		{
			if(idEdge(id1,id2) == NULL)
			{
				printf("bad edge %i %i\n", id1, id2);
				exit(1);
			}

			proj = idEdge(id1,id2)->project_segment(coords[3]);
			if(min_id != NULL)
			{
				if(coords[3] < 0.5)
					min_id[0] = id1;
				else
					min_id[0] = id2;
			}
		}
		else
		{
			printf("in Solid:: apply_coords edge id = 0, impossible\n");
			exit(1);
		}
	}
	else if(type == 3)
	{
		id1 = round_to_int(coords[1]);
		double cds[2] = {coords[2],coords[3]};

		if(id1 != 0)
		{
			
			if(idFace(id1) == NULL)
			{
				printf("bad face %i\n", id1);
				exit(1);
			}

			proj = idFace(id1)->apply_coords(&cds[0]);		

			int i = 0;
			double min_dist;

			if(min_id != NULL)
			{
				Face *f = idFace(id1);

				for(FaceVertexIterator fviter(f); !fviter.end(); ++fviter)
				{
					Vertex *v = *fviter;
					if(i == 0)
					{
						min_dist = (v->point() - proj).norm();
						min_id[0] = v->id();
					}
					else if(min_dist > (v->point() - proj).norm())
					{
						min_dist = (v->point() - proj).norm();
						min_id[0] = v->id();
					}
				}
			}
		}
		else
		{
			printf("in Solid:: apply_coords face id = 0, impossible\n");
			exit(1);
		}
	}
	else
	{
		printf("in Solid:: apply_coords unknown projection type %i\n",type);
		exit(1);
	}


	return proj;
}


Vertex * Solid::closest_vert_to_point(Point p)
{
	int i = 0, min_id = 0;
	double min_dist = 1e5;
	AVL::TreeIterator<Vertex> viter( m_verts );
    for ( ; !viter.end() ; ++ viter )
	{
		Vertex *v = *viter;
		double dist = (v->point() - p).norm();

		if(i == 0)
		{
			min_dist = dist;
			min_id = v->id();			
		}

		if(min_dist > dist)
		{
			min_dist = dist;
			min_id = v->id();			
		}
		i++;
	}

	return idVertex(min_id);
}

void Solid::segment_regular_by_curve(int bw, Point * curve_orig, int N, float * mask, bool dilate_by_edge_length)
{
	/*double av_length = 0.0;
	for(SolidEdgeIterator eiter(this); !eiter.end(); ++eiter)
	{
		Edge *e = *eiter;
		av_length += e->length();
	}
	av_length /= (double)numEdges();*/

	double ratio = 3.0;
	mask[4*bw*bw] = 0.0;
	Solid sphere;
	sphere.generate_regular_sphere(bw);
	sphere.close_boundaries();
	close_boundaries();
	
	int * Ncrossings = new int[2*bw];
	Point * curve = new Point[N];	
	Point ** crossings = new Point*[2*bw];
	int ** cross_ids = new int*[2*bw];
	int ** curve_ids = new int*[2*bw];
	int even = 0, odd = 0;


	for(int i = 0; i < 2*bw; i++)
	{
		crossings[i] = new Point[N];
		cross_ids[i] = new int[N];
		curve_ids[i] = new int[N];
	}
	
	for(int n = 0; n < N; n++)
	{		
		Vector4 coords = project_point(curve_orig[n],false);
		curve[n] = sphere.apply_coords(coords);
	}

	//for(int i = 0; i < 4*bw*bw+2; i++)
	//	mask[i] = 0.0;
	

	
	Point Sval(0.0,0.0,-1.0);
	Point Nval(0.0,0.0,1.0);
	
	/*Point Sval;
	Point Nval;	
	for(int i = 0; i < 2*bw; i++)
		Sval += idVertex(i+1)->point();
	Sval = Sval/(double)(2*bw);

	for(int i = 4*bw*bw - 2*bw; i < 4*bw*bw; i++)
		Nval += idVertex(i+1)->point();
	Nval = Nval/(double)(2*bw);*/



	for(int i = 0; i < 2*bw; i++)
	{
		Ncrossings[i] = 0;
		for(int j = 0; j <= 2*bw; j++)
		{
			Point p1 = Sval;
			int ind1 = 4*bw*bw;
			if(j > 0)
			{
				ind1 = i + 2*bw*(j-1);
				p1 = sphere.idVertex(ind1 + 1)->point();
			}
			Point p2 = Nval - p1;
			int ind2 = 4*bw*bw + 1;
			if(j <= 2*bw - 1)
			{
				ind2 = i + 2*bw*(j);
				p2 = sphere.idVertex(ind2 + 1)->point() - p1;
			}

			int Ncrossings_loc = 0;

			for(int k = 0; k < N; k++)
			{
				Point w1 = curve[k];
				Point w2 = curve[(k+1)%N] - w1;

				double t,s;
				double dir_dot = w2*p2;

				//s - for mesh edge
				//t - for curve segment

				 // curve speeds orthogonal
				//	s = p2*(w1-p1)/p2.norm2();
				//	t = w2*(p1-w1)/w2.norm2();
				
				double t_num = p2*(p1-w1)*dir_dot + w2*(w1-p1)*p2.norm2();
				double t_denom = dir_dot*dir_dot - w2.norm2()*p2.norm2();
		

				if(fabs(t_denom)/(dir_dot*dir_dot) < 1e-4)//parallel segments
				{		
					t = -1.0;
					s = -1.0;
				}
				else
				{
					t = t_num/t_denom;
					s = (t*(dir_dot) + p2*(w1-p1))/p2.norm2();
				}
				

				if(s < 1.0 /*+ FLT_EPSILON*/ && s >= 0.0 - FLT_EPSILON 
					&& t < 1.0 /*+ FLT_EPSILON */&& t >= 0.0 - FLT_EPSILON)
				{
					Point c_pt = w1 + w2*t;
					Point m_pt = p1 + p2*s;

					if((c_pt - m_pt).norm() < ratio*p2.norm())//std::min(w2.norm(),p2.norm()))
					{
						

						bool too_close = false;
						if(Ncrossings[i] > 0)
							for(int ii = 0; ii < Ncrossings[i]; ii++)
							//	if((m_pt - crossings[i][ii]).norm() < 0.25*M_PI/(double)(4*bw))
								if((curve_ids[i][ii] == k || (curve_ids[i][ii] == (k+1)%N)) ||
									(m_pt - crossings[i][ii]).norm() < 0.25*M_PI/(double)(4*bw))//up to now 2-sphere was not assumed, now, it is.
								{
									too_close = true;
									break;
								}
						if(!too_close)
						{
							crossings[i][Ncrossings[i]] = m_pt;
							cross_ids[i][Ncrossings[i]] = j;
							curve_ids[i][Ncrossings[i]] = k;
							Ncrossings[i]++;
							Ncrossings_loc++;
						}
					}
				}
			}

			if(Ncrossings_loc%2 == 1)
				mask[ind2] = 1.0 - mask[ind1];
			else
				mask[ind2] = mask[ind1];	
		}

		/*if(Ncrossings[i]%2 == 1)
		{
			printf("Ncrossings = %i\n",Ncrossings);
			printf("i = %i\n",i);
		}*/
		if(Ncrossings[i]%2 == 1)
			odd++;
		else
			even++;
	}

	for(int i = 0; i < 2*bw; i++)
	{
		if((odd > even && Ncrossings[i]%2 == 0 && Ncrossings[i] > 0) || (odd < even && Ncrossings[i]%2 == 1 && Ncrossings[i] > 0))
		{
		//	printf("at id %i, crossings = %i\n",i,Ncrossings[i]);
			int erase_int_id = 0;
			double max_int_separation = 2.0;
			for(int j = 1; j < Ncrossings[i]; j++)
			{
				double sep = (crossings[i][j-1] - crossings[i][j]).norm();
		//		printf("id0, id1, sep: %i %i %lf\n",cross_ids[i][j-1],cross_ids[i][j],sep);
				if(sep < max_int_separation /*&& cross_ids[i][j] == cross_ids[i][j-1]*/)
				{
					max_int_separation = sep;
					erase_int_id = j;
				}
			}

			if(erase_int_id == Ncrossings[i] - 1)
				Ncrossings[i]--;
			else
			{
				for(int j = erase_int_id; j < Ncrossings[i] - 1; j++)
				{
					cross_ids[i][j] = cross_ids[i][j+1];
					crossings[i][j] = crossings[i][j+1];
				}
				Ncrossings[i]--;
			}

	//		printf("\nafter correction\n\n");
// 			for(int j = 1; j < Ncrossings[i]; j++)
// 			{
// 				double sep = (crossings[i][j-1] - crossings[i][j]).norm();
// 		//		printf("id0, id1, sep: %i %i %lf\n",cross_ids[i][j-1],cross_ids[i][j],sep);
// 			}

			for(int j = 0; j <= 2*bw; j++)
			{
		//		Point p1 = Sval;
				int ind1 = 4*bw*bw;
				if(j > 0)
				{
					ind1 = i + 2*bw*(j-1);
		//			p1 = sphere.idVertex(ind1 + 1)->point();
				}
		//		Point p2 = Nval - p1;
				int ind2 = 4*bw*bw + 1;
				if(j <= 2*bw - 1)
				{
					ind2 = i + 2*bw*(j);
		//			p2 = sphere.idVertex(ind2 + 1)->point() - p1;
				}

				int Ncrossings_loc = 0;

				for(int k = 0; k < Ncrossings[i]; k++)
				{
					//if(crossings[i][k][2] > p1[2] && crossings[i][k][2] <= (p1+p2)[2])
					if(cross_ids[i][k] == j)
						Ncrossings_loc++;
				}

				if(Ncrossings_loc%2 == 1)
					mask[ind2] = 1.0 - mask[ind1];
				else
					mask[ind2] = mask[ind1];	
			}
		}
	}

	reset_id2s();

	if(dilate_by_edge_length)
	{
		bool * border_vert = new bool[4*bw*bw + 2];
		float * tmp_mask = new float[4*bw*bw+2];
		double av_border_edge_length = 0.0;
		int Nedges = 0;
		double edge_ratio = 4.0;

		for(int i = 0; i < 4*bw*bw+2; i++)
			tmp_mask[i] = mask[i];

		for(SolidVertexIterator viter(this); !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			if(mask[v->id2()-1] == 1.0)
			for(VertexVertexIterator vviter(v); !vviter.end(); ++vviter)
			{
				Vertex *vv = *vviter;
				if(mask[vv->id2()-1] == 0.0)
					border_vert[vv->id2()-1] = true;
			}
		}

		for(SolidVertexIterator viter(this); !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			if(border_vert[v->id2()-1])
			for(VertexVertexIterator vviter(v); !vviter.end(); ++vviter)
			{
				Vertex *vv = *vviter;
				if(mask[vv->id2()-1] == 0.0)
				{
					av_border_edge_length += (v->point() - vv->point()).norm();
					Nedges++;
				}
			}
		}
		av_border_edge_length /= (double)Nedges;


		for(SolidVertexIterator viter(this); !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			if(border_vert[v->id2()-1])
			for(VertexVertexIterator vviter(v); !vviter.end(); ++vviter)
			{
				Vertex *vv = *vviter;
				if(av_border_edge_length*edge_ratio < (v->point() - vv->point()).norm())
				{
					tmp_mask[v->id2()-1] = 1.0;
				//	tmp_mask[v->id2()-1] = 1.0;
				}
			}
		}

		for(int i = 0; i < 4*bw*bw+2; i++)
			mask[i] = tmp_mask[i];

		delete [] tmp_mask;
		delete [] border_vert;
	}

	for(int i = 0; i < 2*bw; i++)
	{
		delete [] crossings[i];
		delete [] cross_ids[i];
		delete [] curve_ids[i];
	}
	delete[] curve_ids;
	delete[] cross_ids;
	delete[] crossings;
	delete[] curve;
	delete[] Ncrossings;
	
}

//void Solid::flow_onto_mesh(Solid * target, int grid_ratio, int iters, double step_size, double reg_weight)
//{
//	Point euler = target->align_main_axes();
//	Volume 

Vertex * Solid::closest_vert_to_point(Point p, int init_id, int s_rad)
{
	Vertex *guess = idVertex(init_id);
	bool _min = false;
	double min_dist = (p-guess->point()).norm();
	int rad = 0;

	while(rad < s_rad && !_min)
	{
		double cur_min = min_dist;
		int best_id = guess->id();
		for(VertexVertexIterator vviter(guess); !vviter.end(); ++vviter)
		{
			Vertex * v = *vviter;
			double dist = (p-v->point()).norm();

			if(dist < cur_min)
			{
				dist = cur_min;
				best_id = v->id();
			}
		}

		if(cur_min == min_dist)
			_min = true;
		else
			guess = idVertex(best_id);

		rad++;
	}

	return guess;
}

//assumes homologous meshes, initial guess given by correspondence
void Solid::Hausdorff(Solid * mesh, int s_rad, int fix_iters, int * ids)
{
		
	for(SolidVertexIterator viter(mesh); !viter.end(); ++viter)
	{
		Vertex *v = *viter;
		ids[v->id2() - 1] = closest_vert_to_point(v->point(),v->id(), s_rad)->id();
	}

	for(int j = 0; j < fix_iters; j++)
	{
		int num_fixes = 0;
		for(SolidVertexIterator viter(mesh); !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			Vertex *v_pr = idVertex(ids[v->id2() - 1]);
			Point p_pr = v_pr->point();
			Point p = v->point();


			for(VertexVertexIterator vviter(v); !vviter.end(); ++vviter)
			{
				Vertex *vv = *vviter;
				Vertex *vv_pr = idVertex(ids[vv->id2() - 1]);
				Point pp_pr = vv_pr->point();
				Point pp = vv->point();

				if( (pp-p).norm() +(p - p_pr).norm() < (pp_pr - pp).norm() )//wrong
				{
					ids[vv->id2() - 1] = closest_vert_to_point(v->point(),ids[v->id2() - 1], s_rad)->id();
					num_fixes++;
				}
			}
		}

		printf("at fix iter %i, %i fixes\n",j,num_fixes);
		if(num_fixes == 0)
			break;
	}
}



				

		


	



void Solid::write( const char * output)
{
	float *mask = new float[numVertices()];

	for(int i = 0; i < numVertices(); i++)
		mask[i] = 1.0;

	write( output, mask, 0.0 );

	delete [] mask;
}

		

void Solid::write( const char * output, float * mask, float thresh )
{
//	printf("resetting id2s in write\n");
	reset_id2s();
//	printf("done id2s in write\n");

    FILE * _os = fopen( output,"w");
    if ( ! _os ) throw FException();

    /* Set IO stream buffer */
    if ( setvbuf(_os,NULL,_IOFBF,65536) != 0 )
        printf("ERROR: Setting IO stream buffer failed on line: %d\n",__LINE__);

	fprintf(_os,"# CCBBM 1.0\n");


    //remove vertices
    AVL::TreeIterator<Vertex> viter( m_verts );
    for ( ; !viter.end() ; ++ viter )
    {
        tVertex v = *viter;

		if(mask[v->id2()-1] > thresh)
		{

        fprintf(_os,"Vertex %d ", v->id() );

        for ( int i = 0; i < 3; i ++ )
        {
            fprintf(_os,"%g ",v->point()[i]);
        }
        if ( v->string().size() > 0 )
        {
            fprintf(_os,"{%s}", v->string().c_str() );
        }
        fprintf(_os,"\n");
		}

    }

//	printf("done writing verts\n");

    for ( AVL::TreeIterator<Face> fiter(m_faces); !fiter.end(); ++ fiter )
    {
        tFace f = *fiter;
		bool to_print = true;

	//	printf("face %i\n",f->id());
		
		for(FaceVertexIterator fviter(f); !fviter.end(); ++fviter)
		{
			Vertex * v = *fviter;
		//	printf("id %i\n",v->id2());
			if(mask[v->id2()-1] <= thresh)
			{
				to_print = false;
				break;
			}
		}

	//	printf("done with verts\n");
		if(to_print)
		{


        fprintf(_os,"Face %d",f->id());

        tHalfEdge he = f->halfedge();
        do
        {
            fprintf(_os," %d", he->target()->id() );
            he = he->he_next();
        }while ( he != f->halfedge() );

        if ( f->string().size() > 0 )
        {
            fprintf(_os," {%s}", f->string().c_str() );
        }

        fprintf(_os,"\n");
		}
    }


//	printf("done writing faces\n");


    for ( AVL::TreeIterator<Edge> eiter(m_edges); !eiter.end(); ++ eiter )
    {
        tEdge e = *eiter;

		bool to_print = true;

		Vertex *v1, *v2;
		e->get_vertices(v1, v2);

		if(mask[v1->id2()-1] <= thresh)
			to_print = false;

		if(mask[v2->id2()-1] <= thresh)
			to_print = false;

        if ( e->string().size() > 0 && to_print)
        {
            fprintf(_os,"Edge %d %d {%s} \n", edgeVertex1(e)->id(), edgeVertex2(e)->id(), e->string().c_str() );
        }

    }

    for ( AVL::TreeIterator<Face> pfiter(m_faces); !pfiter.end(); ++ pfiter )
    {
        tFace f = *pfiter;

		bool to_print = true;


		for(FaceVertexIterator fviter(f); !fviter.end(); ++fviter)
		{
			Vertex * v = *fviter;
			if(mask[v->id2()-1] <= thresh)
			{
				to_print = false;
				break;
			}
		}

		if(to_print)
		{

        tHalfEdge he = f->halfedge();
        do
        {
            if ( he->string().size() > 0 )
            {
                fprintf(_os,"Corner %d %d {%s}\n", he->vertex()->id(), f->id(), he->string().c_str() );
            }
            he = he->he_next();
        }while ( he != f->halfedge() );
		}
    }

    fclose(_os);

//	printf("returning from write\n");

};



//method edge length
double Solid::edgeLength( Solid::tEdge e )
{
    Vertex * v1 = edgeVertex1(e);
    Vertex * v2 = edgeVertex2(e);

    return( v1->point() - v2->point() ).norm();
}


void
Solid:: labelBoundaryEdges()
{
    for ( AVL::TreeIterator<Edge> eiter( m_edges ); ! eiter.end() ; ++ eiter )
    {
        Edge  *   edge = *eiter;
        HalfEdge* he[2];

        he[0] = edge->halfedge(0);
        he[1] = edge->halfedge(1);

        //label boundary
        if ( he[1] == NULL )
        {
            he[0]->vertex()->boundary() = true;
            he[0]->he_prev()->vertex()->boundary()  = true;
        }

    }
}


void 
Solid::removeDanglingVertices()
{
    List<Vertex> dangling_verts;
    //Label boundary edges
    for ( AVL::TreeIterator<Vertex> viter( m_verts ); ! viter.end() ; ++ viter )
    {
        Solid::tVertex     v = *viter;
        if ( v->halfedge() != NULL ) continue;
        dangling_verts.Append( v );
    }

    for ( ListIterator<Vertex> iter( dangling_verts ); !iter.end(); ++ iter )
    {
        Solid::tVertex v = *iter;
        m_verts.remove( v );
        delete v;
        v = NULL;
    }
}

void Solid::readBYU( std::istream & is )
{
    char line[MAX_LINE];

    is.getline(line, MAX_LINE);
    std::string s(line );
    string_token_iterator iter(s, " \n");
    std::string str = *(++iter);
    int vnum = atoi( str.c_str() );
    str = *(++iter);
//	int fnum = atoi( str.c_str() );
    is.getline(line, MAX_LINE);
    int vind, find;
    vind = find = 1;

//	printf("starting to read BYU\n");
//	printf("vnum = %i\n",vnum);

    while ( is && !is.eof() && is.getline(line, MAX_LINE) )
    {
        if ( strlen( line ) == 0 ) continue;

        std::string s(line);

        string_token_iterator iter(s, " \n"); 
        std::string str;

        if ( vind <= vnum )
        {
            Point p;
            for ( int i=0; i<3; i++ )
            {
                str = *iter;
                p[i] = atof( str.c_str() );
                iter++;
            }

            tVertex v = createVertex( vind );
            v->point() = p;
            v->id() = vind;
            vind ++;
        }
        else
        {
            int v[3];
            for ( int i = 0; i < 3; i ++ )
            {
                str = *iter;
                v[i] = abs(atoi( str.c_str() ));
                iter++;
            }
           // v[2] = -v[2];

		//	printf("creating face %i %i %i\n",v[0],v[1],v[2]);
            createFace( v, find );
		//	printf("done face\n");

            find ++;
        }
    }

//	printf("geometry read\n");

    labelBoundaryEdges();

    removeDanglingVertices();

	reset_id2s();

//	printf("returning from read\n");

};

void Solid::writeBYU( std::ostream & os )
{
	char buffer_array[65536];
    os.rdbuf()->pubsetbuf(buffer_array, sizeof(buffer_array));
    // first two line
    os << "1 " << numVertices() << " " << numFaces() << " " << numVertices() + numFaces() - 1 << std::endl;
    os << "1 " << numFaces() << std::endl;

    AVL::TreeIterator<Vertex> viter( m_verts );
    for ( ; !viter.end() ; ++ viter )
    {
        tVertex v = *viter;

        for ( int i = 0; i < 3; i ++ )
        {
            os << "    " << std::fixed << v->point()[i];  
        }

        os << std::endl;
    }

    for ( AVL::TreeIterator<Face> fiter(m_faces); !fiter.end(); ++ fiter )
    {
        tFace f = *fiter;
        tHalfEdge he = f->halfedge();
        for ( int i=0; i<3; i ++ )
        {
            if ( i!= 2 )
            {
                os << "    ";
            }
            else
            {
                os << "    -"; 
            }
            os << std::right << he->target()->id2();
            he = he->he_next();
        }

        os << std::endl;
    }
};

//Based very strongly on code found in CT.
void Solid::UpdateNormals(void)
{
    _FaceNormal(this);
    _VertexNormal(this);
//	printf("thru normal update\n");
}


void Solid::_FaceNormal(MeshLib::Solid *pMesh)
{
    for ( MeshLib::SolidFaceIterator fiter( pMesh ); !fiter.end(); ++ fiter )
    {
        MeshLib::Face * f = *fiter;
        if ( f->trait() == NULL )
        {
            MeshLib::FaceNormalTrait * ft = new MeshLib::FaceNormalTrait;
            if ( ft == NULL ) return;
            MeshLib::add_trait<MeshLib::FaceNormalTrait,MeshLib::Face>(f,ft);
        }

        MeshLib::Point p[3];
        int i = 0;

        for ( MeshLib::FaceVertexIterator fviter(f); !fviter.end(); ++ fviter )
        {
            MeshLib::Vertex * v = *fviter;
            p[i++] = v->point();
        }
        MeshLib::Point n = (p[1]-p[0])^(p[2]-p[0]);
		if(n.norm() > 1e-6)
			n /= n.norm();
		else
			n = Point(1.0,0.0,0.0);

        f_normal(f) = n;
    }
}

void Solid::_VertexNormal(MeshLib::Solid *pMesh)
{
#if 0
    for ( MeshLib::SolidVertexIterator viter( pMesh ); !viter.end(); ++ viter )
    {
        MeshLib::Vertex * v = *viter;

        MeshLib::Point n(0,0,0);
        int count = 0;

        for ( MeshLib::VertexFaceIterator vfiter(v); !vfiter.end(); ++ vfiter )
        {
            MeshLib::Face * f = *vfiter;
            n += f_normal(f);
            count ++;
        }
        n /= count;

        if ( n.norm() > 1e-5 )
        {
            n /= n.norm();
        }

        v->normal() = n;
    }

    for ( MeshLib::SolidFaceIterator fiter( pMesh ); !fiter.end(); ++ fiter )
    {
        MeshLib::Face * f = *fiter;
        MeshLib::FaceNormalTrait * pfnt = pTrait<MeshLib::FaceNormalTrait,MeshLib::Face>( f );
        del_trait<MeshLib::FaceNormalTrait,MeshLib::Face>(f, pfnt);
    }
#else
    // using edge lengths as weight on calculating vertex normal
//	printf("I get here\n");
    for ( SolidFaceIterator fIter(this); !fIter.end(); ++fIter )
    {
        Face *pFace = *fIter;
        HalfEdge *pHE[3] = {pFace->halfedge(),
            pFace->halfedge()->he_next(), pFace->halfedge()->he_prev()};
        Point e[3];
        for ( int i = 0; i < 3; ++i )
        {
            e[i] = pHE[i]->target()->point() - pHE[i]->source()->point();
        }
        Point faceNormal = e[0] ^ e[1];
		if(faceNormal.norm() < 1e-6)
			faceNormal = Point(1.0,0.0,0.0);
        for ( int i = 0; i < 3; ++i )
        {
			double div = (e[(i + 1) % 3].norm2() * e[i].norm2());
			if(div > 1e-10)
				pHE[i]->target()->normal() += faceNormal / div;
        }
    }

//	printf("I get here too\n");


    for ( SolidVertexIterator vIter(this); !vIter.end(); ++vIter )
    {
        Vertex *pVer = *vIter;
		if(pVer->normal().norm() > 1e-10)
			pVer->normal() /= pVer->normal().norm();
    }
#endif
}

void Solid::color_solid(float * atr, double min, double max, bool reverse)
{
	for(int i = 0; i < numVertices(); i++)
		atr[i] -= min;
		
	max -= min;
	min = 0;
	
	int	i = 0;
	for(SolidVertexIterator viter(this);!viter.end();++viter)
	{
		Vertex *v_std = * viter;
		double r,g,b, disp = atr[i];
	
		//	printf("disp = %lf\n",disp); exit(0);
		if(disp > max)
		{
			g = 0.0;
			b = 1.0;
			r = 0.0;
		}

		else if(disp < min)
		{
			g = 0.0;
			b = 0.0;
			r = 1.0;
		}

		else
		{	
			if(disp > ((max-min)/2) )
			{				
				g = 1.0 - ((disp - (max/2) )/(max/2) );
				b = 1.0 - g;
				r = 0.0;
			}
			else
			{
				g = disp/(max/2);
				r = 1.0 - g;
				b = 0.0;
			}
		}
		
		char line[1024];
		if(!reverse)
			sprintf(line,"rgb=(%lf %lf %lf)",r,g,b);
		else
			sprintf(line,"rgb=(%lf %lf %lf)",b,g,r);

		v_std->string() = std::string(line);
		i++;
	}
}


void Solid::color_solid(float * atr, bool reverse)
{
	double min = atr[0];
	double max = atr[0];

	for(int i = 0; i < numVertices(); i++)
	{
		if(min > atr[i]) min = atr[i];
		if(max < atr[i]) max = atr[i];
	}

	color_solid(atr, min, max, reverse);
}

void Solid::set_scalar_attribute(float * atr, char* key)
{
	int i = 0;
	for(SolidVertexIterator viter(this); !viter.end(); ++viter)
	{
		Vertex * v = *viter;
		v->set_scalar_attribute(atr[i],key);
		i++;
	}
}

void Solid::get_scalar_attribute(float * atr, char* key, double default_val)
{
	int i = 0;
	for (SolidVertexIterator viter(this); !viter.end(); ++viter)
	{
		Vertex * v = *viter;
		if (!v->get_scalar_attribute(&atr[i], key))
			atr[i] = default_val;
		i++;
	}
}

void Solid::get_vector_attribute(Point * atr, char* key, Point default_val)
{
	int i = 0;
	for(SolidVertexIterator viter(this); !viter.end(); ++viter)
	{
		Vertex * v = *viter;	
		if(!v->get_vector_attribute(&atr[i],key))
			atr[i] = default_val;
		i++;
	}	
}


}//end of namespace MeshLib










