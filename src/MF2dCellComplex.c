#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <MFAtlas.h>
#include <math.h>

struct MFVertexSt;
typedef struct MFVertexSt *MFVertex;

struct MFEdgeSt;
typedef struct MFEdgeSt *MFEdge;

struct MFTriangleSt;
typedef struct MFTriangleSt *MFTriangle;

struct MFVertexSt
 {
  int iPoly;

  int nC;
  double *center;

  int nNeighbors;
  int mNeighbors;
  int *neighbor;

  int nCobnd;
  int mCobnd;
  MFEdge *e;
  int *pm;

  int inverse;

  MFVertex next;
  MFVertex prev;
 };

struct MFEdgeSt
 {
  int iEdge;

  int io;
  int id;

  int pmo;
  MFVertex o;
  int pmd;
  MFVertex d;

  int l0;
  int l1;
  int l2;

  int r0;
  int r1;
  int r2;

  int pml;
  MFTriangle l;
  int pmr;
  MFTriangle r;

  int counter;
  int inverse;

  MFEdge next;
  MFEdge prev;
 };

struct MFTriangleSt
 {
  int iTri;

  int i0;
  int i1;
  int i2;

  int    pm0;
  MFEdge e0;
  int    pm1;
  MFEdge e1;
  int    pm2;
  MFEdge e2;

  int inverse;

  MFTriangle next;
  MFTriangle prev;
 };

#ifdef MFDUALLINKEDLIST

struct MF2dCellComplexSt
 {
  MFVertex firstVertex;
  MFVertex lastVertex;
  MFVertex *vertexIndex;

  MFEdge firstEdge;
  MFEdge lastEdge;
  MFEdge *edgeIndex;

  MFTriangle firstTriangle;
  MFTriangle lastTriangle;
  MFTriangle *triangleIndex;
 };
typedef struct MF2dCellComplexSt *MF2dCellComplex;

#else

struct MF2dCellComplexSt
 {
  MFVertex *vertexList;
  int nVertices;
  int mVertices;

  MFEdge *edgeList;
  int nEdges;
  int mEdges;

  MFTriangle *triangleList;
  int nTriangles;
  int mTriangles;
 };
typedef struct MF2dCellComplexSt *MF2dCellComplex;

#endif

/*
 *   Put lists into a "complex" data structure
 *   Build nPoly long index into the (sorted) edge and triangle lists. This would speed the lookups.
 *      (each index is associated with at most nVert others).
 *      move to linked lists and insert to keep sorted.
 *      Break out the building of the data structure from the export to file.
 *      Once dual is constructed "fix up" the Voronoi to match.
 */

static MFVertex MFAddVertex(MF2dCellComplex,int iPoly,int nC,double *center);
static MFVertex MFFindVertex(MF2dCellComplex,int iPoly);
static void VertexAddNeighbor(MFVertex V,int neighbor);
static void VertexAddCoboundary(MFVertex V,int pm,MFEdge E);
static MFEdge MFAddEdge(MF2dCellComplex,int o,int d, int l0, int l1, int l2, int r0, int r1, int r2);
static void MFRemoveEdge(MF2dCellComplex,MFVertex,MFVertex);
static void MFRemoveEdgeByIndex(MF2dCellComplex,int);
static MFEdge MFFindEdge(MF2dCellComplex,MFVertex o,MFVertex d);
static void MFSetEdgeDuals(MFEdge e, MFTriangle l, MFTriangle r);
static MFTriangle MFAddTriangle(MF2dCellComplex,int i0, int i1, int i2);
static MFTriangle MFAddTriangleAndSetEdges(MF2dCellComplex,MFVertex a, MFEdge ab, MFVertex b, MFEdge bc, MFVertex c, MFEdge ac);
static MFTriangle MFFindTriangle(MF2dCellComplex,MFVertex a,MFVertex b, MFVertex c);
static void MFRemoveTriangle(MF2dCellComplex,MFVertex,MFVertex,MFVertex);
static int VertexOrder(const void *l,const void *r);
static int EdgeOrder(const void *l,const void *r);
static int TriangleOrder(const void *l,const void *r);
static MFVertex MFPivotTriangle(MF2dCellComplex,MFVertex v,MFVertex b0, MFVertex b1);
static void SortLists(MF2dCellComplex);
static void PrintItAll(MF2dCellComplex);
static void PrintVertex(int i, MFVertex v);
static void PrintEdge(int i, MFEdge e);
static void PrintTriangle(int i, MFTriangle t);

static void MFFreeMFVertex(MFVertex v, MFErrorHandler e)
 {
  free(v->center);
  free(v->e);
  free(v);
  return;
 }

static void MFFreeMFEdge(MFEdge edge, MFErrorHandler e)
 {
  free(edge);
  return;
 }

static void MFFreeMFTriangle(MFTriangle t, MFErrorHandler e)
 {
  free(t);
  return;
 }

static MFVertex MFAddVertex(MF2dCellComplex C, int iPoly,int nC,double *center)
 {
  int i;

  if(C->nVertices<=C->mVertices)
   {
    C->mVertices+=100;
    C->vertexList=(MFVertex*)realloc(C->vertexList,C->mVertices*sizeof(MFVertex));
    if(C->vertexList==NULL)
     {
      fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
      exit(12);
     }
   }
  C->vertexList[C->nVertices]=(MFVertex)malloc(sizeof(struct MFVertexSt));
  if(C->vertexList[C->nVertices]==NULL)
   {
    fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
    exit(12);
   }

  (C->vertexList[C->nVertices])->iPoly=iPoly;
  (C->vertexList[C->nVertices])->nC=nC;
  (C->vertexList[C->nVertices])->center=(double*)malloc(nC*sizeof(double));
  if((C->vertexList[C->nVertices])->center==NULL)
   {
    fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
    exit(12);
   }
  for(i=0;i<nC;i++)(C->vertexList[C->nVertices])->center[i]=center[i];
  (C->vertexList[C->nVertices])->nCobnd=0;
  (C->vertexList[C->nVertices])->mCobnd=5;
  (C->vertexList[C->nVertices])->e =(MFEdge*)malloc(5*sizeof(MFEdge));
  if((C->vertexList[C->nVertices])->e==NULL)
   {
    fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
    exit(12);
   }
  (C->vertexList[C->nVertices])->pm=(int*)malloc(5*sizeof(int));
  if((C->vertexList[C->nVertices])->pm==NULL)
   {
    fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
    exit(12);
   }

  (C->vertexList[C->nVertices])->nNeighbors=0;
  (C->vertexList[C->nVertices])->mNeighbors=10;
  (C->vertexList[C->nVertices])->neighbor =(int*)malloc(10*sizeof(int));
  if((C->vertexList[C->nVertices])->neighbor==NULL)
   {
    fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
    exit(12);
   }

  C->vertexList[C->nVertices]->inverse=C->nVertices;

  C->nVertices++;
  return C->vertexList[C->nVertices-1];
 }

static MFVertex MFFindVertex(MF2dCellComplex C, int iPoly)
 {
  struct MFVertexSt TMPVERTEX;
  MFVertex *tmpVertexPtr;
  MFVertex v;

  v=&TMPVERTEX;
  v->iPoly=iPoly;

  tmpVertexPtr=(MFVertex*)bsearch(&v,C->vertexList,C->nVertices,sizeof(MFVertex),VertexOrder);

  if(tmpVertexPtr!=NULL)return *tmpVertexPtr;
   else return NULL;
 }

static void VertexAddNeighbor(MFVertex V,int neighbor)
 {
  if(V->nNeighbors<=V->mNeighbors)
   {
    V->mNeighbors+=5;
    V->neighbor=(int*)realloc(V->neighbor,V->mNeighbors*sizeof(int));
    if(V->neighbor==NULL)
     {
      fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
      exit(12);
     }
    V->pm=(int*)realloc(V->pm,V->mNeighbors*sizeof(int));
    if(V->pm==NULL)
     {
      fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
      exit(12);
     }
   }
  V->neighbor[V->nNeighbors]=neighbor;

  (V->nNeighbors)++;
  
  return;
 }

static void VertexAddCoboundary(MFVertex V,int pm,MFEdge E)
 {
  if(V->nCobnd<=V->mCobnd)
   {
    V->mCobnd+=5;
    V->e=(MFEdge*)realloc(V->e,V->mCobnd*sizeof(MFEdge));
    if(V->e==NULL)
     {
      fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
      exit(12);
     }

    V->pm=(int*)realloc(V->pm,V->mCobnd*sizeof(int));
    if(V->pm==NULL)
     {
      fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
      exit(12);
     }
   }
  V->e[V->nCobnd]=E;
  V->pm[V->nCobnd]=pm;

  (V->nCobnd)++;
  
  return;
 }

static MFEdge MFAddEdge(MF2dCellComplex C, int o,int d, int l0, int l1, int l2, int r0, int r1, int r2)
 {
  if(o==-1)return NULL;
  if(d==-1)return NULL;

  if(C->nEdges<=C->mEdges)
   {
    C->mEdges+=100;
    C->edgeList=(MFEdge*)realloc(C->edgeList,C->mEdges*sizeof(MFEdge));
    if(C->edgeList==NULL)
     {
      fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
      exit(12);
     }
   }
  C->edgeList[C->nEdges]=(MFEdge)malloc(sizeof(struct MFEdgeSt));
  if(C->edgeList[C->nEdges]==NULL)
   {
    fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
    exit(12);
   }
  (C->edgeList[C->nEdges])->iEdge=C->nEdges;

  (C->edgeList[C->nEdges])->io=o;
  (C->edgeList[C->nEdges])->id=d;

  (C->edgeList[C->nEdges])->pmo=1;
  (C->edgeList[C->nEdges])->o=NULL;
  (C->edgeList[C->nEdges])->pmd=-1;
  (C->edgeList[C->nEdges])->d=NULL;

  if(l0<r0)
   {
    (C->edgeList[C->nEdges])->l0=l0;
    (C->edgeList[C->nEdges])->l1=l1;
    (C->edgeList[C->nEdges])->l2=l2;

    (C->edgeList[C->nEdges])->r0=r0;
    (C->edgeList[C->nEdges])->r1=r1;
    (C->edgeList[C->nEdges])->r2=r2;
   }else{
    (C->edgeList[C->nEdges])->l0=r0;
    (C->edgeList[C->nEdges])->l1=r1;
    (C->edgeList[C->nEdges])->l2=r2;

    (C->edgeList[C->nEdges])->r0=l0;
    (C->edgeList[C->nEdges])->r1=l1;
    (C->edgeList[C->nEdges])->r2=l2;
   }

  (C->edgeList[C->nEdges])->pml=1;
  (C->edgeList[C->nEdges])->l=NULL;
  (C->edgeList[C->nEdges])->pmr=-1;
  (C->edgeList[C->nEdges])->r=NULL;

  C->edgeList[C->nEdges]->inverse=C->nEdges;
  C->edgeList[C->nEdges]->counter=0;

  C->nEdges++;
  return C->edgeList[C->nEdges-1];
 }

static MFEdge MFFindEdge(MF2dCellComplex C, MFVertex o,MFVertex d)
 {
  struct MFEdgeSt TMPEDGE;
  MFEdge *tmpEdgePtr;
  MFEdge e;

  if(o==NULL||d==NULL)return NULL;

/* Assumes that the triangle list is sorted */

  e=&TMPEDGE;

  if(o->iPoly<d->iPoly)
   {
    e->io=o->iPoly;
    e->id=d->iPoly;
   }else{
    e->io=d->iPoly;
    e->id=o->iPoly;
   }
  tmpEdgePtr=(MFEdge*)bsearch(&e,C->edgeList,C->nEdges,sizeof(MFEdge),EdgeOrder);

  if(tmpEdgePtr!=NULL)return *tmpEdgePtr;
   else return NULL;
 }

static void MFRemoveEdgeByIndex(MF2dCellComplex C, int i)
 {
  int j,l,f;
  MFEdge EP;

  EP=C->edgeList[i];


  for(j=i+1;j<C->nEdges;j++)
   {
    (C->edgeList[j]->inverse)--;
    C->edgeList[j-1]=C->edgeList[j];
   }
  C->nEdges--;

  for(j=0;j<C->nVertices;j++)
   {
    for(l=0;l<C->vertexList[j]->nCobnd;l++)
     {
      if(((C->vertexList[j])->e[l])->iEdge==EP->iEdge)
       {
        for(f=l+1;f<C->vertexList[j]->nCobnd;f++)
         {
          ((C->vertexList[j])->e)[f-1]=((C->vertexList[j])->e)[f];
          ((C->vertexList[j])->pm)[f-1]=((C->vertexList[j])->pm)[f];
         }
        (C->vertexList[j]->nCobnd)--;
        l--;
       }
     }
   }

/*free(EP);*/
  EP->iEdge=-1;
  return;
 }

static void MFRemoveEdge(MF2dCellComplex C, MFVertex a,MFVertex b)
 {
  int i,j,l,f;
  MFEdge EP;

  EP=MFFindEdge(C,a,b);
  if(EP==NULL)
   {
    printf("Problem in RemoveEdge: edge (%d,%d) not found\n",a->iPoly,b->iPoly);
    return;
   }


  i=EP->inverse;
  if(0){printf("Remove Edge %5d E[%5d]=(%5d,%5d)\n",i,EP->iEdge,a->iPoly,b->iPoly);fflush(stdout);}
  for(j=i+1;j<C->nEdges;j++)
   {
    (C->edgeList[j]->inverse)--;
    C->edgeList[j-1]=C->edgeList[j];
   }
  C->nEdges--;

  for(j=0;j<C->nVertices;j++)
   {
    for(l=0;l<C->vertexList[j]->nCobnd;l++)
     {
      if(((C->vertexList[j])->e[l])->iEdge==EP->iEdge)
       {
        for(f=l+1;f<C->vertexList[j]->nCobnd;f++)
         {
          ((C->vertexList[j])->e)[f-1]=((C->vertexList[j])->e)[f];
          ((C->vertexList[j])->pm)[f-1]=((C->vertexList[j])->pm)[f];
         }
        (C->vertexList[j]->nCobnd)--;
        l--;
       }
     }
   }

/*free(EP);*/
  EP->iEdge=-1;
  return;
 }
 
static void MFSetEdgeDuals(MFEdge e, MFTriangle l, MFTriangle r)
 {
  if(e==NULL)return;

  e->l=l;
  e->r=r;

  return;
 }

static MFTriangle MFAddTriangle(MF2dCellComplex C, int i0, int i1, int i2)
 {
  int i;
  MFTriangle T;
  MFEdge E;

  if(i1<i0){ i=i0; i0=i1;i1=i;}
  if(i2<i0){ i=i0; i0=i2;i2=i;}
  if(i2<i1){ i=i1; i1=i2;i2=i;}

  if(C->nTriangles<=C->mTriangles)
   {
    C->mTriangles+=100;
    C->triangleList=(MFTriangle*)realloc(C->triangleList,C->mTriangles*sizeof(MFTriangle));
    if(C->triangleList==NULL)
     {
      fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
      exit(12);
     }
   }
  C->triangleList[C->nTriangles]=(MFTriangle)malloc(sizeof(struct MFTriangleSt));
  if(C->triangleList[C->nTriangles]==NULL)
   {
    fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
    exit(12);
   }
  T=C->triangleList[C->nTriangles];
  T->iTri=C->nTriangles;
  T->inverse=C->nEdges;

  C->nTriangles++;

  T->i0=i0;
  T->i1=i1;
  T->i2=i2;

  T->e0=NULL;
  T->e1=NULL;
  T->e2=NULL;

  T->pm0=1;
  T->pm1=1;
  T->pm2=1;

  return T;
 }

static MFTriangle MFAddTriangleAndSetEdges(MF2dCellComplex C, MFVertex a, MFEdge ab, MFVertex b, MFEdge bc, MFVertex c, MFEdge ac)
 {
  int i;
  int i0,i1,i2;
  MFVertex V;
  MFEdge E;
  MFTriangle T;

  if(a==NULL)return;
  if(b==NULL)return;
  if(c==NULL)return;
  if(ab==NULL)return;
  if(bc==NULL)return;
  if(ac==NULL)return;

  i0=a->iPoly;
  i1=b->iPoly;
  i2=c->iPoly;

  if(i1<i0){ i=i0; i0=i1;i1=i;V=a;a=b;b=V;E=bc;bc=ac;ac=E;}
  if(i2<i0){ i=i0; i0=i2;i2=i;V=a;a=c;c=V;E=ab;bc=ab;bc=E;}
  if(i1<i0){ i=i0; i0=i1;i1=i;V=a;a=b;b=V;E=bc;bc=ac;ac=E;}

  if(C->nTriangles<=C->mTriangles)
   {
    C->mTriangles+=100;
    C->triangleList=(MFTriangle*)realloc(C->triangleList,C->mTriangles*sizeof(MFTriangle));
    if(C->triangleList==NULL)
     {
      fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
      exit(12);
     }
   }
  C->triangleList[C->nTriangles]=(MFTriangle)malloc(sizeof(struct MFTriangleSt));
  if(C->triangleList[C->nTriangles]==NULL)
   {
    fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
    exit(12);
   }
  T=C->triangleList[C->nTriangles];
  T->iTri=C->nTriangles;
  T->inverse=C->nTriangles;

  C->nTriangles++;

  T->i0=i0;
  T->i1=i1;
  T->i2=i2;

  T->e0=ab;
  T->e1=bc;
  T->e2=ac;

  if(ab->o==a)T->pm0= 1;
   else       T->pm0=-1;
  if(bc->o==b)T->pm1= 1;
   else       T->pm1=-1;
  if(ac->o==a)T->pm2= 1;
   else       T->pm2=-1;

  if(0){printf("Add    Triangle %5d T[%5d] = (%d,%d,%d)\n",C->nTriangles-1,T->iTri,i0,i1,i2);fflush(stdout);}

  (ab->counter)++;
  (bc->counter)++;
  (ac->counter)++;

  return T;
 }

static MFTriangle MFFindTriangle(MF2dCellComplex C, MFVertex a,MFVertex b, MFVertex c)
 {
  struct MFTriangleSt TMPTRIANGLE;
  MFTriangle *tmpTrianglePtr;
  MFTriangle T;
  int l0,l1,l2,l;

/* Assumes that the triangle list is sorted */

  l0=a->iPoly;
  l1=b->iPoly;
  l2=c->iPoly;
  if(l1<l0){ l=l0; l0=l1;l1=l;}
  if(l2<l0){ l=l0; l0=l2;l2=l;}
  if(l2<l1){ l=l1; l1=l2;l2=l;}

  T=&TMPTRIANGLE;

  T->i0=l0;
  T->i1=l1;
  T->i2=l2;
  tmpTrianglePtr=(MFTriangle*)bsearch(&T,C->triangleList,C->nTriangles,sizeof(MFTriangle),TriangleOrder);

  if(tmpTrianglePtr!=NULL)return *tmpTrianglePtr;
   else return NULL;
 }

static void MFRemoveTriangle(MF2dCellComplex C, MFVertex a,MFVertex b,MFVertex c)
 {
  int i,j,l;
  MFTriangle TP;

  TP=MFFindTriangle(C,a,b,c);
  if(TP==NULL)
   {
    printf("Problem in RemoveTriangle: triangle (%d.%d,%d) not found\n",a->iPoly, b->iPoly, c->iPoly);
    return;
   }

  i=TP->inverse;
  if(0){printf("Remove Triangle %5d T[%5d] = (%d,%d,%d)\n",i,TP->iTri,a->iPoly,b->iPoly,c->iPoly);fflush(stdout);}
  for(j=i+1;j<C->nTriangles;j++)
   {
    (C->triangleList[j]->inverse)--;
    C->triangleList[j-1]=C->triangleList[j];
   }
  C->nTriangles--;

/* Decrement the edge counters. */

 if(TP->e0!=NULL&&(TP->e0)->counter>0)((TP->e0)->counter)--;
 if(TP->e1!=NULL&&(TP->e1)->counter>0)((TP->e1)->counter)--;
 if(TP->e2!=NULL&&(TP->e2)->counter>0)((TP->e2)->counter)--;

/* Must update the duals of the edges if they point to T[i]. */

  for(j=0;j<C->nEdges;j++)
   {
    if((C->edgeList[j])->l!=NULL && ((C->edgeList[j])->l)->iTri==TP->iTri || (C->edgeList[j])->r!=NULL && ((C->edgeList[j])->r)->iTri==TP->iTri)
     {
      for(l=0;l<C->nTriangles;l++)
       {
        if(C->triangleList[l]->e0!=NULL && (C->triangleList[l]->e0)->iEdge==(C->edgeList[j])->iEdge)
         {
          if((C->edgeList[j])->l!=NULL && ((C->edgeList[j])->l)->iTri==TP->iTri)
           {
            (C->edgeList[j])->l=C->triangleList[l];
            (C->edgeList[j])->l0=C->triangleList[l]->i0;
            (C->edgeList[j])->l1=C->triangleList[l]->i1;
            (C->edgeList[j])->l2=C->triangleList[l]->i2;
           }else if((C->edgeList[j])->r!=NULL && ((C->edgeList[j])->r)->iTri==TP->iTri)
           {
            (C->edgeList[j])->r=C->triangleList[l];
            (C->edgeList[j])->r0=C->triangleList[l]->i0;
            (C->edgeList[j])->r1=C->triangleList[l]->i1;
            (C->edgeList[j])->r2=C->triangleList[l]->i2;
           }
         }
        if(C->triangleList[l]->e1!=NULL && (C->triangleList[l]->e1)->iEdge==(C->edgeList[j])->iEdge)
         {
          if((C->edgeList[j])->l!=NULL && ((C->edgeList[j])->l)->iTri==TP->iTri)
           {
            (C->edgeList[j])->l=C->triangleList[l];
            (C->edgeList[j])->l0=C->triangleList[l]->i0;
            (C->edgeList[j])->l1=C->triangleList[l]->i1;
            (C->edgeList[j])->l2=C->triangleList[l]->i2;
           }else if((C->edgeList[j])->r!=NULL && ((C->edgeList[j])->r)->iTri==TP->iTri)
           {
            (C->edgeList[j])->r=C->triangleList[l];
            (C->edgeList[j])->r0=C->triangleList[l]->i0;
            (C->edgeList[j])->r1=C->triangleList[l]->i2;
            (C->edgeList[j])->r2=C->triangleList[l]->i2;
           }
         }
        if(C->triangleList[l]->e2!=NULL && (C->triangleList[l]->e2)->iEdge==(C->edgeList[j])->iEdge)
         {
          if((C->edgeList[j])->l!=NULL && ((C->edgeList[j])->l)->iTri==TP->iTri)
           {
            (C->edgeList[j])->l=C->triangleList[l];
            (C->edgeList[j])->l0=C->triangleList[l]->i0;
            (C->edgeList[j])->l1=C->triangleList[l]->i1;
            (C->edgeList[j])->l2=C->triangleList[l]->i2;
           }else if((C->edgeList[j])->r!=NULL && ((C->edgeList[j])->r)->iTri==TP->iTri)
           {
            (C->edgeList[j])->r=C->triangleList[l];
            (C->edgeList[j])->r0=C->triangleList[l]->i0;
            (C->edgeList[j])->r1=C->triangleList[l]->i1;
            (C->edgeList[j])->r2=C->triangleList[l]->i2;
           }
         }
       }
     }
   }

/*free(TP);*/
  TP->iTri=-1;
  return;
 }

static int VertexOrder(const void *l,const void *r)
 {
  MFVertex vl;
  MFVertex vr;

  int less=0;

  vl=((MFVertex*)l)[0];
  vr=((MFVertex*)r)[0];

  if(vl->iPoly<vr->iPoly)less=-1;
   else if(vl->iPoly>vr->iPoly)less=1;
   else less=0;

  return less;
 }

static int EdgeOrder(const void *l,const void *r)
 {
  MFEdge el;
  MFEdge er;

  int less=0;

  el=((MFEdge*)l)[0];
  er=((MFEdge*)r)[0];

  if(el->io<er->io)less=-1;
   else if(el->io>er->io)less=1;
   else if(el->id<er->id)less=-1;
   else if(el->id>er->id)less=1;
   else less=0;

  return less;
 }

static int TriangleOrder(const void *l,const void *r)
 {
  MFTriangle tl;
  MFTriangle tr;
  int less=0;

  tl=((MFTriangle*)l)[0];
  tr=((MFTriangle*)r)[0];

  if(tl->i0<tr->i0)less=-1;
   else if(tl->i0>tr->i0)less=1;
   else if(tl->i1<tr->i1)less=-1;
   else if(tl->i1>tr->i1)less=1;
   else if(tl->i2<tr->i2)less=-1;
   else if(tl->i2>tr->i2)less=1;
   else less=0;

  return less;
 }

static MFVertex MFPivotTriangle(MF2dCellComplex C, MFVertex v,MFVertex b0, MFVertex b1)
 {
  int i,j;
  int a,b,c,d;
  MFEdge E,e;
  int i0,i1,i2;
  MFVertex v1;

/* Assumes that the triangle list is sorted */

  a=v->iPoly;
  b=b0->iPoly;
  c=b1->iPoly;

/*printf("Pivot: edge b0,b1 is (%d,%d), v=%d\n",b0->iPoly,b1->iPoly,v->iPoly);fflush(stdout);*/

  for(i=0;i<b1->nCobnd;i++)
   {
    E=(b1->e)[i];

    d=E->io;v1=E->o;
    if(d==c){d=E->id;v1=E->d;}

    if(d!=c && d!=b && d!=a)
     {
      if(b<d)e=MFFindEdge(C,b0,v1);
       else  e=MFFindEdge(C,v1,b0);

      if(e!=NULL)
       {
/*      printf("Pivot: edge (b0,d) is E[%d], edge (b1,d) is E[%d]\n",E->iEdge,e->iEdge);fflush(stdout);*/
        return v1;
       }
     }
   }

  return NULL;
 }

static void SortLists(MF2dCellComplex C)
 {
  int i;

  qsort(C->vertexList,C->nVertices,sizeof(MFVertex),VertexOrder);
  qsort(C->edgeList,C->nEdges,sizeof(MFEdge),EdgeOrder);
  qsort(C->triangleList,C->nTriangles,sizeof(MFTriangle),TriangleOrder);

  for(i=0;i<C->nVertices;i++) C->vertexList[i]->inverse  =i;
  for(i=0;i<C->nEdges;i++)    C->edgeList[i]->inverse    =i;
  for(i=0;i<C->nTriangles;i++)C->triangleList[i]->inverse=i;

  return;
 }

MF2dCellComplex MFPlotfileDual(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPlotfileDual"};
  int nC;

  MF2dCellComplex C;

  MFTriangle *T=NULL;
  MFEdge *E=NULL;
  MFVertex V=NULL;

  int i,j,l,n,m,k;
  int iPoly,nv,ne,nf;
  int iv,ie,iE0,je,np,v0,v1;
  int bnd,sing,v,ic,nIndices;
  int f,npoly;
  double c;
  int l0,l1,l2;
  int r0,r1,r2;
  int nEdges0;

  MFVertex va,vb,vc,vd;
  MFEdge ab,bd,cd,ac,ad,bc;
  int  nab,nbd,ncd,nac,nad,nbc;
  int adb,acb,bcd,acd;
  int state=-1;

  MFKVector Center;
  double *center;
  int vertex[100][500];
  double *(vertexC[100]);
  double vertexR[100];
  int edge[100][4];
  int face[100][2];
  MFVertex   VP;
  MFEdge     EP[100];
  MFTriangle TP[100];
  double R;
  double *b0;
  double *b1;
  double t0,t1;
  MFTriangle Ti;
  MFTriangle Tj;
  MFTriangle Tk;
  int iv0,iv1,iv2,iv3;
  MFEdge tmpedge;
  MFEdge *tmpEdgePtr;
  struct MFEdgeSt TMPEDGE;
  MFTriangle tmptri;
  MFTriangle *tmpTriPtr;
  struct MFTriangleSt TMPTRI;
  int verbose=0;

  C=(MF2dCellComplex)malloc(sizeof(struct MF2dCellComplexSt));

  C->nVertices=0;
  C->mVertices=0;
  C->vertexList=NULL;
  C->nEdges=0;
  C->mEdges=0;
  C->edgeList=NULL;
  C->nTriangles=0;
  C->mTriangles=0;
  C->triangleList=NULL;

  fscanf(fid,"Dimension of vertices, %d\n",&nC);
  fscanf(fid,"Dimension of manifold, %d\n",&k);
  if(verbose){printf("Vertices are in IR^%d, M is dimension %d\n",nC,k);fflush(stdout);}

  if(k!=2)
   {
    fprintf(stderr,"Error: %s, k(%2) must be 2\n",RoutineName,k);fflush(stdout);
    return NULL;
   }

  Center=MFCreateKVector(nC,e);
  b0=(double*)malloc(nC*sizeof(double));
  if(b0==NULL)
   {
    fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
    exit(12);
   }
  b1=(double*)malloc(nC*sizeof(double));
  if(b1==NULL)
   {
    fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
    exit(12);
   }
  center=MFKV_CStar(Center,e);
  for(i=0;i<20;i++)
   {
    vertexC[i]=(double*)malloc(nC*sizeof(double));
    if(vertexC[i]==NULL)
     {
      fprintf(stderr,"Out of memory, line %d in file %s\n",__LINE__,__FILE__);
      exit(12);
     }
   }

  while(!feof(fid))
   {

/* Input Example:
 *
 *     Polyhedron 6, 7 vertices, 6 edges, 6 faces, boundary 0, singular 0
 *     Vertex 0 (0.196000,0.969965,-0.155275), 2 [105,113]
 *     Vertex 1 (0.084337,0.989945,-0.127861), 2 [34,53]
 *     Vertex 2 (0.127175,0.989930,-0.085236), 2 [34,45]
 *     Vertex 3 (0.096403,0.980696,-0.180078), 2 [53,123]
 *     Vertex 4 (0.180886,0.980284,-0.098668), 2 [45,113]
 *     Vertex 5 (0.154125,0.970009,-0.196737), 2 [105,123]
 *     Vertex 6 (0.140710,0.979950,-0.141065), 0 [ ]
 *     Edge 0 (0,4), 1 [113]
 *     Edge 1 (0,5), 1 [105]
 *     Edge 2 (1,2), 1 [34]
 *     Edge 3 (1,3), 1 [53]
 *     Edge 4 (2,4), 1 [45]
 *     Edge 5 (3,5), 1 [123]
 *     Face 105 neighbor 15
 *     Face 53 neighbor 8
 *     Face 34 neighbor 1
 *     Face 45 neighbor 7
 *     Face 113 neighbor 16
 *     Face 123 neighbor 17
 */

    fscanf(fid,"Polyhedron %d, R=%lf, %d vertices, %d edges, %d faces,  boundary %d, singular %d\n",&iPoly,&R,&nv,&ne,&nf,&bnd,&sing);
    if(verbose){printf("Polyhedron %d, %d vertices, %d edges, %d faces,  boundary %d, singular %d\n\n",iPoly,nv,ne,nf,bnd,sing);fflush(stdout);}

    for(iv=0;iv<nv-1;iv++)
     {
      fscanf(fid,"Vertex %*d (%lf",&c);
      vertexC[iv][0]=c;
      if(iv==0)b0[0]=-c;
      if(iv==1){b0[0]+=c;b1[0]=-c;}
      if(iv==2)b1[0]+=c;
      if(verbose){printf("Vertex %d (%10.7lf",iv,c);fflush(stdout);}
      for(ic=1;ic<nC;ic++)
       {
        fscanf(fid,",%lf",&c);
        vertexC[iv][ic]=c;
        if(iv==0)b0[ic]=-c;
        if(iv==1){b0[ic]+=c;b1[ic]=-c;}
        if(iv==2)b1[ic]+=c;
        if(verbose){printf(",%10.7lf",c);fflush(stdout);}
       }
      fscanf(fid,"), %d [%d",&nIndices,&f);
      vertex[iv][0]=nIndices;
      vertex[iv][1]=f;
      if(verbose)printf(") [%d",vertex[iv][1]);

      for(ic=1;ic<nIndices;ic++)
       {
        fscanf(fid,",%d",&f);
        vertex[iv][ic+1]=f;
        if(verbose)printf(",%d",vertex[iv][ic+1]);
       }
      fscanf(fid,"]\n");
      if(verbose)printf("]\n");
     }
    if(verbose)printf("\n");

    t0=0.;
    for(i=0;i<nC;i++)t0+=b0[i]*b0[i];
    t0=sqrt(t0);
    t1=0.;
    for(i=0;i<nC;i++){b0[i]=b0[i]/t0;t1+=b0[i]*b1[i];}
    t0=0.;
    for(i=0;i<nC;i++){b1[i]=b1[i]-t1*b0[i];t0+=b1[i]*b1[i];}
    t0=sqrt(t0);
    for(i=0;i<nC;i++)b1[i]=b1[i]/t0;

/* "The" Vertex */

    fscanf(fid,"Vertex %*d (%lf",&c);
    if(verbose)printf("Poly %d, center (%lf",iPoly,c);
    center[0]=c;
    for(i=1;i<nC;i++)
     {
      fscanf(fid,",%lf",&c);
      center[i]=c;
      if(verbose)printf(",%lf",c);
     }
    fscanf(fid,"), %d [ ] \n",&f);
    if(verbose)printf("), %d [ ] \n",f);

    for(j=0;j<nv-1;j++)
     {
      t0=0;t1=0.;
      for(i=0;i<nC;i++){t0+=(vertexC[j][i]-center[i])*b0[i];t1+=(vertexC[j][i]-center[i])*b1[i];}
      vertexR[j]=sqrt(t0*t0+t1*t1);
     }

/* The Edges */

    for(ie=0;ie<ne;ie++)
     {
      fscanf(fid,"Edge %*d (%d,%d), %*d [%d]\n",&(edge[ie][0]),&(edge[ie][1]),&(edge[ie][2]));
      if(verbose)printf("Edge %5d (%5d,%5d), 1 [%5d]\n",ie,edge[ie][0],edge[ie][1],edge[ie][2]);fflush(stdout);
     }
    if(verbose)printf("\n");

/* The Faces */

    for(ie=0;ie<ne;ie++)
     {
      fscanf(fid,"Face %d neighbor %d\n",face[ie]+0,face[ie]+1);
      if(verbose)printf("Face %5d neighbor %5d\n",face[ie][0],face[ie][1]);
     }
    if(verbose)printf("\n");
    fflush(stdout);

/* Order the edges for the co-boundary of the vertex */

    iv=edge[0][1];

    for(f=1;f<ne;f++)
     {
      ie=f;
      while( ie<ne && edge[ie][0]!=iv && edge[ie][1]!=iv)ie++;

      if(edge[ie][1]==iv)
       {
        i=edge[ie][0];
        edge[ie][0]=edge[ie][1];
        edge[ie][1]=i;
       }

      i=edge[f][0];
      edge[f][0]=edge[ie][0];
      edge[ie][0]=i;

      i=edge[f][1];
      edge[f][1]=edge[ie][1];
      edge[ie][1]=i;

      i=edge[f][2];
      edge[f][2]=edge[ie][2];
      edge[ie][2]=i;

      iv=edge[f][1];
     }

    for(ie=0;ie<ne;ie++)
      if(verbose)printf("Edge %5d (%5d,%5d), 1 [%5d]\n",ie,edge[ie][0],edge[ie][1],edge[ie][2]);
    if(verbose){printf("\n");fflush(stdout);}

/* Translate the face indices of the edges into the neighboring chart index. */

    for(ie=0;ie<ne;ie++)
     {
      f=0;
      while(edge[ie][2]!=face[f][0])f++;
      edge[ie][2]=iPoly;
      edge[ie][3]=face[f][1];
     }

    if(verbose)printf("\nAfter translating the edge duals\n");
    for(ie=0;ie<ne;ie++)
      if(verbose)printf("Edge %5d (%5d,%5d), dual: [%5d,%5d]\n",ie,edge[ie][0],edge[ie][1],edge[ie][2],edge[ie][3]);
    if(verbose){printf("\n");fflush(stdout);}

/* Translate the face indices of the vertices into the neighboring chart index. */

    for(iv=0;iv<ne;iv++)
     {
      for(ic=0;ic<vertex[iv][0];ic++)
       {
        f=0;while(vertex[iv][ic+1]!=face[f][0])f++;
        vertex[iv][ic+1]=face[f][1];
       }
      vertex[iv][vertex[iv][0]+1]=iPoly;
      vertex[iv][0]++;

/* Sort */

      for(ic=0;ic<vertex[iv][0];ic++)
       {
        for(i=ic+1;i<vertex[iv][0];i++)
         {
          if(vertex[iv][i+1]<vertex[iv][ic+1])
           {
            f=vertex[iv][i+1];
            vertex[iv][i+1]=vertex[iv][ic+1];
            vertex[iv][ic+1]=f;
           }
         }
       }
     }

    if(verbose)
     {
      printf("\nAfter translating the vertex duals, ne=%d\n",ne);
      for(iv=0;iv<ne;iv++)
       {
        if(verbose)
         {
          printf("Vertex %5d, dual: [%5d",iv,vertex[iv][1]);
          for(ic=1;ic<vertex[iv][0];ic++)
           {
            printf(",%5d",vertex[iv][ic+1]);fflush(stdout);
           }
          printf("]\n");
         }
       }
      printf("Add Vertex\n");fflush(stdout);
     }

    VP=MFAddVertex(C,iPoly,nC,center);
    if(verbose)PrintVertex(iPoly,VP);

    for(i=0;i<ne;i++)
     {
      j=i-1;
      if(j<0)j+=ne;
      l=i+1;
      if(l>ne-1)l-=ne;

      l0=edge[i][2];
      l1=edge[i][3];
      l2=edge[j][3];

      r0=edge[i][2];
      r1=edge[i][3];
      r2=edge[l][3];

      if(0&&verbose){printf("edge %d, prev=(%d,%d) e=(%d,%d) next=(%d,%d)",i,edge[j][2],edge[j][3],edge[i][2],edge[i][3],edge[l][2],edge[l][3],r1);fflush(stdout);}

      if(vertexR[edge[i][0]]<=R || vertexR[edge[i][1]]<=R )
       {
        if(l1<l0){l=l0;l0=l1;l1=l;}
        if(l2<l0){l=l0;l0=l2;l2=l;}
        if(l2<l1){l=l1;l1=l2;l2=l;}

        TP[i]=MFAddTriangle(C,l0,l1,l2);

        if(r1<r0){l=r0;r0=r1;r1=l;}
        if(r2<r0){l=r0;r0=r2;r2=l;}
        if(r2<r1){l=r1;r1=r2;r2=l;}

        if(0&&verbose){printf("  endpoints: [(%5d,%5d,%5d),(%5d,%5d,%5d)]\n",l0,l1,l2,r0,r1,r2);fflush(stdout);}
        if(verbose){printf("Add Triangle ");PrintTriangle(TP[i]->iTri,TP[i]);fflush(stdout);}

        if(edge[i][2]<edge[i][3])EP[i]=MFAddEdge(C,edge[i][2],edge[i][3],l0,l1,l2,r0,r1,r2);
         else                    EP[i]=MFAddEdge(C,edge[i][3],edge[i][2],l0,l1,l2,r0,r1,r2);
        if(verbose){printf("Add Edge     ");PrintEdge(EP[i]->iEdge,EP[i]);fflush(stdout);}
       }else{
        if(verbose){printf("\n");fflush(stdout);}
       }
     }

    for(i=0;i<ne;i++)
     {
      j=i-1;
      if(j<0)j+=ne;
      l0=edge[i][2];
      l1=edge[i][3];
      l2=edge[j][3];

      (TP[i])->e0=EP[i];
      if(l0<l1)(TP[i])->pm0=1.;
        else   (TP[i])->pm0=-1.;

      (TP[i])->e1=NULL;
      (TP[i])->pm1=1.;
      if(l0<l2)(TP[i])->pm1=1.;
        else   (TP[i])->pm1=-1.;

      (TP[i])->e2=EP[j];
      (TP[i])->pm2=-1.;
     }

    if(verbose){printf("\n");fflush(stdout);}
   }

  SortLists(C);

  if(verbose){printf("\nsorted but not compressed\n");PrintItAll(C);}

/* Fill in the third edge of the triangles. */

  for(i=0;i<C->nTriangles;i++)
   {
    if(((C->triangleList[i])->e0)->io==((C->triangleList[i])->e2)->io)
            {va=C->vertexList[((C->triangleList[i])->e0)->id];vb=C->vertexList[((C->triangleList[i])->e2)->id];vc=C->vertexList[((C->triangleList[i])->e2)->io];}
     else if(((C->triangleList[i])->e0)->id==((C->triangleList[i])->e2)->io)
            {va=C->vertexList[((C->triangleList[i])->e0)->io];vb=C->vertexList[((C->triangleList[i])->e2)->id];vc=C->vertexList[((C->triangleList[i])->e2)->io];}
     else if(((C->triangleList[i])->e0)->io==((C->triangleList[i])->e2)->id)
            {va=C->vertexList[((C->triangleList[i])->e0)->id];vb=C->vertexList[((C->triangleList[i])->e2)->io];vc=C->vertexList[((C->triangleList[i])->e2)->id];}
     else                                                             
            {va=C->vertexList[((C->triangleList[i])->e0)->io];vb=C->vertexList[((C->triangleList[i])->e2)->io];vc=C->vertexList[((C->triangleList[i])->e2)->id];}
   
    (C->triangleList[i])->e1=MFFindEdge(C,va,vb);
    if( (C->triangleList[i])->e1!=NULL)
     {
      if( ((C->triangleList[i])->e1)->io==va->inverse)(C->triangleList[i])->pm1=-1;
       else                                        (C->triangleList[i])->pm1= 1;
     }else{
      printf("Couldn't complete triangle %d, adding edge\n",i);

      l0=va->iPoly;l1=vb->iPoly;l2=vc->iPoly;
      if(l1<l0){iv=l0;l0=l1;l1=iv;}
      if(l2<l0){iv=l0;l0=l2;l2=iv;}
      if(l1<l0){iv=l0;l0=l1;l1=iv;}

      (C->triangleList[i])->e1=MFAddEdge(C,va->iPoly,vb->iPoly,l0,l1,l2,r0,r1,r2);

      SortLists(C);

      if( ((C->triangleList[i])->e1)->io==va->inverse)(C->triangleList[i])->pm1=-1;
       else                                        (C->triangleList[i])->pm1= 1;

      PrintTriangle(i,C->triangleList[i]);fflush(stdout);
     }
   }

/* Compress the edge list */

  nEdges0=C->nEdges;
  i=0;
  j=1;
  while(j<C->nEdges)
   {
    while(j<C->nEdges && (C->edgeList[i])->io==(C->edgeList[j])->io && (C->edgeList[i])->id==(C->edgeList[j])->id
                    && ( ( (C->edgeList[i])->l0==(C->edgeList[j])->l0 && (C->edgeList[i])->l1==(C->edgeList[j])->l1 && (C->edgeList[i])->l2==(C->edgeList[j])->l2
                        && (C->edgeList[i])->r0==(C->edgeList[j])->r0 || (C->edgeList[i])->r1==(C->edgeList[j])->r1 || (C->edgeList[i])->r2==(C->edgeList[j])->r2) 
                    || ( (C->edgeList[i])->l0==(C->edgeList[j])->r0 && (C->edgeList[i])->l1==(C->edgeList[j])->r1 && (C->edgeList[i])->l2==(C->edgeList[j])->r2
                        && (C->edgeList[i])->r0==(C->edgeList[j])->l0 && (C->edgeList[i])->r1==(C->edgeList[j])->l1 && (C->edgeList[i])->r2==(C->edgeList[j])->l2) ) )j++;
    i++;
    if(i!=j && j<C->nEdges)
     {
      (C->edgeList[i])->iEdge=(C->edgeList[j])->iEdge;
      (C->edgeList[i])->io=(C->edgeList[j])->io;
      (C->edgeList[i])->id=(C->edgeList[j])->id;
      (C->edgeList[i])->l0=(C->edgeList[j])->l0;
      (C->edgeList[i])->l1=(C->edgeList[j])->l1;
      (C->edgeList[i])->l2=(C->edgeList[j])->l2;
      (C->edgeList[i])->r0=(C->edgeList[j])->r0;
      (C->edgeList[i])->r1=(C->edgeList[j])->r1;
      (C->edgeList[i])->r2=(C->edgeList[j])->r2;
      (C->edgeList[i])->inverse=i;
      (C->edgeList[i])->counter=(C->edgeList[j])->counter;
     }
    j++;
   }
  C->nEdges=i;
  while(C->nEdges>0 && (C->edgeList[C->nEdges-1])->io==(C->edgeList[C->nEdges-2])->io
                 && (C->edgeList[C->nEdges-1])->id==(C->edgeList[C->nEdges-2])->id
                 && ( ( (C->edgeList[C->nEdges-1])->l0==(C->edgeList[C->nEdges-2])->l0
                     && (C->edgeList[C->nEdges-1])->l1==(C->edgeList[C->nEdges-2])->l1
                     && (C->edgeList[C->nEdges-1])->l2==(C->edgeList[C->nEdges-2])->l2
                 && (C->edgeList[C->nEdges-1])->r0==(C->edgeList[C->nEdges-2])->r0
                 || (C->edgeList[C->nEdges-1])->r1==(C->edgeList[C->nEdges-2])->r1
                 || (C->edgeList[C->nEdges-1])->r2==(C->edgeList[C->nEdges-2])->r2) 
                 || ( (C->edgeList[C->nEdges-1])->l0==(C->edgeList[C->nEdges-2])->r0
                   && (C->edgeList[C->nEdges-1])->l1==(C->edgeList[C->nEdges-2])->r1
                   && (C->edgeList[C->nEdges-1])->l2==(C->edgeList[C->nEdges-2])->r2
                 && (C->edgeList[C->nEdges-1])->r0==(C->edgeList[C->nEdges-2])->l0
                 && (C->edgeList[C->nEdges-1])->r1==(C->edgeList[C->nEdges-2])->l1
                 && (C->edgeList[C->nEdges-1])->r2==(C->edgeList[C->nEdges-2])->l2) ) )C->nEdges--;

  for(i=0;i<C->nVertices;i++) C->vertexList[i]->inverse  =i;
  for(i=0;i<C->nEdges;i++)    C->edgeList[i]->inverse    =i;
  for(i=0;i<C->nTriangles;i++)C->triangleList[i]->inverse=i;

  for(i=0;i<C->nEdges;i++)
   {
    C->edgeList[i]->o=MFFindVertex(C,C->edgeList[i]->io);
    C->edgeList[i]->d=MFFindVertex(C,C->edgeList[i]->id);

    VertexAddCoboundary((C->edgeList[i])->o,-1,C->edgeList[i]);
    VertexAddCoboundary((C->edgeList[i])->d, 1,C->edgeList[i]);
   }
  if(verbose){printf("sorted and compressed\n");PrintItAll(C);}

  tmpedge=&TMPEDGE;
  tmptri=&TMPTRI;
  
  for(i=0;i<C->nTriangles;i++)
   {
    l0=(C->triangleList[i])->i0;
    l1=(C->triangleList[i])->i1;
    l2=(C->triangleList[i])->i2;

    tmpedge->io=l0;
    tmpedge->id=l1;
    tmpEdgePtr=(MFEdge*)bsearch(&tmpedge,C->edgeList,C->nEdges,sizeof(MFEdge),EdgeOrder);
    if(tmpEdgePtr==NULL)
     {
      (C->triangleList[i])->e0=NULL;
     }else{
      (C->triangleList[i])->e0=*tmpEdgePtr;
      ((*tmpEdgePtr)->counter)++;
     }
    (C->triangleList[i])->pm0=1;

    tmpedge->io=l1;
    tmpedge->id=l2;
    tmpEdgePtr=(MFEdge*)bsearch(&tmpedge,C->edgeList,C->nEdges,sizeof(MFEdge),EdgeOrder);
    if(tmpEdgePtr==NULL)
     {
      (C->triangleList[i])->e1=NULL;
     }else{
      (C->triangleList[i])->e1=*tmpEdgePtr;
      ((*tmpEdgePtr)->counter)++;
     }
    (C->triangleList[i])->pm1=1;

    tmpedge->io=l0;
    tmpedge->id=l2;
    tmpEdgePtr=(MFEdge*)bsearch(&tmpedge,C->edgeList,C->nEdges,sizeof(MFEdge),EdgeOrder);
    if(tmpEdgePtr==NULL)
     {
      (C->triangleList[i])->e2=NULL;
     }else{
      (C->triangleList[i])->e2=*tmpEdgePtr;
      ((*tmpEdgePtr)->counter)++;
     }
    (C->triangleList[i])->pm2=-1;
   }

  for(i=0;i<C->nEdges;i++)
   {
    if(C->edgeList[i]->counter==0)
     {
      MFRemoveEdgeByIndex(C,i);
      i--;
     }
   }

  i=0;
  j=1;
  while(j<C->nEdges)
   {
    while(j<C->nEdges && C->edgeList[i]->counter==0)j++;
    i++;
    if(i!=j && j<C->nEdges)
     {
      (C->edgeList[i])->iEdge=(C->edgeList[j])->iEdge;
      (C->edgeList[i])->io=(C->edgeList[j])->io;
      (C->edgeList[i])->id=(C->edgeList[j])->id;
      (C->edgeList[i])->l0=(C->edgeList[j])->l0;
      (C->edgeList[i])->l1=(C->edgeList[j])->l1;
      (C->edgeList[i])->l2=(C->edgeList[j])->l2;
      (C->edgeList[i])->r0=(C->edgeList[j])->r0;
      (C->edgeList[i])->r1=(C->edgeList[j])->r1;
      (C->edgeList[i])->r2=(C->edgeList[j])->r2;
      (C->edgeList[i])->inverse=i;
      (C->edgeList[i])->counter=(C->edgeList[j])->counter;
     }
    j++;
   }
  C->nEdges=i;
  while(C->nEdges>0 && C->edgeList[C->nEdges-1]->counter==0)C->nEdges--;

  if(verbose){printf("sorted, compressed and triangle edges filled in\n");PrintItAll(C);}

  if(1||verbose)
   {
    printf("\nSuspect edges: (counts of # triangles each edge is in)\n\n");
    for(i=0;i<C->nEdges;i++)
     {
      if((C->edgeList[i])->counter!=6)
       {
        printf(" %5d E[%5d] (%5d,%5d) [(%5d,%5d,%5d),(%5d,%5d,%5d)] appears %5d times\n",
            i,C->edgeList[i]->iEdge,
            C->edgeList[i]->io,C->edgeList[i]->id,
            C->edgeList[i]->l0,C->edgeList[i]->l1,C->edgeList[i]->l2,
            C->edgeList[i]->r0,C->edgeList[i]->r1,C->edgeList[i]->r2,(C->edgeList[i])->counter);
       }
     }
    printf("\n");fflush(stdout);
   }

  if(0){PrintItAll(C);printf("\n");fflush(stdout);}

  if(1||verbose){printf("Examine and correct suspect triangles\n\n");fflush(stdout);}

  for(i=0;i<C->nTriangles;i++)
   {
    l0=-1;
    l1=-1;
    l2=-1;
    if((C->triangleList[i])->e0!=NULL)
     {
      if( ((C->triangleList[i])->e0)->iEdge<0 || ((C->triangleList[i])->e0)->iEdge>nEdges0-1)
       {
        printf(" checking triangle %d\n",i);fflush(stdout);
        printf("illegal value for Triangle %d->e0->iEdge! %d !in [0,%d)\n",i,((C->triangleList[i])->e0)->iEdge,nEdges0);
        PrintItAll(C);fflush(stdout);
        exit(12);
       }
      l0=((C->triangleList[i])->e0)->counter;
      r0=((C->triangleList[i])->e0)->inverse;
     }
    if((C->triangleList[i])->e1!=NULL)
     {
      if( ((C->triangleList[i])->e1)->iEdge<0 || ((C->triangleList[i])->e1)->iEdge>nEdges0-1)
       {
        printf(" checking triangle %d\n",i);fflush(stdout);
        printf("illegal value for Triangle %d->e1->iEdge! %d !in [0,%d)\n",i,((C->triangleList[i])->e1)->iEdge,nEdges0);
        PrintItAll(C);fflush(stdout);
        exit(12);
       }
      l1=((C->triangleList[i])->e1)->counter;
      r1=((C->triangleList[i])->e1)->inverse;
     }
    if((C->triangleList[i])->e2!=NULL)
     {
      if( ((C->triangleList[i])->e2)->iEdge<0 || ((C->triangleList[i])->e2)->iEdge>nEdges0-1)
       {
        printf(" checking triangle %d\n",i);fflush(stdout);
        printf("illegal value for Triangle %d->e2->iEdge! %d !in [0,%d)\n",i,((C->triangleList[i])->e2)->iEdge,nEdges0);
        PrintItAll(C);fflush(stdout);
        exit(12);
       }
      l2=((C->triangleList[i])->e2)->counter;
      r2=((C->triangleList[i])->e2)->inverse;
     }
    if(l1<l0){j=l0;l0=l1;l1=j;j=r0;r0=r1;r1=j;}
    if(l2<l0){j=l0;l0=l2;l2=j;j=r0;r0=r2;r2=j;}
    if(l2<l1){j=l1;l1=l2;l2=j;j=r1;r1=r2;r2=j;}

    if(l0!=-1&&l1!=-1&&l2!=-1 && (l0!=6||l1!=6||l2!=6) )
     {
      va=NULL;
      vb=NULL;
      vc=NULL;
      vd=NULL;

      ab=NULL;
      ac=NULL;
      ad=NULL;
      bc=NULL;
      bd=NULL;
      cd=NULL;

      if(1||verbose){printf("%5d T[%5d] is suspect (%d,%d,%d)\n",i,(C->triangleList[i])->iTri,l0,l1,l2);fflush(stdout);}

      vb=(C->edgeList[r0])->o;
      vc=(C->edgeList[r0])->d;
      va=(C->edgeList[r1])->o;
      if(va->iPoly==vb->iPoly || va->iPoly==vc->iPoly)va=(C->edgeList[r1])->d;
      vd=MFPivotTriangle(C,va,vb,vc);

      ab=C->edgeList[r1];
      ac=C->edgeList[r2];
      ad=MFFindEdge(C,va,vd);
      bc=C->edgeList[r0];
      bd=MFFindEdge(C,vb,vd);
      cd=MFFindEdge(C,vc,vd);

      if(ad!=NULL)
       {
        nab=ab->counter-3;
        nac=ac->counter-3;
        nad=ad->counter;
        nbc=bc->counter;
        nbd=bd->counter-3;
        ncd=cd->counter-3;

        acb=(nab+nac-nad)/2;
        bcd=nbc-acb;
        adb=nbd-bcd;
        acd=ncd-bcd;
        if(nab==acb+adb && nac==acb+acd && nad==adb+acd && nbc==acb+bcd && nbd==bcd+adb && ncd==bcd+acd)
         {
          if(1||verbose)
           {
            printf(" acb=%d\n",acb);fflush(stdout);
            printf(" bcd=%d\n",bcd);fflush(stdout);
            printf(" adb=%d\n",adb);fflush(stdout);
            printf(" acd=%d\n\n",acd);fflush(stdout);
  
            printf("Check: ab: %d=%d\n",nab,acb+adb);fflush(stdout);
            printf("       ac: %d=%d\n",nac,acb+acd);fflush(stdout);
            printf("       ad: %d=%d\n",nad,adb+acd);fflush(stdout);
            printf("       bc: %d=%d\n",nbc,acb+bcd);fflush(stdout);
            printf("       bd: %d=%d\n",nbd,bcd+adb);fflush(stdout);
            printf("       cd: %d=%d\n",ncd,bcd+acd);fflush(stdout);
  
            printf("\n Before removing and adding triangles\n\n");fflush(stdout);
  
            printf(" va=");PrintVertex(-1,va);
            printf(" vb=");PrintVertex(-1,vb);
            printf(" vc=");PrintVertex(-1,vc);
            printf(" vd=");PrintVertex(-1,vd);
            printf("\n");fflush(stdout);
  
            if(ab!=NULL){printf("ab=");PrintEdge(-1,ab);fflush(stdout);}
            if(ac!=NULL){printf("ac=");PrintEdge(-1,ac);fflush(stdout);}
            if(ad!=NULL){printf("ad=");PrintEdge(-1,ad);fflush(stdout);}
            if(bc!=NULL){printf("bc=");PrintEdge(-1,bc);fflush(stdout);}
            if(bd!=NULL){printf("bd=");PrintEdge(-1,bd);fflush(stdout);}
            if(cd!=NULL){printf("cd=");PrintEdge(-1,cd);fflush(stdout);}
           }
  
          for(i=0;i<adb;i++)MFRemoveTriangle(C,va,vd,vb);
          for(i=0;i<acd;i++)MFRemoveTriangle(C,va,vc,vd);
          for(i=acb;i<3;i++)MFAddTriangleAndSetEdges(C,va,ab,vb,bc,vc,ac);
          for(i=bcd;i<3;i++)MFAddTriangleAndSetEdges(C,vb,bc,vc,cd,vd,bd);
  
          if(1||verbose)
           {
            printf("\n After removing and adding triangles\n\n");fflush(stdout);
            if(ab!=NULL){printf("ab=");PrintEdge(-1,ab);fflush(stdout);}
            if(ac!=NULL){printf("ac=");PrintEdge(-1,ac);fflush(stdout);}
            if(ad!=NULL){printf("ad=");PrintEdge(-1,ad);fflush(stdout);}
            if(bc!=NULL){printf("bc=");PrintEdge(-1,bc);fflush(stdout);}
            if(bd!=NULL){printf("bd=");PrintEdge(-1,bd);fflush(stdout);}
            if(cd!=NULL){printf("cd=");PrintEdge(-1,cd);fflush(stdout);}
            printf("\n");fflush(stdout);
           }

          MFRemoveEdge(C,va,vd);
          printf("done\n");
  
          SortLists(C);
          i=-1;
         }
       }
     }
   }

  if(1||verbose)
   {
    printf("\nSuspect edges after fixing them: (counts of # triangles each edge is in)\n\n");
    for(i=0;i<C->nEdges;i++)
     {
      if((C->edgeList[i])->counter!=6)
       {
        printf(" %5d E[%5d] (%5d,%5d) [(%5d,%5d,%5d),(%5d,%5d,%5d)] appears %5d times\n",
            i,C->edgeList[i]->iEdge,
            C->edgeList[i]->io,C->edgeList[i]->id,
            C->edgeList[i]->l0,C->edgeList[i]->l1,C->edgeList[i]->l2,
            C->edgeList[i]->r0,C->edgeList[i]->r1,C->edgeList[i]->r2,(C->edgeList[i])->counter);
       }
     }
    printf("\n");fflush(stdout);
   }

  if(verbose){printf("Before Compressing C->triangleList\n");PrintItAll(C);fflush(stdout);}

/* Compress the Triangle list */

  i=0;
  j=1;
  while(j<C->nTriangles)
   {
    while(j<C->nTriangles && (C->triangleList[i])->i0==(C->triangleList[j])->i0 && (C->triangleList[i])->i1==(C->triangleList[j])->i1 && (C->triangleList[i])->i2==(C->triangleList[j])->i2)
     {
      if((C->triangleList[j])->e0!=NULL)(((C->triangleList[j])->e0)->counter)--;
      if((C->triangleList[j])->e1!=NULL)(((C->triangleList[j])->e1)->counter)--;
      if((C->triangleList[j])->e2!=NULL)(((C->triangleList[j])->e2)->counter)--;
      j++;
     }
    i++;
    if(i!=j && j<C->nTriangles)
     {
      (C->triangleList[i])->iTri=(C->triangleList[j])->iTri;
      (C->triangleList[i])->i0=(C->triangleList[j])->i0;
      (C->triangleList[i])->i1=(C->triangleList[j])->i1;
      (C->triangleList[i])->i2=(C->triangleList[j])->i2;

      (C->triangleList[i])->e0 =(C->triangleList[j])->e0;
      (C->triangleList[i])->pm0=(C->triangleList[j])->pm0;
      (C->triangleList[i])->e1 =(C->triangleList[j])->e1;
      (C->triangleList[i])->pm1=(C->triangleList[j])->pm1;
      (C->triangleList[i])->e2 =(C->triangleList[j])->e2;
      (C->triangleList[i])->pm2=(C->triangleList[j])->pm2;
      (C->triangleList[i])->inverse=i;

      if((C->triangleList[i])->e0!=NULL)
       {
        if((C->triangleList[i])->pm0>0)((C->triangleList[i])->e0)->r=C->triangleList[i];
          else                      ((C->triangleList[i])->e0)->l=C->triangleList[i];
       }

      if((C->triangleList[i])->e1!=NULL)
       {
        if((C->triangleList[i])->pm1>0)((C->triangleList[i])->e1)->r=C->triangleList[i];
          else                      ((C->triangleList[i])->e1)->l=C->triangleList[i];
       }

      if((C->triangleList[i])->e2!=NULL)
       {
        if((C->triangleList[i])->pm2>0)((C->triangleList[i])->e2)->r=C->triangleList[i];
          else                      ((C->triangleList[i])->e2)->l=C->triangleList[i];
       }
     }
    j++;
   }
  C->nTriangles=i;

  if(verbose){printf("After Compressing C->triangleList\n");PrintItAll(C);fflush(stdout);}

/* Choose orientations of the dual edges so that the vertex has a closed coboundary */

  return C;
 }

static void PrintItAll(MF2dCellComplex C )
 {
  int i,j;

  printf("\nThe vertex list:\n\n");
  for(i=0;i<C->nVertices;i++)PrintVertex(i,C->vertexList[i]);

  printf("\nThe edge list:\n\n");
  for(i=0;i<C->nEdges;i++)PrintEdge(i,C->edgeList[i]);

  printf("\nThe triangle list:\n\n");
  for(i=0;i<C->nTriangles;i++)PrintTriangle(i,C->triangleList[i]);

  return;
 }

static void PrintVertex(int i, MFVertex v)
 {
  int j;

  if(i>-1)printf(" %5d V[%5d]: (%5d) (%10.7lf",i,v->iPoly,v->iPoly,v->center[0]);
   else   printf("V[%5d]: (%5d) (%10.7lf",v->iPoly,v->iPoly,v->center[0]);

  for(j=1;j<v->nC;j++)printf(",%10.7lf",v->center[j]);
  printf("), cobnd=");fflush(stdout);
  for(j=0;j<v->nCobnd;j++)
   {
    if(v->pm[j]>0&&j>0)printf(" + ");
     else if(v->pm[j]>0)printf("   ");
     else printf(" - ");
    printf("E[%5d]",((v->e)[j])->iEdge);
   }
  if(v->nCobnd==0)printf(" 0 ");
  printf("\n");
  fflush(stdout);

  return;
 }

static void PrintEdge(int i, MFEdge e)
 {
  char i0[10],i1[10];
  char j0[10],j1[10];

  strcpy(i0,"**");
  strcpy(i1,"**");
  strcpy(j0,"**");
  strcpy(j1,"**");

  if(e->o!=NULL)sprintf(i0,"%5d",(e->o)->iPoly);
   else sprintf(i0,"    ?");
  if(e->d!=NULL)sprintf(i1,"%5d",(e->d)->iPoly);
   else sprintf(i1,"    ?");
  if(e->l!=NULL)sprintf(j0,"%5d",(e->l)->iTri);
   else sprintf(j0,"    ?");
  if(e->r!=NULL)sprintf(j1,"%5d",(e->r)->iTri);
   else sprintf(j1,"    ?");

  if(i>-1)printf(" %5d E[%5d] (%5d,%5d) [(%5d,%5d,%5d),(%5d,%5d,%5d)] bnd=",i,e->iEdge,e->io,e->id,e->l0,e->l1,e->l2,e->r0,e->r1,e->r2);
   else   printf("E[%5d] (%5d,%5d) [(%5d,%5d,%5d),(%5d,%5d,%5d)] bnd=",e->iEdge,e->io,e->id,e->l0,e->l1,e->l2,e->r0,e->r1,e->r2);

  if(e->pmd<0)printf(" - ");
  printf("V[%s]",i1);
  if(e->pmo>0)printf(" + ");
  if(e->pmo<0)printf(" - ");
  printf("V[%s]",i0);
  printf(" cobnd=");
  if(e->pmr<0)printf(" - ");
  printf("T[%s]",j1);
  if(e->pml>0)printf(" + ");
  if(e->pml<0)printf(" - ");
  printf("T[%s]",j0);
  printf(" %d\n",e->counter);
  fflush(stdout);

  return;
 }


static void PrintTriangle(int i, MFTriangle t)
 {
  char io[10],id[10];
  char jo[10],jd[10];
  char ko[10],kd[10];

  strcpy(io,"**");
  strcpy(id,"**");
  strcpy(jo,"**");
  strcpy(jd,"**");
  strcpy(ko,"**");
  strcpy(kd,"**");

  if(t->e0!=NULL)
   {
    sprintf(io,"%5d",(t->e0)->io);
    sprintf(id,"%5d",(t->e0)->id);
   }else{
    sprintf(io,"    ?");
    sprintf(id,"    ?");
   }
  if(t->e1!=NULL)
   {
    sprintf(jo,"%5d",(t->e1)->io);
    sprintf(jd,"%5d",(t->e1)->id);
   }else{
    sprintf(jo,"    ?");
    sprintf(jd,"    ?");
   }
  if(t->e2!=NULL)
   {
    sprintf(ko,"%5d",(t->e2)->io);
    sprintf(kd,"%5d",(t->e2)->id);
   }else{
    sprintf(ko,"    ?");
    sprintf(kd,"    ?");
   }

  printf(" %5d T[%5d] (%5d,%5d,%5d) bnd =",i,t->iTri,t->i0,t->i1,t->i2);
  fflush(stdout);

  if(t->pm0==-1)printf(" - ");
   else printf("   ");
  if(t->e0!=NULL)printf("E[%5d]",(t->e0)->iEdge);
    else printf("E[    ?]");

  if(t->pm1==-1)printf(" - ");
   else printf(" + ");
  if(t->e1!=NULL)printf("E[%5d]",(t->e1)->iEdge);
    else printf("E[    ?]");

  if(t->pm2==-1)printf(" - ");
   else printf(" + ");
  if(t->e2!=NULL)printf("E[%5d]",(t->e2)->iEdge);
    else printf("E[    ?]");

  printf(" bnd bnd = ");
  fflush(stdout);

  if(t->pm0==-1)printf("V[%s] - V[%s]",io,id);
    else printf("V[%s] - V[%s]",id,io);

  if(t->pm1==-1)printf(" + V[%s] - V[%s]",jo,jd);
    else printf(" + V[%s] - V[%s]",jd,jo);

  if(t->pm2==-1)printf(" + V[%s] - V[%s]",ko,kd);
    else printf(" + V[%s] - V[%s]",kd,ko);

  printf("\n");
  fflush(stdout);

  return;
 }

void MFOutput2dCellComplexAsEasyMesh(MF2dCellComplex C, char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"MFOutput2dCellComplexAsEasyMesh"};
  int i,j,k;
  MFVertex V;
  double c;

  MFEdge     EP[100];
  MFTriangle TP[100];

  char fullname[1000];
  FILE *fout;

/* Now write three files
 *   nodes: name.n      
 *   sides: name.s
 *   elements: name.e
 */

/* nodes: name.n
 *
 *     # first line:       <number of nodes>
 *     # following lines:  <node number:> <x> <y> <marker>
 */

  printf("%s\n",RoutineName);fflush(stdout);

  strcpy(fullname,name);
  strcat(fullname,"Dual.n");
  fout=fopen(fullname,"w");
  fprintf(fout,"%d\n",C->nVertices);
  for(i=0;i<C->nVertices;i++)
   {
    V=C->vertexList[i];
    fprintf(fout,"%d:",i);
    for(j=0;j<V->nC;j++)fprintf(fout," %lf",(V->center)[j]);
    fprintf(fout," 0\n");
   }
  fclose(fout);

/* edges (sides): name.s
 *
 *    # first line:       <number of sides>
 *    # following lines:  <side number:> <c> <d> <ea> <eb> <marker> 
 */

  strcpy(fullname,name);
  strcat(fullname,"Dual.s");
  fout=fopen(fullname,"w");

  fprintf(fout,"%d\n",C->nEdges);

  for(i=0;i<C->nEdges;i++)
   {
    EP[0]=C->edgeList[i];
    if(EP[0]->iEdge!=-1)
     {
      fprintf(fout,"%5d: %5d %5d",i,(EP[0])->io,(EP[0])->id);
      if(EP[0]->l!=NULL)fprintf(fout," %5d",(EP[0]->l)->inverse);
       else             fprintf(fout," -1");
      if(EP[0]->r!=NULL)fprintf(fout," %5d",(EP[0]->r)->inverse);
       else             fprintf(fout," -1");
      fprintf(fout," 0\n");
     }
   }
  fclose(fout);

/* elements: name.e
 *
 *     # first line        <number of elements>
 *     # following lines:  <element number:> <i> <j> <k> <ei> <ej> <ek> <si> <sj> <sk> <xV> <yV> <marker>
 */

  strcpy(fullname,name);
  strcat(fullname,"Dual.e");
  fout=fopen(fullname,"w");

  fprintf(fout,"%5d\n",C->nTriangles);

  for(i=0;i<C->nTriangles;i++)
   {
    TP[0]=C->triangleList[i];
    if(TP[0]->iTri!=-1)
     {
      fprintf(fout,"%5d: %5d %5d %5d",i,(TP[0])->i0,(TP[0])->i1,(TP[0])->i2);

      if((TP[0])->e0!=NULL && (TP[0])->e1!=NULL && (TP[0])->e2!=NULL)
       {
             if(((TP[0])->e0)->io!=(TP[0])->i0 && ((TP[0])->e0)->id!=(TP[0])->i0)EP[0]=(TP[0])->e0;
        else if(((TP[0])->e1)->io!=(TP[0])->i0 && ((TP[0])->e1)->id!=(TP[0])->i0)EP[0]=(TP[0])->e1;
        else if(((TP[0])->e2)->io!=(TP[0])->i0 && ((TP[0])->e2)->id!=(TP[0])->i0)EP[0]=(TP[0])->e2;
             if(((TP[0])->e0)->io!=(TP[0])->i1 && ((TP[0])->e0)->id!=(TP[0])->i1)EP[1]=(TP[0])->e0;
        else if(((TP[0])->e1)->io!=(TP[0])->i1 && ((TP[0])->e1)->id!=(TP[0])->i1)EP[1]=(TP[0])->e1;
        else if(((TP[0])->e2)->io!=(TP[0])->i1 && ((TP[0])->e2)->id!=(TP[0])->i1)EP[1]=(TP[0])->e2;
             if(((TP[0])->e0)->io!=(TP[0])->i2 && ((TP[0])->e0)->id!=(TP[0])->i2)EP[2]=(TP[0])->e0;
        else if(((TP[0])->e1)->io!=(TP[0])->i2 && ((TP[0])->e1)->id!=(TP[0])->i2)EP[2]=(TP[0])->e1;
        else if(((TP[0])->e2)->io!=(TP[0])->i2 && ((TP[0])->e2)->id!=(TP[0])->i2)EP[2]=(TP[0])->e2;
  
        if(EP[0]->l==NULL || ((EP[0])->l)->inverse==i)TP[0]=(EP[0])->r;
         else                                                   TP[0]=(EP[0])->l;
        if(EP[1]->l==NULL || ((EP[1])->l)->inverse==i)TP[1]=(EP[1])->r;
         else                                                   TP[1]=(EP[1])->l;
        if(EP[2]->l==NULL || ((EP[2])->l)->inverse==i)TP[2]=(EP[2])->r;
         else                                                   TP[2]=(EP[2])->l;
   
        fprintf(fout," %5d %5d %5d",EP[0]->inverse,EP[1]->inverse,EP[2]->inverse);
   
        if(TP[0]!=NULL) fprintf(fout," %5d",TP[0]->inverse);
         else           fprintf(fout," -1");
        if(TP[1]!=NULL) fprintf(fout," %5d",TP[1]->inverse);
         else           fprintf(fout," -1");
        if(TP[2]!=NULL) fprintf(fout," %5d",TP[2]->inverse);
          else           fprintf(fout," -1");
/* Circumcenter. */
/*   fprintf(fout,"%lf %lf",);*/

        TP[0]=C->triangleList[i];
        for(j=0;j<(C->vertexList[0])->nC;j++)
         {
          c=0.;
          V=C->vertexList[(TP[0])->i0];
          c+=V->center[j];
          V=C->vertexList[(TP[0])->i1];
          c+=V->center[j];
          V=C->vertexList[(TP[0])->i2];
          c+=V->center[j];
          c=c/3.;
          fprintf(fout," %lf",c);
         }
       }
      fprintf(fout," 0\n");
     }
   }
  fclose(fout);
  return;
 }

void MFFree2dCellComplex(MF2dCellComplex C,MFErrorHandler e)
 {
  int i;

  for(i=0;i<C->nVertices;i++)MFFreeMFVertex((C->vertexList)[i],e);
  free(C->vertexList);
  for(i=0;i<C->nEdges;i++)MFFreeMFEdge((C->edgeList)[i],e);
  free(C->edgeList);
  for(i=0;i<C->nTriangles;i++)MFFreeMFTriangle((C->triangleList)[i],e);
  free(C->triangleList);
  free(C);

  return;
 }

void MFOutput2dCellComplexAsDX(MF2dCellComplex C, char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"MFOutput2dCellComplexAsDX"};
  int i,j,k;
  char fullname[1000];
  FILE *fout;

  strcpy(fullname,name);
  strcat(fullname,"Dual.dx");
  fout=fopen(fullname,"w");

  printf("%s\n",RoutineName);fflush(stdout);

  fprintf(fout,"object 1 class array type float rank 1 shape %d items %d data follows\n",((C->vertexList)[0])->nC,C->nVertices);

  for(i=0;i<C->nVertices;i++)
   {
    for(j=0;j<((C->vertexList)[i])->nC;j++)
     {
      fprintf(fout,"       %lf",((C->vertexList)[i])->center[j]);
     }
    fprintf(fout,"\n");
   }
  fprintf(fout,"\n");

  fprintf(fout,"object 2 class array type int rank 1 shape 3 items %d data follows\n",C->nTriangles);

  for(i=0;i<C->nTriangles;i++)
   {
    fprintf(fout,"       %d       %d       %d\n",((C->triangleList)[i])->i0,((C->triangleList)[i])->i1,((C->triangleList)[i])->i2);
   }
  fprintf(fout,"attribute \"element type\" string \"triangles\"\n");
  fprintf(fout,"attribute \"ref\" string \"positions\"\n");
  fprintf(fout,"\n");

  fprintf(fout,"object \"mesh\" class field\n");
  fprintf(fout,"component \"positions\" value 1\n");
  fprintf(fout,"component \"connections\" value 2\n");
  fprintf(fout,"end\n");

  fclose(fout);
  return;
 }
