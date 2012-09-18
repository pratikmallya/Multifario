/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   February 22, 1999
 *              October 6, 2004   Added ProjectFromCenter
 */

static char *id="@(#) $Id: MFImplicitMF.c,v 1.7 2010/05/21 12:22:04 mhender Exp $";

#include <MFImplicitMF.h>
#include <MFKVector.h>
#include <MFNVector.h>
#include <MFNKMatrix.h>
#include <MFNSpace.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <MFNVector.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

static char MFImplicitMFErrorHandlerMsg[256]="";

void MFIMFSetK(MFImplicitMF,int,MFErrorHandler);

struct MFImplicitMFSt {
                       MFNSpace space;
                       int n;
                       int k;

                       void *data;
                       void (*writedata)(FILE*,void*,MFErrorHandler);
                       void (*freedata)(void*,MFErrorHandler);
                       int (*projectFromCenter)(int,int,MFNVector,MFNKMatrix,MFKVector,MFNVector,void*,int*,MFErrorHandler);
                       int (*project)(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
                       int (*tangent)(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
                       int (*tangentWithGuess)(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler);
                       void (*evaluate)(int,MFNVector,MFNVector,void*,MFErrorHandler);
                       void (*applyJacobian)(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler);
                       void (*applySecDer)(int,int,MFNVector,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
                       double (*scale)(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
                       int (*stop)(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,void*,MFErrorHandler);
                       int (*singular)(int,int,MFNVector,MFNKMatrix,MFNVector,void*,MFErrorHandler);
                       void (*stability)(MFImplicitMF,MFNVector,MFNKMatrix,void*,MFErrorHandler);
                       MFNVector (*vectorFactory)(MFImplicitMF,MFErrorHandler);
                       MFNKMatrix (*matrixFactory)(MFImplicitMF,MFErrorHandler);
                       double R;
                       double RMin;

                       char *id;
                       int nRefs;
                       int (*saver)(MFNVector,double*,void*,MFErrorHandler);
                       int (*drawer)(MFNVector,double*,void*,MFErrorHandler);
                       int (*bbprojector)(MFNVector,double*,void*,MFErrorHandler);
                      };

/*! \fn void MFRefImplicitMF(MFImplicitMF M,MFErrorHandler e);
 *  \brief Adds a references to an MFImplicitMF.
 *
 *  \param M The an ImplicitMF being referenced.
 *  \param e A place to handle errors.
 *  \sa ReferenceCounting MFFreeImplicitMF
 */
void MFRefImplicitMF(MFImplicitMF M,MFErrorHandler e)
 {
  static char RoutineName[]={"MFRefImplicitMF"};

  if(M==NULL)
   {
#ifdef MFNOCONFIDENCE
    sprintf(MFImplicitMFErrorHandlerMsg,"Pointer to ImplicitMF (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
#endif
    return;
   }

  M->nRefs++;
  return;
 }

/*! \fn char *MFImplicitMFId(MFImplicitMF M);
 *  \brief Returns a character string which assigns a type to an ImplicitMF. Do not change or free the string!
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns The type of the ImplicitMF.
 */
char *MFImplicitMFId(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFImplicitMFId"};

  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Pointer to ImplicitMF (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return NULL;
   }

  return M->id;
 }

/*! \fn void MFFreeImplicitMF(MFImplicitMF M);
 *  \brief Frees a references to the ImplicitMF, and deletes the ImplicitMF if there are no references left.
 *
 *  \param M The ImplicitMF being unreferenced.
 *  \param e A place to handle errors.
 *  \sa ReferenceCounting MFRefImplicitMF
 */
void MFFreeImplicitMF(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeImplicitMF"};

  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Pointer to ImplicitMF (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }

  M->nRefs--;

  if(M->nRefs<1)
   {
    if(M->space!=NULL)MFFreeNSpace(M->space,e);
    if(M->freedata!=NULL && M->data!=NULL)(M->freedata)(M->data,e);
    if(M->id!=NULL)free(M->id);
    free(M);
   }
  return;
 }

/*! \fn int MFIMF_N(MFImplicitMF M,MFErrorHandler e);
 *  \brief Gets the embedding space dimension of an ImplicitMF.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns The dimension of the embedding space (n).
 */
int MFIMF_N(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMF_N"};

  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Pointer to ImplicitMF (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return -1;
   }

  return M->n;
 }

/*! \fn int MFIMF_K(MFImplicitMF M,MFErrorHandler e);
 *  \brief Gets the dimension of an ImplicitMF.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns The dimension of the manifold (k).
 */
int MFIMF_K(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMF_K"};

  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Pointer to ImplicitMF (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return -1;
   }

  return M->k;
 }

/*! \fn MFNSpace MFIMFNSpace(MFImplicitMF M,MFErrorHandler e);
 *  \brief Gets the space in which the manifold is embedded.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns The embedding space.
 */
MFNSpace MFIMFNSpace(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFNSpace"};

  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Pointer to ImplicitMF (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return NULL;
   }

  return M->space;
 }

/*! \fn int MFIMFProjectFromCenter(MFImplicitMF M,MFNVector u0,MFNKMatrix Phi0,MFKVector s,MFNVector u, MFErrorHandler e);
 *  \brief This is the generalizatio of the prediction and correction operations in pseudo-arclength continuation.
 *         The prediction is u0+Phi0.s, and the correction is from thisIDM point to M orthogonal to Phi0.
 *
 *  \param M The implicitly defined manifold.
 *  \param u0 The point to be projected.
 *  \param Phi0 The tangent space.
 *  \param s The point in the tangent space that is to be projected.
 *  \param u The projected point, which lies on M and Phi0^T(u-u0)-0.
 *  \param e A place to handle errors.
 *  \returns TRUE if the projection was sucessful
 */
int MFIMFProjectFromCenter(MFImplicitMF M,MFNVector u0,MFNKMatrix Phi,MFKVector s, MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFProjectFromCenter"};
  int rc;
  int index;
  if(M->projectFromCenter==NULL)
   {
    MFNVector up;
    MFNVector v;
    up=MFCloneNVector(u0,e);
    v=MFCloneNVector(u,e);
    MFMVMul(MFIMFNSpace(M,e),Phi,s,v,e);
    MFNSpaceAdd(MFIMFNSpace(M,e),u0,v,up,e);
    rc=M->project(M->n,M->k,up,Phi,u,M->data,&index,e);
    MFNVSetIndex(u,index,e);
    MFFreeNVector(up,e);
    MFFreeNVector(v,e);
    return rc;
   }

  rc=M->projectFromCenter(M->n,M->k,u0,Phi,s,u,M->data,&index,e);
  MFNVSetIndex(u,index,e);

  return rc;
 }

/*! \fn int MFIMFProject(MFImplicitMF M,MFNVector u0 ,MFNKMatrix Phi0,MFNVector u);
 *  \brief Projects a point u0 onto an ImplicitMF orthogonal to a tangent space Phi0
 *
 *  \param M The implicitly defined manifold.
 *  \param u0 The point to be projected.
 *  \param Phi0 The tangent space.
 *  \param u The projected point, which lies on M and Phi0^T(u-u0)-0.
 *  \param e A place to handle errors.
 *  \returns TRUE if the projection was sucessful
 */
int MFIMFProject(MFImplicitMF M,MFNVector u0,MFNKMatrix Phi,MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFProject"};
  int rc;
  int index;

  if(M->project==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"No Project Routine has been supplied");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
   }

  rc=M->project(M->n,M->k,u0,Phi,u,M->data,&index,e);
  MFNVSetIndex(u,index,e);

  return rc;
 }

/*! \fn MFNKMatrix MFIMFTangentSpace(MFImplicitMF M,MFNVector u);
 *  \brief Computes and returns an orthonormal basis for the tangent space of M at the point u on M.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point at which the tangent space is needed.
 *  \param e A place to handle errors.
 *  \returns An orthonormal basis for the tangent space of M at u on M.
 */
MFNKMatrix MFIMFTangentSpace(MFImplicitMF M,MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFTangentSpace"};
  MFNVector *col;
  MFNKMatrix Phi;
  int i;
  int rc;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  col=(MFNVector*)malloc(M->k*sizeof(MFNVector));

#ifndef MFNOSAFETYNET
  if(col==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",M->k*sizeof(MFNVector));
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  if(M->tangent==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"No Tangent Routine has been supplied");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
   }

  for(i=0;i<M->k;i++)col[i]=MFCloneNVector(u,e);
  Phi=MFCreateNKMatrix(M->k,col,e);
  for(i=0;i<M->k;i++)MFFreeNVector(col[i],e);
  free(col);

  rc=M->tangent(M->n,M->k,u,Phi,M->data,e);
  if(!rc)
   {
    MFFreeNKMatrix(Phi,e);
    return NULL;
   }

  if(verbose)
   {
    int j,jj;
    MFNVector phi0;
    MFNVector phi1;
    {printf("       check basis for tangent space\n");fflush(stdout);}
    for(jj=0;jj<M->k;jj++)
     for(j=0;j<M->k;j++)
      {
       phi0=MFMColumn(Phi,jj,e);
       phi1=MFMColumn(Phi,j,e);
  printf("Line %d, File %s\n",__LINE__,__FILE__);
       printf(" <phi_%d,phi_%d>=%le\n",jj,j,MFNSpaceInner(M->space,phi0,phi1,e));
       MFFreeNVector(phi0,e);
       MFFreeNVector(phi1,e);
      }
    fflush(stdout);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return Phi;
 }

/*! \fn MFNKMatrix MFIMFTangentSpaceWithGuess(MFImplicitMF M,MFNVector u,MFNKMatrix guess,MFErrorHandler e);
 *  \brief Computes and returns an orthonormal basis for the tangent space of M at the point u on M.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point at which the tangent space is needed.
 *  \param guess An approximate tangent space. The new basis will align roughly with the basis in guess.
 *  \param e A place to handle errors.
 *  \returns An orthonormal basis for the tangent space of M at u on M.
 */
MFNKMatrix MFIMFTangentSpaceWithGuess(MFImplicitMF M,MFNVector u,MFNKMatrix Phi0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFTangentSpaceWithGuess"};
  MFNKMatrix Phi;
  int verbose=0;
  int rc;

  if(M->tangentWithGuess==NULL)
   {
    return MFIMFTangentSpace(M,u,e);
   }

  Phi=MFCloneNKMatrix(Phi0,e);
  rc=M->tangentWithGuess(M->n,M->k,u,Phi0,Phi,M->data,e);
  if(!rc)
   {
    MFFreeNKMatrix(Phi,e);
    return NULL;
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    int j,jj;
    MFNVector phi0;
    MFNVector phi1;
    printf("%s\n",RoutineName);
    {printf("       check basis for tangent space\n");fflush(stdout);}
    for(jj=0;jj<M->k;jj++)
     for(j=0;j<M->k;j++)
      {
       phi0=MFMColumn(Phi,jj,e);
       phi1=MFMColumn(Phi,j,e);
       printf(" <phi_%d,phi_%d>=%le\n",jj,j,MFNSpaceInner(M->space,phi0,phi1,e));
       MFFreeNVector(phi0,e);
       MFFreeNVector(phi1,e);
      }
    fflush(stdout);
   }
#endif

  return Phi;
 }

/*! \fn double MFIMFScale(MFImplicitMF M,MFNVector u,MFNKMatrix Tan);
 *  \brief Estimates the radius of a spherical ball in the tangent space of M at u on M.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point at which the tangent space is needed.
 *  \param Tan The tangent space of M at u.
 *  \param e A place to handle errors.
 *  \returns The radius.
 */
double MFIMFScale(MFImplicitMF M,MFNVector u, MFNKMatrix Phi, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFScale"};
  double result;

  if(M->scale==NULL&&M->R>0.)
    result=M->R;
   else if(M->scale==NULL&&M->R<0.)
    result=1.;
   else{
    result=M->scale(M->n,M->k,u,Phi,M->data,e);
    if(M->R>0.&&result>M->R)result=M->R;
   }
  return result;
 }

/*! \fn void MFWriteImplicitMF(FILE* fid,MFImplicitMF M);
 *  \brief Writes a ImplicitMF to a file.
 *
 *  \param fid The file to write to.
 *  \param M The ImplicitMF being written.
 *  \param e A place to handle errors.
 */
void MFWriteImplicitMF(FILE *fid,MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteImplicitMF"};
  int i;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  fprintf(fid,"%s\n","ImplicitMF");
  fprintf(fid,"%d\n",strlen(M->id));
  fprintf(fid,"%s\n",M->id);fflush(fid);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf(" tag is -->%s<--\n","ImplicitMF");
    printf(" strlen(id) is %d\n",strlen(M->id));
    printf(" id is -->%s<--\n",M->id);fflush(fid);
   }
#endif

  M->writedata(fid,M->data,e);
  fprintf(fid,"%lf %d\n",M->R,M->nRefs);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf(" R is %lf\n",M->R);fflush(fid);
    printf(" nRefs is %d\n",M->nRefs);fflush(fid);
   }
#endif

  return;
 }

#ifdef TPBVPIMF
MFImplicitMF MFReadTPBVP(FILE*,MFErrorHandler);
#endif
#ifdef PENDULAIMF
MFImplicitMF MFReadPendula(FILE*,MFErrorHandler);
#endif
#ifdef POLYGONIN3SPACEIMF
MFImplicitMF MFReadPolygonIn3Space(FILE*,MFErrorHandler);
#endif
#ifdef CSTRIMF
MFImplicitMF MFReadCSTR(FILE*,MFErrorHandler);
#endif
#ifdef FLATIMF
MFImplicitMF MFReadFlat(FILE*,MFErrorHandler);
#endif
#ifdef NSPACEIMF
MFImplicitMF MFReadNSpaceMF(FILE*,MFErrorHandler);
#endif
#ifdef PLANEIMF
MFImplicitMF MFReadPlane(FILE*,MFErrorHandler);
#endif
#ifdef SPHEREIMF
MFImplicitMF MFReadSphere(FILE*,MFErrorHandler);
#endif
#ifdef CIRCLEIMF
MFImplicitMF MFReadCircle(FILE*,MFErrorHandler);
#endif
#ifdef SWALLOWIMF
MFImplicitMF MFReadSwallow(FILE*,MFErrorHandler);
#endif
#ifdef EXPRESSIONIMF
MFImplicitMF MFReadExpression(FILE*,MFErrorHandler);
#endif
#ifdef INVARIANTMF
MFImplicitMF MFReadInvariantMF(FILE*,MFErrorHandler);
#endif
#ifdef PEITGENIMF
MFImplicitMF MFReadPeitgenSystem(FILE*,MFErrorHandler);
#endif
#ifdef TORUSIMF
MFImplicitMF MFReadTorus(FILE*,MFErrorHandler);
#endif

/*! \fn MFImplicitMF MFReadImplicitMF(FILE* fid, MFErrorHandler e);
 *  \brief Reads a ImplicitMF from a file.
 *
 *  \param fid The file to read from.
 *  \param e A place to handle errors.
 */
MFImplicitMF MFReadImplicitMF(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadImplicitMF"};
  int i;
  MFImplicitMF M;
  char tag[100]="";
  char *id;
  int n=0;
  int nRefs=0;
  double R=0.;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  fscanf(fid,"%s\n",tag);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" tag is -->%s<--\n",tag);fflush(stdout);}
#endif

  if(strcmp(tag,"ImplicitMF"))
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Next Object is not a ImplicitMF! (%s)\n",RoutineName,tag);
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return M;
   }
  M=NULL;
  fscanf(fid,"%d\n",&n);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" strlen(id) is %d\n",n);fflush(stdout);}
#endif

  id=(char*)malloc((n+1)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(id==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(n+1)*sizeof(char));
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  fscanf(fid,"%s\n",id);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" id is -->%s<--, length %d\n",id,strlen(id));fflush(stdout);}
#endif

  if(0) {}
#ifdef TORUSIMF
   else if(!strcmp(id,"Torus"))
    M=MFReadTorus(fid);
#endif
#ifdef PENDULAIMF
   else if(!strcmp(id,"Pendula"))
    M=MFReadPendula(fid);
#endif
#ifdef POLYGONIN3SPACEIMF
   else if(!strcmp(id,"PolygonIn3Space"))
    M=MFReadPolygonIn3Space(fid);
#endif
#ifdef CSTRIMF
   else if(!strcmp(id,"CSTR"))
    M=MFReadCSTR(fid);
#endif
#ifdef FLATIMF
   else if(!strcmp(id,"Flat"))
    M=MFReadFlat(fid);
#endif
#ifdef NSPACEIMF
   else if(!strcmp(id,"NSpaceMF"))
    M=MFReadNSpaceMF(fid);
#endif
#ifdef PLANEIMF
   else if(!strcmp(id,"Plane"))
    M=MFReadPlane(fid);
#endif
#ifdef SPHEREIMF
   else if(!strcmp(id,"Sphere"))
    M=MFReadSphere(fid);
#endif
#ifdef CIRCLEIMF
   else if(!strcmp(id,"Circle"))
    M=MFReadCircle(fid);
#endif
#ifdef SWALLOWIMF
   else if(!strcmp(id,"Swallow"))
    M=MFReadSwallow(fid);
#endif
#ifdef EXPRESSIONIMF
   else if(!strcmp(id,"Expression"))
    M=MFReadExpression(fid);
#endif
#ifdef INVARIANTMF
   else if(!strcmp(id,"InvariantMF"))
    M=MFReadInvariantMF(fid);
#endif
#ifdef PEITGENMF
   else if(!strcmp(id,"PeitgenSystem"))
    M=MFReadPeitgenSystem(fid);
#endif
#ifdef TPBVPIMF
   else if(!strcmp(id,"TPBVP"))
    M=MFReadTPBVP(fid);
#endif
   else 
    {
     free(id);
     return M;
    }
  free(id);

  if(M!=NULL)
   {
    M->space=NULL;
    fscanf(fid,"%lf %d\n",&(M->R),&(M->nRefs));

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf(" R is %lf\n",M->R);fflush(stdout);
      printf(" nRefs is %d\n",M->nRefs);fflush(stdout);
     }
#endif

   }else
    fscanf(fid,"%lf %d\n",&R,&nRefs);
  
  return M;
 }

void MFIMFSetSpace(MFImplicitMF thisIDM,MFNSpace space, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetSpace"};

  if(thisIDM==NULL)return;
  if(space!=NULL)MFRefNSpace(space,e);
  if(thisIDM->space!=NULL)MFFreeNSpace(thisIDM->space,e);
  thisIDM->space=space; 
  return;
 }

void MFIMFEvaluate(MFImplicitMF thisIDM,MFNVector u,MFNVector f, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFEvaluate"};

  if(thisIDM->evaluate==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"No Routine to evaluate the function");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }

  thisIDM->evaluate(thisIDM->n,u,f,thisIDM->data,e);

  return;
 }

void MFIMFApplyJacobian(MFImplicitMF thisIDM,MFNVector u,MFNKMatrix Phi0,MFNKMatrix Phi1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFApplyJacobian"};
  int n ,k ;
  int n0,k0;
  int n1,k1;

  n=thisIDM->n;
  k=thisIDM->k;
  n0=MFNKMatrixN(Phi0,e);
  k0=MFNKMatrixK(Phi0,e);
  n1=MFNKMatrixN(Phi1,e);
  k1=MFNKMatrixK(Phi1,e);

  if(n0!=n||n1!=n||k1!=k0)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Incompatible matrices for apply jacobian (%d,%d) = (%d,%d).(%d,%d)",n1,k1,n,k,n0,k0);
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }

  if(thisIDM->applyJacobian==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"No Routine to apply the Jacobian has been supplied");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }

  thisIDM->applyJacobian(thisIDM->n,k0,u,Phi0,Phi1,thisIDM->data,e);

  return;
 }

void MFIMFSetK(MFImplicitMF thisIDM,int k, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetK"};

  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Pointer to ImplicitMF (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }

  thisIDM->k=k;

  return;
 }

int MFIMFStop(MFImplicitMF thisIDM, MFNVector u0, MFNKMatrix Phi0, MFNVector u1, MFNKMatrix Phi1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFStop"};

  if(thisIDM->stop!=(NULL))
    return thisIDM->stop(thisIDM,u0,Phi0,u1,Phi1,thisIDM->data,e);
   else
    return 0;
 }

/*! \fn int MFIMFProjectToSave(MFImplicitMF M,MFNVector u,double *x, MFerror e);
 *  \brief Projects a point so that it can be saved to file (used for the centerfile).
 *
 *  The convention is that if x is NULL, the routine returns the dimension of x, so that the user can allocate storage.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point that is to be projected, or NULL if the dimension is needed.
 *  \param x The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns The dimension of the projected point if x is NULL.
 */
int MFIMFProjectToSave(MFImplicitMF M, MFNVector u, double *x, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFProjectToSave"};

  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Pointer to ImplicitMF (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return -1;
   }

  if(M->saver==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Project for Save is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return -1;
   }

  return (M->saver)(u,x,M->data,e);
 }

/*! \fn int MFIMFProjectToDraw(MFImplicitMF M,MFNVector u,double *x, MFerror e);
 *  \brief Projects a point so that it can be drawn (used for the plotfile).
 *
 *  The convention is that if x is NULL, the routine returns the dimension of x, so that the user can allocate storage.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point that is to be projected, or NULL if the dimension is needed.
 *  \param x The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns The dimension of the projected point if x is NULL.
 */
int MFIMFProjectToDraw(MFImplicitMF M, MFNVector u, double *x, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFProjectToDraw"};

  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Pointer to ImplicitMF (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return -1;
   }

  if(M->drawer==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Project for Draw is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return -1;
   }

  return (M->drawer)(u,x,M->data,e);
 }

/*! \fn int MFIMFProjectToBB(MFImplicitMF M,MFNVector u,double *x, MFerror e);
 *  \brief Projects a point so that it can be placed in the hierarchical bounding boxes that are used to speed up
 *         the location of neighboring charts.
 *
 *  The convention is that if x is NULL, the routine returns the dimension of x, so that the user can allocate storage.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point that is to be projected, or NULL if the dimension is needed.
 *  \param x The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns The dimension of the projected point if x is NULL.
 */
int MFIMFProjectToBB(MFImplicitMF M, MFNVector u, double *x, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFProjectToBB"};

  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Pointer to ImplicitMF (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return -1;
   }

  if(M->bbprojector==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Project for BB is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return -1;
   }

  return (M->bbprojector)(u,x,M->data,e);
 }

MFImplicitMF MFIMFCreateBaseClass(int n, int k, char *id, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateBaseClass"};
  MFImplicitMF thisIDM;

  if(id==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Id (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(k<0)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"k (%d) (argument 3) is invalid, must be positive",k);
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return NULL;
   }

  thisIDM=(MFImplicitMF)malloc(sizeof(struct MFImplicitMFSt));

#ifndef MFNOSAFETYNET
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFImplicitMFSt));
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  thisIDM->id=(char*)malloc((strlen(id)+1)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(thisIDM->id==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(strlen(id)+1)*sizeof(char));
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisIDM);
    return NULL;
   }
#endif
  strcpy(thisIDM->id,id);

  thisIDM->space=NULL;
  thisIDM->n=n;
  thisIDM->k=k;
  thisIDM->nRefs=1;
  thisIDM->data=NULL;
  thisIDM->writedata=NULL;
  thisIDM->freedata=NULL;
  thisIDM->project=NULL;
  thisIDM->projectFromCenter=NULL;
  thisIDM->tangent=NULL;
  thisIDM->tangentWithGuess=NULL;
  thisIDM->evaluate=NULL;
  thisIDM->applyJacobian=NULL;
  thisIDM->applySecDer=NULL;
  thisIDM->vectorFactory=NULL;
  thisIDM->matrixFactory=NULL;
  thisIDM->scale=NULL;
  thisIDM->stop=NULL;
  thisIDM->stability=NULL;
  thisIDM->singular=NULL;
  thisIDM->R=-1;
  thisIDM->RMin=-1;
  thisIDM->saver=NULL;
  thisIDM->drawer=NULL;
  thisIDM->bbprojector=NULL;
  thisIDM->data=NULL;

  return thisIDM;
 }

void MFIMFSetData(MFImplicitMF thisIDM,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetData"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->data=data;
  return;
 }

/*! \fn void *MFIMFGetData(MFImplicitMF M, MFerror e);
 *  \brief Returns a pointer to the internal data used by the manifold. Be very careful using thisIDM. You need to know what you're
 *         doing.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns A pointer to the internal data.
 */
void *MFIMFGetData(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFGetData"};

#ifdef MFNOCONFIDENCE
  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return M->data;
 }

void MFIMFSetWriteData(MFImplicitMF thisIDM,void (*writedata)(FILE*,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetWriteData"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->writedata=writedata;
  return;
 }

void MFIMFSetFreeData(MFImplicitMF thisIDM,void (*freedata)(void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetFreeData"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->freedata=freedata;
  return;
 }

void MFIMFSetVectorFactory(MFImplicitMF thisIDM,MFNVector (*factory)(MFImplicitMF,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetVectorFactory"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->vectorFactory=factory;
  return;
 }

void MFIMFSetMatrixFactory(MFImplicitMF thisIDM,MFNKMatrix (*factory)(MFImplicitMF,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetMatrixFactory"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->matrixFactory=factory;
  return;
 }

void MFIMFSetProjectFromCenter(MFImplicitMF thisIDM,int (*projectFromCenter)(int,int,MFNVector,MFNKMatrix,MFKVector,MFNVector,void*,int*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetProjectFromCenter"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->projectFromCenter=projectFromCenter;
  return;
 }

void MFIMFSetProject(MFImplicitMF thisIDM,int (*project)(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetProject"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->project=project;
  return;
 }

void MFIMFSetTangent(MFImplicitMF thisIDM,int (*tangent)(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetTangent"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->tangent=tangent;
  return;
 }

void MFIMFSetTangentWithGuess(MFImplicitMF thisIDM,int (*tangentWithGuess)(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetTangentWithGuess"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->tangentWithGuess=tangentWithGuess;
  return;
 }

void MFIMFSetEvaluate(MFImplicitMF thisIDM,void (*evaluate)(int,MFNVector,MFNVector,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetEvaluate"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->evaluate=evaluate;
  return;
 }

void MFIMFSetApplyJacobian(MFImplicitMF thisIDM,void (*applyJacobian)(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetApplyJacobian"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->applyJacobian=applyJacobian;
  return;
 }

void MFIMFSetScale(MFImplicitMF thisIDM,double (*scale)(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetScale"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->scale=scale;
  return;
 }

void MFIMFSetStop(MFImplicitMF thisIDM,int (*stop)(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetStop"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->stop=stop;
  return;
 }

void MFIMFSetSingular(MFImplicitMF thisIDM,int (*singular)(int,int,MFNVector,MFNKMatrix,MFNVector,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetSingular"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->singular=singular;
  return;
 }

/*! \fn double MFIMFGetR(MFImplicitMF M, MFerror e);
 *  \brief Returns the radius currently associated with a manifold. It provides a guess if no other information about the
 *          scale is available.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns The radius.
 */
double MFIMFGetR(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFGetR"};

#ifdef MFNOCONFIDENCE
  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return M->R;
 }

/*! \fn void MFIMFSetR(MFImplicitMF M,double R, MFerror e);
 *  \brief Associates a radius with a manifold. It provides a guess if no other information about the scale is available.
 *
 *  \param M The implicitly defined manifold.
 *  \param R The radius.
 *  \param e A place to handle errors.
 */
void MFIMFSetR(MFImplicitMF M,double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetR"};

#ifdef MFNOCONFIDENCE
  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  M->R=R;
  return;
 }

/*! \fn double MFIMFGetRMin(MFImplicitMF M, MFerror e);
 *  \brief Returns the current value of the radius that the user has provided as a smallest radius for the manifold.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns The radius.
 */
double MFIMFGetRMin(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFGetRMin"};

#ifdef MFNOCONFIDENCE
  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return M->RMin;
 }

/*! \fn void MFIMFSetRMin(MFImplicitMF M,double R, MFerror e);
 *  \brief Allows the user to set a smallest radius for a manifold. A continuation method may ignore it, but it is there.
 *
 *  \param M The implicitly defined manifold.
 *  \param R The radius.
 *  \param e A place to handle errors.
 */
void MFIMFSetRMin(MFImplicitMF M,double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetRMin"};

#ifdef MFNOCONFIDENCE
  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  M->RMin=R;
  return;
 }

void MFIMFSetProjectForSave(MFImplicitMF M,int (*saver)(MFNVector,double*,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetProjectForSave"};

#ifdef MFNOCONFIDENCE
  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  M->saver=saver;
  return;
 }

void MFIMFSetProjectForDraw(MFImplicitMF thisIDM,int (*drawer)(MFNVector,double*,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetProjectForDraw"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->drawer=drawer;
  return;
 }

void MFIMFSetProjectForBB(MFImplicitMF thisIDM,int (*bbprojector)(MFNVector,double*,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetProjectForBB"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->bbprojector=bbprojector;
  return;
 }

int MFIMFSingular(MFImplicitMF thisIDM,MFNVector u,MFNKMatrix Phi,MFNVector w, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSingular"};

  if(thisIDM->singular!=(NULL))
    return thisIDM->singular(thisIDM->n,thisIDM->k,u,Phi,w,thisIDM->data,e);
   else
    return 0;
 }

void MFIMFApplySecDer(MFImplicitMF thisIDM,MFNVector u,MFNVector v0,MFNVector v1,MFNVector w, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFApplySecDer"};
  int n ,k ;
  int n0,k0;
  int n1,k1;

  n=thisIDM->n;
  k=thisIDM->k;
  n0=MFNV_NC(v0,e);
  n1=MFNV_NC(v1,e);

#ifdef MFNOCONFIDENCE
  if(n0!=n||n1!=n)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Incompatible vectors for apply scond derivative %d = %d,%d",n,n0,n1);
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }

  if(thisIDM->applySecDer==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"No Routine to apply the second derivative has been supplied");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->applySecDer(thisIDM->n,k0,u,v0,v1,w,thisIDM->data,e);

  return;
 }

void MFIMFSetApplySecDer(MFImplicitMF thisIDM,void (*applySecDer)(int,int,MFNVector,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetApplyJacobian"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->applySecDer=applySecDer;
  return;
 }


void MFIMFSetStability(MFImplicitMF thisIDM, MFNVector u0, MFNKMatrix Phi0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetStability"};

  if(thisIDM->stability!=(NULL)) thisIDM->stability(thisIDM,u0,Phi0,thisIDM->data,e);

  return;
 }

void MFIMFSetSetStability(MFImplicitMF thisIDM,void (*stability)(MFImplicitMF,MFNVector,MFNKMatrix,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSetSetStability"};

#ifdef MFNOCONFIDENCE
  if(thisIDM==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisIDM->stability=stability;
  return;
 }

/*! \fn MFNKMatrix MFIMFMatrixFactory(MFImplicitMF M, MFerror e);
 *  \brief This is a factory to create an NKMatrix that is the approriate type for a point on thisIDM manifold. Unless
 *         the manifold uses a dense array of doubles the columns of the matrix will be NVectors of the appropriate type.
 *  \param M The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns A clean and shiny new NKMatrix of the right type for the tangent space of thisIDM manifold.
 */
MFNKMatrix MFIMFMatrixFactory(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFMatrixFactory"};
  MFNKMatrix Phi;
  int i;
  int rc;
  int verbose=0;

#ifdef MFNOCONFIDENCE
  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  if(M->matrixFactory==NULL)return NULL;
   else return M->matrixFactory(M,e);
 }

/*! \fn MFNVector MFIMFVectorFactory(MFImplicitMF M, MFerror e);
 *  \brief This is a factory to create an NVector that is the approriate type for a point on thisIDM manifold. All allocation
 *         of vectors in multifario is by cloning, so the type of thisIDM vector is important. NVectors and ImplicitMF's are
 *         both base classes, so the user has no other way (besides the documentation) of knowing the type of NVector to
 *         use as an initial point.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to handle errors.
 *  \returns A clean and shiny new NVector of the right type for computations of thisIDM manifold.
 */
MFNVector MFIMFVectorFactory(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFVectorFactory"};
  MFNVector Phi;
  int i;
  int rc;
  int verbose=0;

#ifdef MFNOCONFIDENCE
  if(M==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  if(M->vectorFactory==NULL)return NULL;
   else return M->vectorFactory(M,e);
 }

MFNVector MFNVectorFactory(MFImplicitMF thisIDM, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorFactory"};

  return MFCreateNVector(MFIMF_N(thisIDM,e),e);
 }

MFNKMatrix MFNKMatrixFactory(MFImplicitMF thisIDM, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNKMatrixFactory"};
  MFNKMatrix A;
  int i,k,n;
  double *data;

  k=MFIMF_K(thisIDM,e);
  n=MFIMF_N(thisIDM,e);

  data=(double*)malloc(n*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFImplicitMFErrorHandlerMsg,"Out of memory trying to allocate %d doubles",n*k);
    MFSetError(e,4,RoutineName,MFImplicitMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<n*k;i++)data[i]=0.;
  A=MFCreateNKMatrixWithData(n,k,data,e);
  free(data);
  return A;
 }

int (*MFIMFGetProjectForSave(MFImplicitMF this,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler)
 {
  return this->saver;
 }

int (*MFIMFGetProjectForDraw(MFImplicitMF this,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler)
 {
  return this->drawer;
 }

int (*MFIMFGetProjectForBB  (MFImplicitMF this,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler)
 {
  return this->bbprojector;
 }

/*! @} */

#ifdef __cplusplus
}
#endif
