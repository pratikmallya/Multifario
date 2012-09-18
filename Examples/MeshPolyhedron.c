/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 */

static char *id="@(#) $Id: ComputeSphere.c,v 1.2 2007/06/08 21:01:20 mhender Exp $";
/*
    %W%
    %D% %T%

    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <MFAtlas.h>
#include <MFMultifariosMethod.h>
#include <MFEnumDualPolytope.h>
#include <MFMesh.h>

int main(int argc, char *argv[])
 {
  MFErrorHandler e;
  MFEnumDualPolytope P;
  int i;

  int nv      =  8;
  double *v[8];
  double v0[3]=  {0.,0.,0.};
  double v1[3]=  {1.,0.,0.};
  double v2[3]=  {1.,1.,0.};
  double v3[3]=  {0.,1.,0.};
  double v4[3]=  {0.,0.,1.};
  double v5[3]=  {1.,0.,1.};
  double v6[3]=  {1.,1.,1.};
  double v7[3]=  {0.,1.,1.};

  int ne        =  12;
  int *ev[12];
  int ev00[2]= {0,1};
  int ev01[2]= {1,2};
  int ev02[2]= {2,3};
  int ev03[2]= {3,0};
  int ev04[2]= {4,5};
  int ev05[2]= {5,6};
  int ev06[2]= {6,7};
  int ev07[2]= {7,4};
  int ev08[2]= {0,4};
  int ev09[2]= {1,5};
  int ev10[2]= {2,6};
  int ev11[2]= {3,7};

  int *evs[12];
  int evs00[2]={1,-1};
  int evs01[2]={1,-1};
  int evs02[2]={1,-1};
  int evs03[2]={1,-1};
  int evs04[2]={1,-1};
  int evs05[2]={1,-1};
  int evs06[2]={1,-1};
  int evs07[2]={1,-1};
  int evs08[2]={1,-1};
  int evs09[2]={1,-1};
  int evs10[2]={1,-1};
  int evs11[2]={1,-1};

  int nf=6;
  int nfe[6]   =  {4,4,4,4,4,4};

  int *fe[6];
  int *fes[6];

  int fe0 [4] = { 3, 2, 1, 0};
  int fes0[4] = {-1,-1,-1,-1};

  int fe1 [4] = { 4, 5, 6, 7};
  int fes1[4] = { 1, 1, 1, 1};

  int fe2 [4] = { 0, 9, 4, 8};
  int fes2[4] = { 1, 1,-1,-1};

  int fe3 [4] = {10, 2,11, 6};
  int fes3[4] = {-1, 1, 1,-1};

  int fe4 [4] = { 1,10, 5, 9};
  int fes4[4] = { 1, 1,-1,-1};

  int fe5 [4] = {11, 3, 8, 7};
  int fes5[4] = {-1, 1, 1,-1};

  e=MFCreateErrorHandler();

  v[0]=v0;
  v[1]=v1;
  v[2]=v2;
  v[3]=v3;
  v[4]=v4;
  v[5]=v5;
  v[6]=v6;
  v[7]=v7;

  ev[ 0]=ev00;
  ev[ 1]=ev01;
  ev[ 2]=ev02;
  ev[ 3]=ev03;
  ev[ 4]=ev04;
  ev[ 5]=ev05;
  ev[ 6]=ev06;
  ev[ 7]=ev07;
  ev[ 8]=ev08;
  ev[ 9]=ev09;
  ev[10]=ev10;
  ev[11]=ev11;

  evs[ 0]=evs00;
  evs[ 1]=evs01;
  evs[ 2]=evs02;
  evs[ 3]=evs03;
  evs[ 4]=evs04;
  evs[ 5]=evs05;
  evs[ 6]=evs06;
  evs[ 7]=evs07;
  evs[ 8]=evs08;
  evs[ 9]=evs09;
  evs[10]=evs10;
  evs[11]=evs11;

  fe[0]=fe0;
  fe[1]=fe1;
  fe[2]=fe2;
  fe[3]=fe3;
  fe[4]=fe4;
  fe[5]=fe5;

  fes[0]=fes0;
  fes[1]=fes1;
  fes[2]=fes2;
  fes[3]=fes3;
  fes[4]=fes4;
  fes[5]=fes5;

  P=MFGenerateMeshInPolyhedron(                     nv,         v,    ne,      ev,       evs,     nf,      nfe,       fe,      fes,       .1,               e);

/*
  MFFreeEnumDualPolytope(P,e);
 */
  MFFreeErrorHandler(e);

  return 0;
 }
