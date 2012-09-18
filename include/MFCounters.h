/* 
    @(#)MFCounters.h	1.2
    02/04/19 14:40:16

    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

int MFNProjectionsToM=0;
int MFNGetTangentSpaces=0;
double MFMinEps=0.;
double MFMaxEps=0.;

/*  KVectors */

   int MFNNVectorsCreated=0;
   int MFNNVectorsFreed=0;
   int MFNMaxNVectorsDefined=0;
   int MFNVectorSetC=0;
   int MFNVectorGetC=0;
   int MFNVectorLength=0;
   int MFNVectorSums=0;
   int MFNVectorDiffs=0;

/*  KVectors */

   int MFNKVectorsCreated=0;
   int MFNKVectorsFreed=0;
   int MFKMaxVectorsDefined=0;
   int MFKVectorSetC=0;
   int MFKVectorGetC=0;
   int MFKVectorLength=0;
   int MFKVectorDots=0;
   int MFKVectorNorms=0;
   int MFKVectorScale=0;
   int MFKVectorScaleMul=0;
   int MFKVectorDiffs=0;
   int MFKVectorAdds=0;

/* Tangent Spaces */

   int MFNNMatrixsCreated=0;
   int MFNNKMatrixsFreed=0;
   int MFNKMatrixMaxDefined=0;
   int MFNKMatrixGetColumns=0;
   int MFNKMatrixGetRows=0;
   int MFNKMatrixaMult=0;
   int MFNKMatrixaMultT=0;

int MFNNonIntCharts=0;
int MFNBChartChk=0;
double MFMaxR=0.;
double MFMinR=0.;
int MFNResizes=0;
