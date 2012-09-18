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

static char *id="@(#) $Id: MFListOfCharts.c,v 1.3 2007/02/13 01:22:33 mhender Exp $";

static char MFListOfChartsErrorMsg[256]="";

#include <MFListOfCharts.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFListOfChartsSt {
                       int n;
                       int *charts;
                      };

MFListOfCharts MFCreateListOfCharts(int n,int *charts, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateListOfCharts"};
  MFListOfCharts thisList;
  int i,j,k;

  thisList=(MFListOfCharts)malloc(sizeof(struct MFListOfChartsSt));

#ifndef MFNOSAFETYNET
  if(thisList==NULL)
   {
    sprintf(MFListOfChartsErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFListOfChartsSt));
    MFSetError(e,12,RoutineName,MFListOfChartsErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  thisList->n=n;
  thisList->charts=charts;

  for(i=0;i<n-1;i++)
   for(j=i+1;j<n;j++)
    {
     if(thisList->charts[i]>thisList->charts[j])
      {
       k=thisList->charts[j];
       thisList->charts[j]=thisList->charts[i];
       thisList->charts[i]=k;
      }
    }

  return thisList;
 }

int MFNumberOfIntersectingCharts(MFListOfCharts L, MFErrorHandler r)
 {
  static char RoutineName[]={"MFNumberOfIntersectingCharts"};

  return L->n;
 }

int MFIntersectingChart(MFListOfCharts L,int chart, MFErrorHandler r)
 {
  static char RoutineName[]={"MFIntersectingChart"};

  return L->charts[chart];
 }

void MFFreeListOfIntersectingCharts(MFListOfCharts L, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeListOfIntersectingCharts"};

  if(L->charts!=NULL)free(L->charts);
  free(L);

  return;
 }

#ifdef __cplusplus
}
#endif
