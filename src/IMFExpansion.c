/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   November 11, 1997
 *              February 2, 1999   Ported to C
 *              April 29, 2002     Made into a base class
 */

static char *id="@(#) $Id: IMFExpansion.c,v 1.7 2011/07/21 17:42:46 mhender Exp $";

static char IMFExpansionErrorMsg[256];

#include <MFAtlas.h>
#include <IMFExpansion.h>
#include <IMFFlow.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <IMF.h>

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifdef __cplusplus
 extern "C" {
#endif

struct IMFExpansionSt {int n;
                       int k;
                       double *u;
                       double *du;
                       double *ddu;
                       double *dddu;
                       int dataLn;
                       double *data;
                       int nRefs;
                      };

IMFExpansion IMFCreateExpansion(int n, int k,MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateExpansion"};
  IMFExpansion result;

  result=(struct IMFExpansionSt*)malloc(sizeof(struct IMFExpansionSt));

#ifndef MFNOSAFETYNET
  if(result==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->n=n;
  result->k=k;

  result->u=NULL;
  result->du=NULL;
  result->ddu=NULL;
  result->dddu=NULL;
  result->dataLn=0;
  result->data=NULL;

  result->nRefs=1;

  return result;
 }


/*! \fn void IMFWriteExpansion(FILE* fid, IMFExpansion E, MFErrorHandler e);
 *  \brief Writes an expansion to a file.
 *
 *  \param fid The file to write to.
 *  \param E The expansion being queried.
 *   \param e A place to handle exceptions.
 */
void IMFWriteExpansion(FILE* fid,IMFExpansion E,MFErrorHandler e);

/*! \fn IMFExpansion IMFReadExpansion(FILE* fid,MFErrorHandler e);
 *  \brief Reads an expansion from a file.
 *
 *  \param fid The file to write to.
 *  \param e A place to handle exceptions.
 *  \returns The expansion.
 */
IMFExpansion IMFReadExpansion(FILE* fid,MFErrorHandler e);


/*! \fn void IMFFreeExpansion(IMFExpansion E, MFErrorHandler e);
 *  \brief Frees a reference to an expansion, and deletes the expansion if there are no references left.
 *
 *  \param E The expansion being unreferenced.
 *  \param e A place to handle exceptions.
 *  \sa ReferenceCounting MFRefExpansion(MFErrorHandler)
 */
void IMFFreeExpansion(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFFreeExpansion"};

  E->nRefs--;

  if(E->nRefs<=0)
   {
    if(E->data!=NULL){free(E->data);E->data=NULL;}
    E->dataLn=0;
    E->u=NULL;
    E->du=NULL;
    E->ddu=NULL;
    E->dddu=NULL;
    E->n=0;
    E->k=0;

    free(E);
   }

  return;
 }

/*! \fn void IMFRefExpansion(IMFExpansion E, MFErrorHandler e);
 *  \brief Adds a reference to a expansion.
 *
 *  \param E The expansion being referenced.
 *  \param e A place to handle exceptions.
 *  \sa ReferenceCounting MFFreeExpansion()
 */
void IMFRefExpansion(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFRefExpansion"};

  E->nRefs++;
  return;
 }

/*! \fn void IMFEvaluateExpansion(IMFExpansion E,MFKVector s,MFNVector u, MFErrorHandler e);
 *  \brief Evaluates an expansion.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param u Provided by the user, a place to put the result.
 *  \param e A place to handle exceptions.
 */
void IMFEvaluateExpansion(IMFExpansion E,MFKVector s,MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateExpansion"};
  double *ds;
  double *du;
  int i,p,q,r;

  ds=MFKV_CStar(s,e);
  du=MFNV_CStar(u,e);

  for(i=0;i<E->n;i++)
   {
    du[i]=E->u[i];
    if(E->du!=NULL)
     {
      for(p=0;p<E->k;p++)
       {
        du[i]+=E->du[i+E->n*p]*ds[p];
        if(E->ddu!=NULL)
         {
          for(q=0;q<E->k;q++)
           {
            du[i]+=.5*E->ddu[i+E->n*(p+E->k*q)]*ds[p]*ds[q];
            if(E->dddu!=NULL)
             {
              for(r=0;r<E->k;r++)
               {
                du[i]+=E->dddu[i+E->n*(p+E->k*(q+E->k*r))]*ds[p]*ds[q]*ds[r]/6.;
               }
             }
           }
         }
       }
     }
   }
  return;
 }

/*! \fn void IMFEvaluateExpansionDirectionalDerivative(IMFExpansion E,MFKVector s,MFKVector ds,MFNVector du, MFErrorHandler e)
 *  \brief Evaluates an expansion.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param ds A direction in which to find the partial derivative.
 *  \param du Provided by the user, a place to put the result.
 *  \param e A place to handle exceptions.
 */
void IMFEvaluateExpansionDirectionalDerivative(IMFExpansion E,MFKVector s,MFKVector ds,MFNVector du, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateExpansionDirectionalDerivative"};
  double *As;
  double *Ads;
  double *Adu;
  int i,p,q,r;

  As=MFKV_CStar(s,e);
  Ads=MFKV_CStar(ds,e);
  Adu=MFNV_CStar(du,e);

  for(i=0;i<E->n;i++)
   {
    Adu[i]=0.;
    if(E->du!=NULL)
     {
      for(p=0;p<E->k;p++)
       {
        Adu[i]+=E->du[i+E->n*p]*Ads[p];
        if(E->ddu!=NULL)
         {
          for(q=0;q<E->k;q++)
           {
            Adu[i]+=E->ddu[i+E->n*(p+E->k*q)]*Ads[p]*As[q];
            if(E->dddu!=NULL)
             {
              for(r=0;r<E->k;r++)
                Adu[i]+=.5*E->dddu[i+E->n*(p+E->k*(q+E->k*r))]*Ads[p]*As[q]*As[r];
             }
           }
         }
       }
     }
   }
  return;
 }

/*! \fn void IMFEvaluateExpansionSecondDirectionalDerivative(IMFExpansion E,MFKVector s,MFKVector ds0,MFKVector ds1,MFNVector ddu, MFErrorHandler e);
 *  \brief Evaluates an expansion.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param ds0 The first direction.
 *  \param ds1 The second direction.
 *  \param ddu Provided by the user, a place to put the result.
 *  \param e A place to handle exceptions.
 */
void IMFEvaluateExpansionSecondDirectionalDerivative(IMFExpansion E,MFKVector s,MFKVector ds0,MFKVector ds1,MFNVector ddu, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateExpansionSecondDirectionalDerivative"};
  double *As;
  double *Ads0,*Ads1;
  double *Addu;
  int i,p,q,r;

  As=MFKV_CStar(s,e);
  Ads0=MFKV_CStar(ds0,e);
  Ads1=MFKV_CStar(ds1,e);
  Addu=MFNV_CStar(ddu,e);

  for(i=0;i<E->n;i++)
   {
    Addu[i]=0.;
    if(E->du!=NULL)
     {
      for(p=0;p<E->k;p++)
       {
        if(E->ddu!=NULL)
         {
          for(q=0;q<E->k;q++)
           {
            Addu[i]+=E->ddu[i+E->n*(p+E->k*q)]*Ads0[p]*Ads1[q];

            if(E->dddu!=NULL)
             {
              for(r=0;r<E->k;r++)
                Addu[i]+=E->dddu[i+E->n*(p+E->k*(q+E->k*r))]*Ads0[p]*Ads1[q]*As[r];
             }
           }
         }
       }
     }
   }
  return;
 }

/*! \fn void IMFEvaluateExpansionDerivative(IMFExpansion E,MFKVector s,int d,MFNVector du, MFErrorHandler e)
 *  \brief Evaluates the derivative of an expansion along a coordinate direction.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param d A coordinate axis in which to find the partial derivative.
 *  \param du Provided by the user, a place to put the result.
 *  \param e A place to handle exceptions.
 */
void IMFEvaluateExpansionDerivative(IMFExpansion E,MFKVector s,int d,MFNVector du, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateExpansionDerivative"};
  double *As;
  double *Adu;
  int i,p,q;

  As=MFKV_CStar(s,e);
  Adu=MFNV_CStar(du,e);

  if(d<0||d>E->k){fprintf(stderr,"d out of range in EvaluateExpansion Derivative, 0<=%d<%d\n",d,E->k);fflush(stderr);MFErrorHandlerOutOfMemory(e);}

  for(i=0;i<E->n;i++)
   {
    Adu[i]=0.;
    if(E->du!=NULL)  
     {
      Adu[i]+=E->du[i+E->n*d];
      for(p=0;p<E->k;p++)
       {
        if(E->ddu!=NULL)
         {
          Adu[i]+=E->ddu[i+E->n*(p+E->k*d)]*As[p];
          if(E->dddu!=NULL)
           {
            for(q=0;q<E->k;q++)
             {
              Adu[i]+=E->dddu[i+E->n*(p+E->k*(q+E->k*d))]*As[p]*As[q]/2.;
             }
           }
         }
       }
     }
   }
  return;
 }

/*! \fn void IMFEvaluateExpansionSecondDerivative(IMFExpansion E,MFKVector s,int d0,int d1,MFNVector ddu, MFErrorHandler e)
 *  \brief Evaluates the second derivative of an expansion along coordinate directions.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param d0 A coordinate axis in which to find the partial derivative.
 *  \param d1 A second coordinate axis in which to find the partial derivative.
 *  \param ddu Provided by the user, a place to put the result.
 *  \param e A place to handle exceptions.
 */
void IMFEvaluateExpansionSecondDerivative(IMFExpansion E,MFKVector s,int d0,int d1,MFNVector ddu, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateExpansionSecondDerivative"};
  double *As;
  double *Addu;
  int i,p;

  As=MFKV_CStar(s,e);
  Addu=MFNV_CStar(ddu,e);

  for(i=0;i<E->n;i++)
   {
    Addu[i]=0.;
    if(E->du!=NULL && E->ddu!=NULL)
     {
      Addu[i]=E->ddu[i+E->n*(d0+E->k*d1)];
      if(E->dddu!=NULL)
       {
        for(p=0;p<E->k;p++)
          Addu[i]+=E->dddu[i+E->n*(p+E->k*(d0+E->k*d1))]*As[p];
       }
     }
   }
  return;
 }

/*! \fn void IMFEvaluateExpansionThirdDerivative(IMFExpansion E,MFKVector s,int d0,int d1,int d2,MFNVector dddu, MFErrorHandler e);
 *  \brief Evaluates the third derivative of an expansion along coordinate directions.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param d0 A coordinate axis in which to find the partial derivative.
 *  \param d1 A second coordinate axis in which to find the partial derivative.
 *  \param d2 A third coordinate axis in which to find the partial derivative.
 *  \param dddu Provided by the user, a place to put the result.
 *  \param e A place to handle exceptions.
 */
void IMFEvaluateExpansionThirdDerivative(IMFExpansion E,MFKVector s,int d0,int d1,int d2,MFNVector dddu, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateExpansionThirdDerivative"};
  double *As;
  double *Adddu;
  int i;

  As=MFKV_CStar(s,e);
  Adddu=MFNV_CStar(dddu,e);

  for(i=0;i<E->n;i++)
   {
    if(E->du!=NULL && E->ddu!=NULL&&E->dddu!=NULL)
      Adddu[i]=E->dddu[i+E->n*(d0+E->k*(d1+E->k*d2))];
     else
      Adddu[i]=0.;
   }
  return;
 }

/*! \fn int IMFExpansionN(IMFExpansion E, MFErrorHandler e)
 *  \brief Returns the dimension of the range of an expansion.
 *
 *  \param E An expansion.
 *  \param e A place to handle exceptions.
 *  \returns The dimension of the range of the expansion.
 */
int IMFExpansionN(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionN"};

  return E->n;
 }

/*! \fn int IMFExpansionK(IMFExpansion E, MFErrorHandler e)
 *  \brief Returns the dimension of the domain of an expansion.
 *
 *  \param E An expansion.
 *  \param e A place to handle exceptions.
 *  \returns The dimension of the domain of the expansion.
 */
int IMFExpansionK(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionK"};

  return E->k;
 }

/*! \fn int IMFExpansionOrder(IMFExpansion E, MFErrorHandler e)
 *  \brief Returns the order of an expansion (the degree of the highest non-zero term in the Taylor series).
 *
 *  \param E An expansion.
 *  \param e A place to handle exceptions.
 *  \returns The order of the expansion.
 */
int IMFExpansionOrder(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionOrder"};
  int O;

  O=0;
  if(E->du!=NULL)
   {
    O=1;
    if(E->ddu!=NULL)
     {O=2;if(E->dddu!=NULL)O=3;}
   }

  return O;
 }

double *IMFExpansionU(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionU"};

  return E->u;
 }

double *IMFExpansionDu(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionDu"};

  return E->du;
 }

double *IMFExpansionDDu(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionDDu"};

  return E->ddu;
 }

double *IMFExpansionDDDu(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionDDDu"};

  return E->dddu;
 }

double *IMFExpansionData(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionData"};

  return E->data;
 }

int IMFExpansionDataLn(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionDataLn"};

  return E->dataLn;
 }

void IMFExpansionSetDerivatives(IMFExpansion E,double *u,
                                               double *du,
                                               double *ddu,
                                               double *dddu, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionSetDerivatives"};
  int i;

  if(E->data!=NULL)free(E->data);

  E->dataLn=E->n;
  if(du!=NULL)
   {
    E->dataLn+=E->n*E->k;
    if(ddu!=NULL)
     {
      E->dataLn+=E->n*E->k*E->k;
      if(dddu!=NULL)
        E->dataLn+=E->n*E->k*E->k*E->k;
     }
   }

  E->data=(double*)malloc(E->dataLn*sizeof(double));

#ifndef MFNOSAFETYNET
  if(E->data==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",E->dataLn*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  if(E->data==NULL && E->dataLn!=0)
   {
    MFErrorHandlerOutOfMemory(e);
   }

  E->u=E->data;
  if(du!=NULL)
   {
    E->du=E->u+E->n;
    if(ddu!=NULL)
     {
      E->ddu=E->du+E->n*E->k;
      if(dddu!=NULL)
       {
        E->dddu=E->ddu+E->n*E->k*E->k;
       }else{
        E->dddu=NULL;
       }
     }else{
      E->ddu=NULL;
      E->dddu=NULL;
     }
   }else{
    E->du=NULL;
    E->ddu=NULL;
    E->dddu=NULL;
   }

  if(u!=NULL)
   {
    for(i=0;i<E->n;i++)(E->u)[i]=u[i];
    if(du!=NULL)
     {
      for(i=0;i<E->n*E->k;i++)(E->du)[i]=du[i];
      if(ddu!=NULL)
       {
        for(i=0;i<E->n*E->k*E->k;i++)(E->ddu)[i]=ddu[i];
        if(dddu!=NULL)
         {
          for(i=0;i<E->n*E->k*E->k*E->k;i++)(E->dddu)[i]=dddu[i];
         }
       }
     }
   }

  return;
 }

char *IMFC(double a, char *suffix,int *leading)
 {
  static char RoutineName[]={"IMFC"};
  static char result[1024];
  int ia;

  result[0]=0x0;
  ia=round(a);
  if(suffix[0]!=0x0 && !(*leading))
   {
    if(fabs(a-1)<1.e-7)
      sprintf(result,"+%s",suffix);
     else if(fabs(a+1)<1.e-7)
      sprintf(result,"-%s",suffix);
     else if(fabs(a-ia)<1.e-7 && ia > 0)
      sprintf(result,"+ %d*%s",ia,suffix);
     else if(fabs(a-ia)<1.e-7 && ia < 0)
      sprintf(result,"- %d*%s",-ia,suffix);
     else if(fabs(a-ia-1)<1.e-7 && ia+1 > 0)
      sprintf(result,"+ %d*%s",ia+1,suffix);
     else if(fabs(a-ia-1)<1.e-7 && ia+1 < 0)
      sprintf(result,"- %d*%s",-ia-1,suffix);
     else if(fabs(a-ia+1)<1.e-7 && ia-1 > 0)
      sprintf(result,"+ %d*%s",ia-1,suffix);
     else if(fabs(a-ia+1)<1.e-7 && ia-1 < 0)
      sprintf(result,"- %d*%s",-ia+1,suffix);
     else if(a>1.e-7)
      sprintf(result,"+ %le*%s",a,suffix);
     else if(a<-1.e-7)
      sprintf(result,"- %le*%s",-a,suffix);
     else
      sprintf(result,"");
   }else if(suffix[0]==0x0 && !(*leading))
   {
    if(fabs(a-1)<1.e-7)
      sprintf(result,"+1");
     else if(fabs(a+1)<1.e-7)
      sprintf(result,"-1");
     else if(fabs(a-ia)<1.e-7 && ia > 0)
      sprintf(result,"+ %d",ia);
     else if(fabs(a-ia)<1.e-7 && ia < 0)
      sprintf(result,"- %d",-ia);
     else if(fabs(a-ia-1)<1.e-7 && ia+1 > 0)
      sprintf(result,"+ %d",ia+1);
     else if(fabs(a-ia-1)<1.e-7 && ia+1 < 0)
      sprintf(result,"- %d",-ia-1);
     else if(fabs(a-ia+1)<1.e-7 && ia-1 > 0)
      sprintf(result,"+ %d",ia-1);
     else if(fabs(a-ia+1)<1.e-7 && ia-1 < 0)
      sprintf(result,"- %d",-ia+1);
     else if(a>1.e-7)
      sprintf(result,"+ %le",a);
     else if(a<-1.e-7)
      sprintf(result,"- %le",-a);
     else
      sprintf(result,"");
   }else if(suffix[0]!=0x0 && *leading)
   {
    if(fabs(a-1)<1.e-7)
      {sprintf(result,"%s",suffix);*leading=0;}
     else if(fabs(a+1)<1.e-7)
      {sprintf(result,"-%s",suffix);*leading=0;}
     else if(fabs(a-ia)<1.e-7 && ia > 0)
      {sprintf(result," %d*%s",ia,suffix);*leading=0;}
     else if(fabs(a-ia)<1.e-7 && ia < 0)
      {sprintf(result,"- %d*%s",-ia,suffix);*leading=0;}
     else if(fabs(a-ia-1)<1.e-7 && ia+1 > 0)
      {sprintf(result," %d*%s",ia+1,suffix);*leading=0;}
     else if(fabs(a-ia-1)<1.e-7 && ia+1 < 0)
      {sprintf(result,"- %d*%s",-ia-1,suffix);*leading=0;}
     else if(fabs(a-ia+1)<1.e-7 && ia-1 > 0)
      {sprintf(result," %d*%s",ia-1,suffix);*leading=0;}
     else if(fabs(a-ia+1)<1.e-7 && ia-1 < 0)
      {sprintf(result,"- %d*%s",-ia+1,suffix);*leading=0;}
     else if(a>1.e-7)
      {sprintf(result," %le*%s",a,suffix);*leading=0;}
     else if(a<-1.e-7)
      {sprintf(result,"- %le*%s",-a,suffix);*leading=0;}
     else
      sprintf(result,"");
   }else if(suffix[0]==0x0 && *leading)
   {
    if(fabs(a-1)<1.e-7)
      {sprintf(result,"1");*leading=0;}
     else if(fabs(a+1)<1.e-7)
      {sprintf(result,"-1");*leading=0;}
     else if(fabs(a-ia)<1.e-7 && ia > 0)
      {sprintf(result," %d",ia);*leading=0;}
     else if(fabs(a-ia)<1.e-7 && ia < 0)
      {sprintf(result,"- %d",-ia);*leading=0;}
     else if(fabs(a-ia-1)<1.e-7 && ia+1 > 0)
      {sprintf(result," %d",ia+1);*leading=0;}
     else if(fabs(a-ia-1)<1.e-7 && ia+1 < 0)
      {sprintf(result,"- %d",-ia-1);*leading=0;}
     else if(fabs(a-ia+1)<1.e-7 && ia-1 > 0)
      {sprintf(result," %d",ia-1);*leading=0;}
     else if(fabs(a-ia+1)<1.e-7 && ia-1 < 0)
      {sprintf(result,"- %d",-ia+1);*leading=0;}
     else if(a>1.e-7)
      {sprintf(result," %le",a);*leading=0;}
     else if(a<-1.e-7)
      {sprintf(result,"- %le",-a);*leading=0;}
     else
      sprintf(result,"");
   }
  return result;
 }

/*! \fn void IMFPrintExpansion(FILE *fid,IMFExpansion E, MFErrorHandler e);
 *  \brief Prints a \"pretty\" version of the expansion.
 *
 *  \param fid The file.
 *  \param E The expansion.
 *  \param e A place to handle exceptions.
 */
void IMFPrintExpansion(FILE *fid,IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFPrintExpansion"};
  int i;
  int leading;

  fprintf(fid,"IMFPrintExpansion, data=0x%8.8x, length %d\n",E->data,E->dataLn);
  fprintf(fid,"                   u=0x%8.8x\n",E->u);
  fprintf(fid,"                  du=0x%8.8x\n",E->du);
  fprintf(fid,"                 ddu=0x%8.8x\n",E->ddu);
  fprintf(fid,"                dddu=0x%8.8x\n",E->dddu);

  switch(E->k)
   {
    case 1:
     printf("u(s) = (");
     for(i=0;i<E->n;i++)
      {
       leading=1;
       if(i>0)fprintf(fid,"       ");
       fprintf(fid,"%s",IMFC(E->u[i],"",&leading));
       fprintf(fid,"%s",IMFC(E->du[i],"s",&leading));
       fprintf(fid,"%s",IMFC(E->ddu[i]/2.,"s*s",&leading));
       fprintf(fid,"%s",IMFC(E->dddu[i]/6.,"s*s*s",&leading));
       if(i<E->n-1)fprintf(fid,",\n");
        else fprintf(fid,") + O(s^4)\n");
      }

     break;
    case 2:
     fprintf(fid,"u(s,t) = (");
     for(i=0;i<E->n;i++)
      {
       leading=1;
       if(i>0)fprintf(fid,"         ");
       fprintf(fid,"%s",IMFC(E->u[i],"",&leading));
       if(E->du!=NULL)
        {
         fprintf(fid,"%s",IMFC(E->du[i],"s",&leading));
         fprintf(fid,"%s",IMFC(E->du[i+E->n],"t",&leading));
         if(E->ddu!=NULL)
          {
           fprintf(fid,"%s",IMFC(E->ddu[i]/2,"s*s",&leading));
           fprintf(fid,"%s",IMFC(E->ddu[i+E->n*(1+E->k*0)],"s*t",&leading));
           fprintf(fid,"%s",IMFC(E->ddu[i+E->n*(1+E->k*1)],"t*t",&leading));
           if(E->dddu!=NULL)
            {
             fprintf(fid,"%s",IMFC(E->dddu[i+E->n*(0+E->k*(0+E->k*0))]/6,"s*s*s",&leading));
             fprintf(fid,"%s",IMFC(E->dddu[i+E->n*(0+E->k*(0+E->k*1))]/2,"s*s*t",&leading));
             fprintf(fid,"%s",IMFC(E->dddu[i+E->n*(0+E->k*(1+E->k*1))]/2,"s*t*t",&leading));
             fprintf(fid,"%s",IMFC(E->dddu[i+E->n*(1+E->k*(1+E->k*1))]/6,"t*t*t",&leading));
            }
          }
        }
       if(i<E->n-1)fprintf(fid,",\n");
        else fprintf(fid,") + O(s^4,s^3t,s^2t^2,s t^3,t^4)\n");
      }
     break;
   }
  return;
 }

/*! \fn IMFExpansion IMFCloneExpansion(IMFExpansion E, MFErrorHandler e);
 *  \brief Creates a deep copy of an expansion.
 *
 *  \param E The expansion.
 *  \param e A place to handle exceptions.
 *  \returns A clone.
 */
IMFExpansion IMFCloneExpansion(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCloneExpansion"};
  IMFExpansion result;
  int i;

  result=IMFCreateExpansion(E->n,E->k,e);
  result->data=(double*)malloc(E->dataLn*sizeof(double));

#ifndef MFNOSAFETYNET
  if(result->data==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",E->dataLn*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->dataLn=E->dataLn;
  result->u=NULL;
  result->du=NULL;
  result->ddu=NULL;
  result->dddu=NULL;

  if(E->dataLn>0&&E->u!=NULL)
   {
    result->u=result->data;
    for(i=0;i<E->n;i++)(result->u)[i]=(E->u)[i];
   }
  if(E->dataLn>0&&E->du!=NULL)
   {
    result->du=result->u+result->n;
    for(i=0;i<E->n*E->k;i++)(result->du)[i]=(E->du)[i];
   }
  if(E->dataLn>0&&E->ddu!=NULL)
   {
    result->ddu=result->du+result->n*result->k;
    for(i=0;i<E->n*E->k*E->k;i++)(result->ddu)[i]=(E->ddu)[i];
   }
  if(E->dataLn>0&&E->dddu!=NULL)
   {
    result->dddu=result->ddu+result->n*result->k*result->k;
    for(i=0;i<E->n*E->k*E->k*E->k;i++)(result->dddu)[i]=(E->dddu)[i];
   }

  return result;
 }

/*! \fn IMFExpansion IMFRectify(IMFExpansion E0, MFErrorHandler e)
 *  \brief Changes the metric of the expansion to the identity. 
 *
 *  \param fid The file.
 *  \param E The expansion.
 *  \param e A place to handle exceptions.
 */
IMFExpansion IMFRectify(IMFExpansion E0, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFRectify"};
  double t;

  double *g0;
  double *G0;

  double *G;

  double *ng;
  double *nG;

  double *dX;
  double *ddX;

  double *A;
  double *b;
  int i,j,k,p,q,r,w;
  int verbose=0;
  IMFExpansion E;
  int n,m;

  double *u;
  double *du;
  double *ddu;
  double *dddu;

  n=E0->n;
  m=E0->k;

  g0=(double*)malloc(m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(g0==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  G0=(double*)malloc(m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(G0==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  G =(double*)malloc(m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(G==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  ng=(double*)malloc(m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ng==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  nG=(double*)malloc(m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(nG==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  dX=(double*)malloc(m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dX==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  ddX=(double*)malloc(m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ddX==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  A=(double*)malloc(m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  b=(double*)malloc(m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(b==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  u=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(u==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  du=(double*)malloc(n*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(du==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  ddu=(double*)malloc(n*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ddu==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  dddu=(double*)malloc(n*m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dddu==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("IMFRectify n=%d, k=%d\n",n,m);fflush(stdout);}
#endif

  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      g0[i+m*j]=0.;for(p=0;p<n;p++)g0[i+m*j]+=E0->du[p+n*i]*E0->du[p+n*j];
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("\n");
    printf("Metric in old coordinates\n\n");
    for(i=0;i<m;i++)
     {
      printf("   [%lf",g0[i+m*0]);
      for(j=1;j<m;j++)printf(" %lf",g0[i+m*j]);
      printf("]\n");
     }
    if(m==2)printf(" determinant %le\n",g0[0+m*0]*g0[1+m*1]-g0[1+m*0]*g0[0+m*1]);
    fflush(stdout);
   }
#endif

/* Gram-Schmidt */

  for(j=0;j<m;j++)
   {
    for(i=0;i<n;i++)du[i+n*j]=E0->du[i+n*j];

    for(k=0;k<j;k++)
     {
      t=0.;for(i=0;i<n;i++)t+=du[i+n*j]*du[i+n*k];
      for(i=0;i<n;i++)du[i+n*j]=du[i+n*j]-t*du[i+n*k];
     }
  
    t=0.;for(i=0;i<n;i++)t+=du[i+n*j]*du[i+n*j];t=1./sqrt(t);
    for(i=0;i<n;i++)du[i+n*j]=t*du[i+n*j];
   }

  for(j=0;j<m;j++)
   {
    for(k=0;k<m;k++)
     {
      t=0.;for(i=0;i<n;i++)t+=du[i+n*j]*du[i+n*k];
      ng[j+m*k]=t;
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("\n");
    printf("Metric in new coordinates\n\n");
    for(i=0;i<m;i++)
     {
      printf("   [%lf",ng[i+m*0]);
      for(j=1;j<m;j++)printf(" %lf",ng[i+m*j]);
      printf("]\n");
     }
    if(m==2)printf(" new determinant %le\n",ng[0+m*0]*ng[1+m*1]-ng[1+m*0]*ng[0+m*1]);
    fflush(stdout);
   }
#endif

  for(k=0;k<m;k++)
   {
    for(i=0;i<m;i++)for(j=0;j<m;j++)A[i+m*j]=g0[j+m*i];
    for(j=0;j<m;j++){b[j]=0.;for(p=0;p<n;p++)b[j]+=E0->du[p+n*j]*du[p+n*k];}
    if(!MFSolveFull(m,A,b,e))exit(12);

    for(j=0;j<m;j++)dX[j+m*k]=b[j];
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("\n");
    printf("Check of first derivative of coordinate transformation X^i_{,j}\n\n");
    for(j=0;j<m;j++)
     {
      printf("   u_{,%d}             = [%9.6lf",j,du[0+n*j]);
      for(i=1;i<n;i++)printf(",%9.6lf",du[i+n*j]);
      printf("]\n");
  
      printf("   X^p_{,%d} u_{,\\xi_p}= [",j);
      for(i=0;i<n;i++)
       {
        t=0.;for(p=0;p<m;p++)t+=dX[p+m*j]*E0->du[i+n*p];
        if(i>0)printf(",");
        printf("%9.6lf",t);
       }
      printf("]\n");
     }
    fflush(stdout);
   }
#endif

  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        G0[i+m*(j+m*k)]=0.;for(p=0;p<n;p++)G0[i+m*(j+m*k)]+=E0->du[p+n*i]*E0->ddu[p+n*(j+m*k)];
       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("\n");
    printf("Connection in old coordinates\n\n");
    for(i=0;i<m;i++)
     {
      for(j=0;j<m;j++)
       {
        if(j==0)printf("  %d [",i);
         else printf("    [");
        for(k=0;k<m;k++)
         {
          if(k>0)printf(" ");
          printf("%9.6lf",G0[i+m*(j+m*k)]);
         }
        printf("]\n");
       }
      printf("\n");
     }
    fflush(stdout);
   }
#endif

  for(i=0;i<m*m*m;i++)G[i]=0.;

  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        ddX[i+m*(j+m*k)]=0.;
        for(p=0;p<m;p++)
         {
          ddX[i+m*(j+m*k)]+=dX[i+m*p]*G[p+m*(j+m*k)];
          for(q=0;q<m;q++)
           {
            for(r=0;r<m;r++)
             {
              for(w=0;w<m;w++)
               {
                ddX[i+m*(j+m*k)]-=dX[i+m*w]*dX[p+m*w]*dX[q+m*j]*dX[r+m*k]*G0[p+m*(q+m*r)];
               }
             }
           }
         }
       }
     }
   }


  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        nG[i+m*(j+m*k)]=0.;
        for(p=0;p<m;p++)
         {
          for(q=0;q<m;q++)
           {
            nG[i+m*(j+m*k)]+=dX[p+m*i]*ddX[q+m*(j+m*k)]*g0[p+m*q];
            for(r=0;r<m;r++)
             {
              nG[i+m*(j+m*k)]+=dX[p+m*i]*dX[q+m*j]*dX[r+m*k]*G0[p+m*(q+m*r)];
             }
           }
         }
       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("\n");
    printf("Connection in new coordinates (requested) \n\n");
    for(i=0;i<m;i++)
     {
      for(j=0;j<m;j++)
       {
        if(j==0)printf("  %d [",i);
         else printf("    [");
        for(k=0;k<m;k++)
         {
          if(k>0)printf(" ");
          printf("%9.6lf",G[i+m*(j+m*k)]);
         }
        printf("]\n");
       }
      printf("\n");
     }
  
    printf("\n");
    printf("Connection in new coordinates (transformed)\n\n");
    for(i=0;i<m;i++)
     {
      for(j=0;j<m;j++)
       {
        if(j==0)printf("  %d [",i);
         else printf("    [");
        for(k=0;k<m;k++)
         {
          if(k>0)printf(" ");
          printf("%9.6lf",nG[i+m*(j+m*k)]);
         }
        printf("]\n");
       }
      printf("\n");
     }
    fflush(stdout);
   }
#endif
  
  for(i=0;i<n;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        ddu[i+n*(j+m*k)]=0.;
        for(p=0;p<m;p++)
         {
          ddu[i+n*(j+m*k)]+=ddX[p+m*(j+m*k)]*E0->du[i+n*p];
          for(q=0;q<m;q++)
           ddu[i+n*(j+m*k)]+=dX[p+m*j]*dX[q+m*k]*E0->ddu[i+n*(p+m*q)];
         }
       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("\n");
    printf("Connection in new coordinates (direct)\n\n");
    for(i=0;i<m;i++)
     {
      for(j=0;j<m;j++)
       {
        if(j==0)printf("  %d [",i);
         else printf("    [");
        for(k=0;k<m;k++)
         {
          if(k>0)printf(" ");
          t=0.;for(p=0;p<n;p++)t+=du[p+n*i]*ddu[p+n*(j+m*k)];
          printf("%9.6lf",t);
         }
        printf("]\n");
       }
      printf("\n");
     }
    fflush(stdout);
   }
#endif

/* Transform third derivatives too ? */

  free(g0);
  free(G0);
  free(ng);
  free(G);
  free(nG);
  free(dX);
  free(ddX);
  free(A);
  free(b);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done IMFRectify\n");fflush(stdout);fflush(stdout);}
#endif

  E=IMFCreateExpansion(n,m,e);
  for(i=0;i<n;i++)u[i]=(E0->u)[i];
  IMFExpansionSetDerivatives(E,u,du,ddu,NULL,e);
  free(u);
  free(du);
  free(ddu);
  free(dddu);

  return E;
 }

/*! \fn IMFExpansion IMFInflateExpansionWithFlow(IMFExpansion E,IMFFlow f, MFKVector p, MFErrorHandler e)
 *  \brief Creates the product of an expansion with a trjectory of a flow.
 *
 *  \param E The starting expansion.
 *  \param f The flow. (At the origin of the expansion should be transverse to the surface defined by the expansion).
 *  \param p The flow's parameter values for the trajectory.
 *  \param e A place to handle exceptions.
 *  \returns A new expansion.
 */
IMFExpansion IMFInflateExpansionWithFlow(IMFExpansion E,IMFFlow f, MFKVector p, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFInflateExpansionWithFlow"};
  double *u;
  double *du;
  double *ddu;
  double *F;
  double *dF;
  double d;
  IMFExpansion tmp;
  IMFExpansion result;
  int n,m,k;
  int i,j,l;
  MFNVector v;

  n=IMFExpansionN(E,e);
  k=IMFExpansionK(E,e);
  m=k+1;

  tmp=IMFCreateExpansion(n,m,e);
  u=E->u;
  du=(double*)malloc(n*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(du==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  ddu=(double*)malloc(n*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ddu==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  F=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(F==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  dF=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dF==NULL)
   {
    sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  v=MFCreateWrappedNVector(n,E->u,e);
  IMFEvaluateFlow(f,v,p,F,e);
  IMFEvaluateDerivativeOfFlow(f,v,p,dF,e);
  MFFreeNVector(v,e);

  d=0.;for(i=0;i<n;i++)d+=F[i]*F[i];d=sqrt(d);
  for(i=0;i<n;i++)F[i]=F[i]/d;

  for(i=0;i<n;i++)
   {
    for(j=0;j<k;j++)
     {
      du[i+n*j]=E->du[i+n*j];
     }
    du[i+n*k]=F[i];
   }

  for(i=0;i<n;i++)
   {
    for(j=0;j<k;j++)
     {
      for(l=0;l<k;l++)
       {
        ddu[i+n*(j+m*l)]=E->ddu[i+n*(j+k*l)];
       }
      ddu[i+n*(j+m*k)]=0.;
      ddu[i+n*(k+m*j)]=0.;
     }
    ddu[i+n*(k+m*k)]=0.;
    for(j=0;j<n;j++)ddu[i+n*(k+m*k)]+=dF[i+n*j]*F[j]/d; /* Hadn't divided by d before REMOVEME?? */
   }

  IMFExpansionSetDerivatives(tmp,u,du,ddu,NULL,e);
  free(du);
  free(ddu);
  free(F);
  free(dF);
  result=IMFRectify(tmp,e);
  IMFFreeExpansion(tmp,e);
  return result;
 }

/*! \fn double IMFExpansionR(IMFExpansion E, double epsilon, MFErrorHandler e)
 *  \brief Estimates a ball such that within the ball the second order terms are bounded in absolute value by epsilon.
 *
 *  \param fid The file.
 *  \param E The expansion.
 *  \returns The radius of the ball.
 */
double IMFExpansionR(IMFExpansion E, double epsilon, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionR"};
  int i,j;
  double A;
  double t;

  A=0.;
  for(j=0;j<E->k*E->k;j++)
   {
    t=0.;
    for(i=0;i<E->n;i++)t+=pow(E->ddu[i+E->n*j],2.);
    t=sqrt(t);
    if(t<A||j==0)A=t;
   }
  return sqrt(2*epsilon/A);
 }

/*! \fn MFNKMatrix IMFExpansionTS(IMFExpansion E, MFErrorHandler e)
 *  \brief Returns a new MFNKMatrix with an orthonormal basis for the tangent space at the origin of the expansion.
 *
 *  \param E The expansion.
 *  \param e A place to handle exceptions.
 *  \returns The basis in an MFNKMatrix.
 */
MFNKMatrix IMFExpansionTS(IMFExpansion E, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionTS"};

  return MFCreateNKMatrixWithData(E->n,E->k,E->du,e);
 }

void IMFOLDExpansionSetDerivatives(IMFExpansion E,double *u,
                                               double *du,
                                               double *ddu,
                                               double *dddu, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFOLDExpansionSetDerivatives"};
  int i;

  if(E->u!=NULL)free(E->u);
  if(u!=NULL)
   {
    E->u=(double*)malloc(E->n*sizeof(double));

#ifndef MFNOSAFETYNET
    if(E->u==NULL)
     {
      sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",E->n*sizeof(double));
      MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=0;i<E->n;i++)(E->u)[i]=u[i];
   }else
    E->u=NULL;

  if(E->du!=NULL)free(E->du);
  if(du!=NULL)
   {
    E->du=(double*)malloc(E->n*E->k*sizeof(double));if(E->du==NULL)MFErrorHandlerOutOfMemory(e);

#ifndef MFNOSAFETYNET
    if(E->du==NULL)
     {
      sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",E->n*E->k*sizeof(double));
      MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=0;i<E->n*E->k;i++)(E->du)[i]=du[i];
   }else
    E->du=NULL;

  if(E->ddu!=NULL)free(E->ddu);
  if(ddu!=NULL)
   {
    E->ddu=(double*)malloc(E->n*E->k*E->k*sizeof(double));if(E->ddu==NULL)MFErrorHandlerOutOfMemory(e);

#ifndef MFNOSAFETYNET
    if(E->ddu==NULL)
     {
      sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",E->n*E->k*E->k*sizeof(double));
      MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=0;i<E->n*E->k*E->k;i++)(E->ddu)[i]=ddu[i];
   }else
    E->ddu=NULL;

  if(E->dddu!=NULL)free(E->dddu);
  if(dddu!=NULL)
   {
    E->dddu=(double*)malloc(E->n*E->k*E->k*E->k*sizeof(double));if(E->dddu==NULL)MFErrorHandlerOutOfMemory(e);

#ifndef MFNOSAFETYNET
    if(E->dddu==NULL)
     {
      sprintf(IMFExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",E->n*E->k*E->k*E->k*sizeof(double));
      MFSetError(e,12,RoutineName,IMFExpansionErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=0;i<E->n*E->k*E->k*E->k;i++)(E->dddu)[i]=dddu[i];
   }else
    E->dddu=NULL;

  return;
 }

void IMFEvaluateExpansionDirectionalDerivativeE(IMFExpansion E,double *s,double *ds,double *du, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateExpansionDirectionalDerivativeE"};
  int i,p,q,r;
  int j,k,l,m;

  for(i=0;i<E->n;i++)
   {
    du[i]=0.;
    if(E->du!=NULL)
     {
      for(p=0;p<E->k;p++)
       {
        du[i]+=E->du[i+E->n*p]*ds[p];
        if(E->ddu!=NULL)
         {
          for(q=0;q<E->k;q++)
           {
            du[i]+=E->ddu[i+E->n*(p+E->k*q)]*ds[p]*s[q];
            if(E->dddu!=NULL)
             {
              for(r=0;r<E->k;r++)
                du[i]+=.5*E->dddu[i+E->n*(p+E->k*(q+E->k*r))]*ds[p]*s[q]*s[r];
             }
           }
         }
       }
     }

    for(j=0;j<E->k;j++)
     {
      du[i+E->n*j]=0.;
      if(E->du!=NULL)
       {
        for(p=0;p<E->k;p++)
         {
          if(E->ddu!=NULL)
           {
            du[i+E->n*j]+=E->ddu[i+E->n*(p+E->k*j)]*ds[p];
            for(l=0;l<E->k;l++)
             {
              if(E->dddu!=NULL)
                du[i+E->n*(j+E->k*l)]+=E->dddu[i+E->n*(p+E->k*(l+E->k*j))]*ds[p];
             }
           }
         }
       }
     }

   }
  return;
 }

#ifdef __cplusplus
}
#endif
