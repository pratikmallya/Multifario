/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   June 20, 2005
 */


static char *id="@(#) $Id: MFErrorHandler.c,v 1.3 2007/02/13 01:22:33 mhender Exp $";

#include <MFErrorHandler.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFErrorHandlerSt {
                  int nErrors;
                  int mErrors;
                  int *severity;
                  char **routine;
                  char **message;
                  char **filename;
                  int *line;

                  int nRefs;
                 };

MFErrorHandler MFCreateErrorHandler()
 {
  static char RoutineName[]="MFCreateErrorHandler";
  MFErrorHandler thisErrHand;

  thisErrHand=(MFErrorHandler)malloc(sizeof(struct MFErrorHandlerSt));
  if(thisErrHand==(MFErrorHandler)NULL)
   {
    fprintf(stderr,"malloc of an MFErrorHandlerSt failed. Length is %d bytes\n",sizeof(struct MFErrorHandlerSt));
    fflush(stderr);
    abort();
   }

  thisErrHand->nErrors=0;
  thisErrHand->mErrors=0;
  thisErrHand->severity=(int*)NULL;
  thisErrHand->routine=(char**)NULL;
  thisErrHand->message=(char**)NULL;
  thisErrHand->filename=(char**)NULL;
  thisErrHand->line=(int*)NULL;
  thisErrHand->nRefs=1;

  return thisErrHand;
 }

void MFRefError(MFErrorHandler e)
 {
  static char RoutineName[]="MFRefError";

  e->nRefs++;
  return;
 }

void MFFreeErrorHandler(MFErrorHandler e)
 {
  static char RoutineName[]="MFFreeError";
  int i;

  e->nRefs--;
  if(e->nRefs>0)return;

  if(e->routine!=(char**)NULL)
   {
    for(i=0;i<e->nErrors;i++)
     if(e->routine[i]!=(char*)NULL)free(e->routine[i]);
   }
  if(e->message!=(char**)NULL)
   {
    for(i=0;i<e->nErrors;i++)
     if(e->message[i]!=(char*)NULL)free(e->message[i]);
   }
  if(e->filename!=(char**)NULL)
   {
    for(i=0;i<e->nErrors;i++)
     if(e->filename[i]!=(char*)NULL)free(e->filename[i]);
   }
  if(e->severity!=(int*)NULL)free(e->severity);
  if(e->line!=(int*)NULL)free(e->line);

  free(e);
  return;
 }

void MFSetError(MFErrorHandler e, int sev,char *routine,char *msg,int line,char *file)
 {
  static char RoutineName[]="MFSetError";

  if(!(e->nErrors<e->mErrors))
   {
    e->mErrors++;
    e->severity=(int*)realloc((void*)e->severity,e->mErrors*sizeof(int));
    if(e->severity==(int*)NULL)
     {
      printf("Catastrophic error in MFSetError, out of memory trying to alloc %d bytes for e->severity, line %d in file %s\n",
             e->mErrors*sizeof(int),__LINE__,__FILE__);
      printf("  Error message in stack is MFSetError(%d,\"%s\",\"%s\",%d,\"%s\");\n",sev,routine,msg,line,file);
      fflush(stdout);
      abort();
     }
    e->routine=(char**)realloc((void*)e->routine,e->mErrors*sizeof(char*));
    if(e->routine==(char**)NULL)
     {
      printf("Catastrophic error in MFSetError, out of memory trying to alloc %d bytes for e->routine, line %d in file %s\n",
             e->mErrors*sizeof(int),__LINE__,__FILE__);
      printf("  Error message in stack is MFSetError(%d,\"%s\",\"%s\",%d,\"%s\");\n",sev,routine,msg,line,file);
      fflush(stdout);
      abort();
     }
    e->message=(char**)realloc((void*)e->message,e->mErrors*sizeof(char*));
    if(e->message==(char**)NULL)
     {
      printf("Catastrophic error in MFSetError, out of memory trying to alloc %d bytes for e->message, line %d in file %s\n",
             e->mErrors*sizeof(int),__LINE__,__FILE__);
      printf("  Error message in stack is MFSetError(%d,\"%s\",\"%s\",%d,\"%s\");\n",sev,routine,msg,line,file);
      fflush(stdout);
      abort();
     }
    e->filename=(char**)realloc((void*)e->filename,e->mErrors*sizeof(char*));
    if(e->filename==(char**)NULL)
     {
      printf("Catastrophic error in MFSetError, out of memory trying to alloc %d bytes for e->filename, line %d in file %s\n",
             e->mErrors*sizeof(int),__LINE__,__FILE__);
      printf("  Error message in stack is MFSetError(%d,\"%s\",\"%s\",%d,\"%s\");\n",sev,routine,msg,line,file);
      fflush(stdout);
      abort();
     }
    e->line=(int*)realloc((void*)e->line,e->mErrors*sizeof(int));
    if(e->line==(int*)NULL)
     {
      printf("Catastrophic error in MFSetError, out of memory trying to alloc %d bytes for line, line %d in file %s\n",
             e->mErrors*sizeof(int),__LINE__,__FILE__);
      printf("  Error message in stack is MFSetError(%d,\"%s\",\"%s\",%d,\"%s\");\n",sev,routine,msg,line,file);
      fflush(stdout);
      abort();
     }
   }
  e->severity[e->nErrors]=sev;
  e->routine[e->nErrors]=(char*)malloc((strlen(routine)+1)*sizeof(char));
  strcpy(e->routine[e->nErrors],routine);
  e->message[e->nErrors]=(char*)malloc((strlen(msg)+1)*sizeof(char));
  strcpy(e->message[e->nErrors],msg);
  e->filename[e->nErrors]=(char*)malloc((strlen(file)+1)*sizeof(char));
  strcpy(e->filename[e->nErrors],file);
  e->line[e->nErrors]=line;

  if(0)
   {
    printf("  severity[%d]=%d\n",e->nErrors,e->severity[e->nErrors]);fflush(stdout);
    printf("  routine[%d]=\"%s\"\n",e->nErrors,e->routine[e->nErrors]);fflush(stdout);
    printf("  message[%d]=\"%s\"\n",e->nErrors,e->message[e->nErrors]);fflush(stdout);
    printf("  filename[%d]=\"%s\"\n",e->nErrors,e->filename[e->nErrors]);fflush(stdout);
    printf("  line[%d]=%d\n",e->nErrors,e->line[e->nErrors]);fflush(stdout);
    fflush(stdout);
    abort();
   }

  if(sev>0)
   {
    fprintf(stderr,"  %s(%s:%d): %s\n",e->routine[e->nErrors],e->filename[e->nErrors],e->line[e->nErrors],e->message[e->nErrors]);
    fflush(stderr);
   }
  e->nErrors++;

  return;
 }

int MFGetNErrors(MFErrorHandler e)
 {
  static char RoutineName[]="MFGetNErrors";

  return e->nErrors;
 }

int MFErrorHandlerGetSev(MFErrorHandler e, int n)
 {
  static char RoutineName[]="MFErrorHandlerGetSev";

  if(n<0||!(n<e->nErrors))return 0;
  return e->severity[n];
 }

char *MFErrorHandlerGetRoutine(MFErrorHandler e, int n)
 {
  static char RoutineName[]="MFErrorHandlerGetRoutine";

  if(n<0||!(n<e->nErrors))return (char*)NULL;
  return e->routine[n];
 }

char *MFErrorHandlerGetMsg(MFErrorHandler e, int n)
 {
  static char RoutineName[]="MFErrorHandlerGetMsg";

  if(n<0||!(n<e->nErrors))return (char*)NULL;
  return e->message[n];
 }

int MFErrorHandlerGetLine(MFErrorHandler e, int n)
 {
  static char RoutineName[]="MFErrorHandlerGetLine";

  if(n<0||!(n<e->nErrors))return 0;
  return e->line[n];
 }

char *MFErrorHandlerGetFile(MFErrorHandler e, int n)
 {
  static char RoutineName[]="MFErrorHandlerGetFile";

  if(n<0||!(n<e->nErrors))return (char*)NULL;
  return e->filename[n];
 }

int MFIsError(MFErrorHandler e)
 {
  static char RoutineName[]="MFIsError";

  return(e->nErrors!=0);
 }

void MFClearErrors(MFErrorHandler e)
 {
  static char RoutineName[]="MFClearErrors";

  int i;
  for(i=0;i<e->nErrors;i++)
   {
    if(e->routine[i]!=(char*)NULL)free(e->routine[i]);
    if(e->message[i]!=(char*)NULL)free(e->message[i]);
    if(e->filename[i]!=(char*)NULL)free(e->filename[i]);
   }
  if(e->severity!=(int*)NULL)free(e->severity);e->severity=(int*)NULL;
  if(e->routine!=(char**)NULL)free(e->routine);e->routine=(char**)NULL;
  if(e->message!=(char**)NULL)free(e->message);e->message=(char**)NULL;
  if(e->filename!=(char**)NULL)free(e->filename);e->filename=(char**)NULL;
  if(e->line!=(int*)NULL)free(e->line);e->line=(int*)NULL;
  e->mErrors=0;
  e->nErrors=0;

  return;
 }

void MFErrorHandlerOutOfMemory(MFErrorHandler e)
 {
  static char RoutineName[]="MFErrorHandlerOutOfMemory";

  exit(12);
 }

#ifdef __cplusplus
}
#endif
