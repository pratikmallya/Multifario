/* @(#)CmpExpr.c	1.2
   02/05/03 09:47:19 
 
    PROGRAM NAME:  multifario
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

    author: mhender@watson.ibm.com               */

#include <ExpCmpI.h>
 
#ifdef __cplusplus
extern "C" {
#endif
int ECCompileExpression(char *sourceCode,struct ECObjectCode **object)
 {
  struct ECTokenizedCode *tokenized;
  int rc;
  int i;
  int Destination;
 
  tokenized=ECTokenizeExpression(sourceCode);
  if((*tokenized).successful!=0)
   {
    if(ECPrintErrorMessages)fprintf(stderr,"ECCompileExpression: tokenization failed\n");
    rc=(*tokenized).successful;
    ECFreeTokenizedCode(tokenized);
    return(rc);
   }

  *object=ECParseTokenizedCode(tokenized);
  if((*object)->successful!=0)
   {
    if(ECPrintErrorMessages)fprintf(stderr,"ECCompileExpression: parsing failed\n");
    rc=(*object)->successful;
    ECFreeTokenizedCode(tokenized);
    ECFreeObjectCode(object);
    return(rc);
   }
  ECFreeTokenizedCode(tokenized);

  for(i=0;i<(*object)->nStatements;i++)
   {
    Destination=((*object)->ExecutableCode[i]).Destination;
    if(Destination>(*object)->nRegisters)(*object)->nRegisters=Destination+1;
   }

  return((*object)->successful);
 }
#ifdef __cplusplus
 }
#endif
