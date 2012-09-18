/* 
    @(#)sh.h	1.4
    02/04/19 14:43:11
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory
 
        author: Mike Henderson mhender@watson.ibm.com 
*/

#ifndef __SH_H__
#define __SH_H__

#ifdef __cplusplus
 extern "C" {
#endif
void shinit(int*,int*,int*);
void shview(float*,float*,float*,
                   float*,float*,float*,float*,float*,float*);
void shcube(float*,float*,float*,float*,float*,float*);
void shsrfp(float*,float*,float*,int*);
void shpause(void);
void shrefsh(void);
void shclr(void);
void shend(void);

void shpl(int*,float*,float*,float*);
void shdraw(float*,float*,float*,int*);
void shline(float*,float*,float*,float*,float*,float*);
void shlinc(int*,int*,int*);

void shpg(int*,float*,float*,float*);
void shpgcol(int*,float*,float*,float*,int*,int*,int*);
void shpgnrm(int*,float*,float*,float*,float*);
void shpgc(int*,int*,int*);
void shpec(int*,int*,int*);
void shbgc(int*,int*,int*);
void shbec(int*,int*,int*);

void shpnt(float*,float*,float*);
void shpntc(int*,int*,int*);

void shsave(char*,int*,int);
void shstr(float,float,float,char*);
void shstrs(int,int,char*);
void shmask(int*);
void shSetOutputResolution(int,int);
void shSetOutputFilename(char*);
void shSetOutputFormat(char*);

void sharrow(float*,float*,float*,int*,int*,float*,float*,float*,int*,int*);
void shchar(int*,int*,float*,char*,int*,int*,int*,int*);
void shcs(int*,int*,int*);
void shqlinc(int*,int*,int*);
void shline2s(int*,int*,int*,int*,int*,int*);
void shlit(int*,float*,float*,float*,int*,int*,int*);
void shlnonrm(float*,float*,float*,float*,float*,float*,float*,float*,int*,int*,int*);
void shnlit(int*);
void shnpln(int*);
void shvw(float*,float*,float*,float*,float*,float*,float*,float*,float*);
void shqe(float*,float*,float*);
void shsi(float*,float*,float*,float*,float*);
void shpg2(int*,double*,double*,double*);
void shqpgc(int*,int*,int*);
void shqpec(int*,int*,int*);
void shqbgc(int*,int*,int*);
void shqbec(int*,int*,int*);
void shpln(int*,float*,float*,float*,float*,float*,float*,int*);
void shplnrm(int*,float*,float*,float*);
void shqpntc(int*,int*,int*);
void shqlit(int*,float*,float*,float*,int*,int*,int*);
void shqnpln(int*);
void shscal(float*,float*,float*,float*,float*,float*);
void shsphere(float*,float*,float*,float*);
void shsync();
void shtri(float*,float*,float*,float*,float*,float*,float*,float*,float*);
void shtric(int*,int*,int*);

void shSetLineWidth(int);
#ifdef __cplusplus
 }
#endif
#endif
