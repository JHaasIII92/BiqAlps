/*
 *        BOB version 1.0:  Branch and Bound Optimization LiBrary
 *                    PNN Team of PRiSM laboratory
 *             University of Versailles St-Quentin en Yvelines.
 *      Authors:  M. Benaichouche, V. Cung, S. Dowaji, B. Le Cun
 *                      T. Mautor, C. Roucairol.
 *                    (C) 1995 All Rights Reserved
 *
 *                              NOTICE
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose and without fee is hereby granted
 * provided that the above copyright notice appear in all copies and
 * that both the copyright notice and this permission notice appear in
 * supporting documentation.
 *
 * Neither the institutions (Versailles University, PRiSM Laboratory, 
 * the PNN Team), nor the Authors make any representations about the 
 * suitability of this software for any purpose.  This software is 
 * provided ``as is'' without express or implied warranty.
 *
 */

/*
 *  File   : goprocs.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: initializations for ARTP=SHARED.
 */


#include <stdio.h>
#include <signal.h>
#include <sys/syscall.h>
#include "../Include/bb.h"

BobTPsSt  BobPsSt[MAXPROCS];
extern int EVALN;
BusyData busy[MAXPROCS];

/*------------------------------------------------------------------*/
void Bob_PsStInit() {
int i;

  for (i=0;i<BobCt.NbProcs;i++) {
     BobPsSt[i].Evl=0;
     BobPsSt[i].NoExp=0;
     BobPsSt[i].Ins=0;
     BobPsSt[i].Del=0;
     BobPsSt[i].Delg=0;
     Bob_PRINULL(BobPsSt[i].LastPri);
     Bob_PsArchStInit(i);
     busy[i].Val = 0;
  }
}

Bob_STEVL(){ BobPsSt[setget_bobId(0)].Evl++; EVALN++;}
int Bob_EVL(){ return BobPsSt[setget_bobId(0)].Evl;}
Bob_STNEXP(){ BobPsSt[setget_bobId(0)].NoExp++;}
int Bob_NEXP(){ return BobPsSt[setget_bobId(0)].NoExp;}

/*------------------------------------------------------------------*/


void Bob_PsStPrint(io)
FILE *io;
{
int i;
int nbinstot, nbexptot,nbdemtot,nbdelgtot,noinstot;

  nbinstot=0;
  nbdemtot=0;
  nbexptot=0;
  noinstot=0;
  nbdelgtot=0;

  fprintf(io,"--PROCESS ---Stat:----Bob_UpBd:%ld\n",Bob_ULBGet());
  fprintf(io,"          Evl   NoExp     Ins     Dem    Delg       Pri");
  Bob_PsArchStPrint(io,-1);
  for (i=0;i<BobCt.NbProcs;i++) {
    fprintf(io,"%3d    %6d  %6d  %6d  %6d  %6d ",i,
           BobPsSt[i].Evl,BobPsSt[i].NoExp,BobPsSt[i].Ins,BobPsSt[i].Del,
           BobPsSt[i].Delg);Bob_PRIPRINT(io,BobPsSt[i].LastPri);
    Bob_PsArchStPrint(io,i);
    nbinstot += BobPsSt[i].Ins;
    nbdemtot += BobPsSt[i].Del;
    nbdelgtot += BobPsSt[i].Delg;
    nbexptot += BobPsSt[i].Evl;
    noinstot += BobPsSt[i].NoExp;
  }
  fprintf(io,"Tot    %6d  %6d  %6d  %6d  %6d\n",nbexptot,noinstot,nbinstot,nbdemtot,nbdelgtot);
}

/*------------------------------------------------------------------*/
int Bob_NbEvalNode() {
  // Faster : use the BiqCrunch global variable EVALN
  /*
  int i, nbexptot;
  nbexptot=0;
  for (i=0;i<BobCt.NbProcs;i++) {
    nbexptot += BobPsSt[i].Evl;
  }
  return nbexptot;
  */

  return EVALN;
}

/*-------------------------- Data Structure STAT --------------------*/
/*-------------------------------------------------------------------*/
BobTDsSt dst[MAXPROCS];

/*----------------------------------------------------------------------*/

Bob_DsStIns(i) {
    ARCH_LCKON(dst[i]);
    dst[i].NbNode++;
    dst[i].MaxNode = ( dst[i].MaxNode<dst[i].NbNode ? 
                    dst[i].NbNode : dst[i].MaxNode );
    ARCH_LCKOFF(dst[i]);
}

/*----------------------------------------------------------------------*/

Bob_DsStDel(i,Pri)
BobTPri Pri;
{ 
    ARCH_LCKON(dst[i]);
    dst[i].LastPri=Pri;
    dst[i].NbNode--;
    ARCH_LCKOFF(dst[i]);
}
/*----------------------------------------------------------------------*/

void Bob_DsStInit() {
int i;

  for (i=0;i<BobCt.NbDS;i++) {
    ARCH_LCKINIT(dst[i]);
    dst[i].NbNode  = 0;
    dst[i].MaxNode = 0;
    Bob_PRIINFI(dst[i].LastPri);
    Bob_DsArchStInit(i);
  }
#ifdef DOMINANCE
  Bob_ClStInit();
#endif
}

/*----------------------------------------------------------------------*/

Bob_DsStPrint(io)
FILE *io;
{
  int i=0;

  fprintf(io,"-- DSTR %d STAT ----------\n",BobCt.NbDS);
  fprintf(io,"  PQ  NbNode MaxNode   LPri");
  Bob_DsArchStPrint(io,-1);
  for (i=0;i<BobCt.NbDS;i++) {
     fprintf(io," %2d  %6d  %6d ",i,dst[i].NbNode,dst[i].MaxNode);
     Bob_PRIPRINT(io,dst[i].LastPri);
     Bob_DsArchStPrint(io,i);
  }
#ifdef DOMINANCE
  Bob_ClStPrint(io);
#endif

}



/*----------------------------------------------------------------------*/
Bob_Main(n,v,bt) 
int n;
char **v;
BobTTime *bt;
{
  char **Para;
  int NbPara,IndPara;
  
  IndPara = Bob_PrsParam(n,v);
  Para = v+IndPara;
  NbPara = n-IndPara;

  Bob_PsStInit();
  Bob_DsStInit();

  Bob_Init(NbPara,Para);

//  if ( BobSTAT1 ) 
//     printf("Lower Bound :%d  Upper Bound :%d\n",
//             Bob_PRIEVAL(BobRoot->Pri),Bob_ULBGet());

  BobCt.PbSize = BobPbSize;
  Bob_GPQAlloc();
  Bob_CtInit();
  Bob_TimeComp(bt);

  signal(SIGINT,Bob_StPrint);

  Bob_GoJob(Bob_Algo,v[0],&(BobCt.NbProcs),v+1);
  
  Bob_TimeEnd(bt);

  Bob_End();

  Bob_GPQFree();
}

/*----------------------------------------------------------------------*/
int Bob_GoJob(fonct,binary,nbproc,v) 
int *nbproc;
char *binary;
int (*fonct)();
char **v;
{

#if defined(__ksr__)
  pthread_t *pth;
  int i;
  void *return_pointer;
  int Np;
   
  pth = (pthread_t *)malloc((*nbproc)*sizeof(pthread_t));
  if ( pth==NULL ) { 
     fprintf(stderr,"Can't allocate jobs structures\n");
     exit(1);
  }
  for (i=1 ; i<*nbproc ; i++ ) {
     pthread_create(pth+i,NULL,fonct,(void *)i);
  }
  fonct(0);

  for (i=1 ; i<*nbproc ; i++ ) {
     pthread_join(pth[i],&return_pointer);
  }
  free(pth);
  return 1;
  
#elif defined(ns32000)

  while(m_set_procs(*nbproc)==-1) (*nbproc)--;

  m_fork(fonct);
  m_kill_procs();

#else
  pthread_t *pth;
  int i;
  void *return_pointer;
  int Np;
   
  pth = (pthread_t *)malloc((*nbproc)*sizeof(pthread_t));
  if ( pth==NULL ) { 
     fprintf(stderr,"Can't allocate jobs structures\n");
     exit(1);
  }
  for (i=1 ; i<*nbproc ; i++ ) {
     pthread_create(pth+i,NULL,fonct,(void *)i);
  }
  fonct(0);

  for (i=1 ; i<*nbproc ; i++ ) {
     pthread_join(pth[i],&return_pointer);
  }
  free(pth);
  return 1;

#endif
}

