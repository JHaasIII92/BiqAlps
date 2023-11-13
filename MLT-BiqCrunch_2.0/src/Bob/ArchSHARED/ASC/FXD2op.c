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
 *  File   : FXDop.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: Global PQ for ARTP=SHARED, MDTP=ASC, DSTP=FXD.
 */

#include <stdio.h>
#include "../../Include/bb.h"
#include "stat.h"

/*----------------------------------------------------------------------*/
VAR_SHARED BobTData NbBlk;

/*------------------------ STAT ----------------------------------------*/
VAR_SHARED BobTDsArchSt dast[MAXPROCS];
VAR_SHARED BobTPsArchSt BobPsArchSt[MAXPROCS];

/*----------------------------------------------------------------------*/
void Bob_PsArchStInit() {
int i,j;

  for (i=0; i<BobCt.NbProcs; i++) {
     BobPsArchSt[i].DelA=0;
     for(j=0;j<MAXPROCS;j++) BobPsArchSt[i].LbDS[j]=0;
     BobPsArchSt[i].LbDS[1]=-1;
  }
}

/*----------------------------------------------------------------------*/
void Bob_PsArchStPrint(io,i)
FILE *io;
int i;
{
  if ( i==-1 ) {
    fprintf(io,"   DelA  Is Dl\n");
  } else {
    fprintf(io," %6d  %2d %2d\n",
           BobPsArchSt[i].DelA,BobPsArchSt[i].InsDS,BobPsArchSt[i].DelDS);
  }
}



/*------------------Load Balancing == 3 ou 4 ---------------------------*/

#if LBTP==3 || LBTP==4 || LBTP==5

typedef struct {
  int Ub,Lb;
} BobInfoBd;

VAR_SHARED BobInfoBd bib[MAXPROCS];

typedef struct { 
   ARCHLock Lock;
   int dsi;
   int NbNode;
   BobTPri LastLB;
} DSI;

VAR_SHARED DSI gldsi;
VAR_PRIVATE int LNb=0;

void DSI_Init()
{
int i;

  for (i=0;i<BobCt.NbProcs;i++) {
     bib[i].Lb = Bob_PRIEVAL(BobRoot->Pri);
     bib[i].Ub = Bob_ULBGet();
  }
  gldsi.dsi = 0;
  gldsi.NbNode = 1;
  Bob_PRIINFI(gldsi.LastLB);
printf("Cou DSI :%d:%d:%d\n",gldsi.dsi,gldsi.NbNode,gldsi.LastLB);
  ARCH_LCKINIT(gldsi);
}

DSI_Test(dsi,BobId)
int dsi,BobId;
{
int ret;
int Nb;

  if ( (BobId%BobCt.NbMast)!=0 ) {
     return BobPsArchSt[BobId].InsDS;
  }
  ARCH_LCKON(gldsi);
#if LBTP==5
  ARCH_LCKON(dst[BobPsArchSt[BobId].InsDS]);
  Nb = dst[BobPsArchSt[BobId].InsDS].NbNode;
  ARCH_LCKOFF(dst[BobPsArchSt[BobId].InsDS]);
  
  ret = (( Nb==0 || ( gldsi.NbNode/2>Nb || 
          Bob_PRIEVAL(gldsi.LastLB)+(bib[BobId].Ub-bib[BobId].Lb)/10 <
              Bob_PRIEVAL(BobPsSt[BobId].LastPri))) ?
                   gldsi.dsi : BobPsArchSt[BobId].InsDS );
#else
  ret = gldsi.dsi;
  if (ret==BobPsArchSt[BobId].InsDS) {
     ret = gldsi.dsi = (gldsi.dsi+1)%BobCt.NbDS;
     gldsi.NbNode =0;
     Bob_PRIINFI(gldsi.LastLB);
  }
#endif
  ARCH_LCKOFF(gldsi);
  return ret;
}

DSI_Set(dsi,BobId,Pri,Nb)
int dsi,BobId,Nb;
BobTPri Pri;
{
  ARCH_LCKON(gldsi);
  if ( gldsi.dsi == BobPsArchSt[BobId].InsDS ) {
     gldsi.NbNode = Nb;
     gldsi.LastLB = Pri;
  } else {
#if LBTP==3 
     if ( gldsi.NbNode<=Nb) {
#elif LBTP==4
     if ( gldsi.LastLB>=Pri ) {
#elif LBTP==5
     if ( gldsi.LastLB>=Pri || gldsi.NbNode<=Nb ) {
#endif
        gldsi.NbNode= Nb;
        gldsi.LastLB= Pri;
        gldsi.dsi = dsi;
     }
  }
  LNb = Nb;
  ARCH_LCKOFF(gldsi);
}
#endif

/*----------------------------------------------------------------------*/

void Bob_DsArchStInit() {
int i;
   for (i=0;i<BobCt.NbDS;i++) {
       dast[i].Capa = 0;
       dast[i].ExpCapa = 0;
       dast[i].NbProcs= 0;
       dast[i].UpBd= 0;
   }
}
void Bob_DsArchStPrint(io,i)
FILE *io;
int i; {
  if ( i==-1 ) {
#if (LBTP==3) || (LBTP==4) || (LBTP==5) 
     fprintf(io," Capa ECapa NbPs UpBd DSI:%d:%d:%d\n", gldsi.dsi,gldsi.NbNode,
             Bob_PRIEVAL(gldsi.LastLB));
#else
     fprintf(io," Capa ECapa NbPs   UpBd\n");
#endif
  } else {
     fprintf(io," %4d %5d %4d %6d\n",dast[i].Capa,dast[i].ExpCapa,
                   dast[i].NbProcs,dast[i].UpBd);
  }
}


/*----------------------------------------------------------------*/
Bob_GPQTestNbNode()
{
  int current_BobId = setget_bobId(0);
  return dst[BobPsArchSt[current_BobId].DelDS].NbNode >= dast[BobPsArchSt[current_BobId].InsDS].NbProcs;
}

/*---------------------------------------------------------*/
void Bob_GPQAlloc()
{
int i,j;
  BobCt.NbDS=(BobCt.NbDS==-1?BobCt.NbProcs:BobCt.NbDS);
  for (i=0;i<BobCt.NbProcs;i++) { 

    BobPsArchSt[i].InsDS=BobPsArchSt[i].DelDS=(i*BobCt.NbDS)/BobCt.NbProcs; 
    dast[(i*BobCt.NbDS)/BobCt.NbProcs].NbProcs++; 
  }
#if LBTP==0
    for (i=0;i<BobCt.NbProcs;i++) { 
      BobPsArchSt[i].LbDS[0]=BobPsArchSt[i].InsDS; 
      BobPsArchSt[i].LbDS[1]=(BobPsArchSt[i].InsDS==0?BobCt.NbDS-1:BobPsArchSt[i].InsDS-1); 
      BobPsArchSt[i].LbDS[2]=(BobPsArchSt[i].InsDS==BobCt.NbDS-1?0:BobPsArchSt[i].InsDS+1); 
      BobPsArchSt[i].LbDS[3]=-1; 
    } 
#endif
#if LBTP==1
    for (i=0;i<BobCt.NbProcs;i++) { 
      BobPsArchSt[i].InsDS = (BobPsArchSt[i].DelDS+1)%BobCt.NbDS;
      BobPsArchSt[i].LbDS[0]=BobPsArchSt[i].DelDS; 
      BobPsArchSt[i].LbDS[1]=(BobPsArchSt[i].DelDS==0?BobCt.NbDS-1:BobPsArchSt[i].DelDS-1); 
      BobPsArchSt[i].LbDS[2]=-1; 
    } 
#endif
#if LBTP==2
    for (i=0;i<BobCt.NbProcs;i++) { 
      j = BobPsArchSt[i].InsDS;
      BobPsArchSt[i].LbDS[0]=j; 
      BobPsArchSt[i].LbDS[1]=(j==0?BobCt.NbDS-1:j/2); 
      BobPsArchSt[i].LbDS[2]=(j*2<BobCt.NbDS?j*2:(j+j/2)%BobCt.NbDS); 
      BobPsArchSt[i].LbDS[3]=(j*2+1<BobCt.NbDS?j*2+1:(j+j/2)%BobCt.NbDS); 
      BobPsArchSt[i].LbDS[4]=-1; 
    } 
#endif
#if LBTP==3 || LBTP==4 || LBTP==5
    for (i=0;i<BobCt.NbProcs;i++) { 
      BobPsArchSt[i].LbDS[0]=BobPsArchSt[i].InsDS; 
      BobPsArchSt[i].LbDS[1]=-1; 
    } 
#endif

  Bob_PQAlloc(BobCt.NbDS,BobRoot->Pri,Bob_ULBGet());

  Bob_DataInit(&NbBlk);
#if LBTP==3 || LBTP==4 || LBTP==5
  DSI_Init(Bob_ULBGet());
#endif
}

void Bob_GPQFree() {
  Bob_PQFree();
}

void Bob_GPQIns(p)
BobNode *p;
{
	int current_BobId = setget_bobId(0);
	int dsi=BobPsArchSt[current_BobId].InsDS;

  BobPsSt[current_BobId].Ins++;
  Bob_DsStIns(dsi);
  Bob_PQIns(dsi,p);
}

Bob_CountBusy()
{
     int i;
     int b = 0;
     for( i = 0 ; i < BobCt.NbProcs;i++) {
         if( busy[i].Val == 1 ) {
             b++;
         }
     }
     return b;
}

Bob_CountSpare()
{
    return BobCt.NbProcs - Bob_CountBusy();
}
 
BobNode *Bob_GPQDel() 
{
  int dsi,i,NbNode;
  int Blck=0;
  BobNode *p;
  int Val;
  int current_BobId = setget_bobId(0);
  p=NULL;
  while(1) {
#if LBTP==5
    dsi = DSI_Test(dsi,current_BobId);
    if ( (p=Bob_PQDel(dsi))!=NULL ) {
      break;
    }
#else
  for(i=0;(dsi=BobPsArchSt[current_BobId].LbDS[i])>-1;i++) {
      if ( (p=Bob_PQDel(dsi))!=NULL ) {
        break;
      }
  }
  if ( p!=NULL ) break;
  #if LBTP==3
    dsi = DSI_Test(dsi,current_BobId);
    if ( (p=Bob_PQDel(dsi))!=NULL ) {
      break;
    }
  #endif
  #if LBTP==4
    dsi = DSI_Test(dsi,current_BobId);
    if ( (p=Bob_PQDel(dsi))!=NULL ) {
       break;
    }
  #endif
#endif
    if ( Blck==0 ) {
       Blck=1;
       ARCH_LCKON(NbBlk);
       NbBlk.Val++;
       ARCH_LCKOFF(NbBlk);
    }
    ARCH_LCKON(NbBlk);
    if ( NbBlk.Val==BobCt.NbProcs ) {
       ARCH_LCKOFF(NbBlk);
       return NULL;
    } else { ARCH_LCKOFF(NbBlk); }
    ARCH_LCKON( busy[current_BobId] );
    busy[current_BobId].Val = 0;
    ARCH_LCKOFF( busy[current_BobId] );
  }
  if ( Blck==1 ) {
     ARCH_LCKON(NbBlk);
     NbBlk.Val--;
     ARCH_LCKOFF(NbBlk);
  }
  if ( dsi!=BobPsArchSt[current_BobId].DelDS ) BobPsArchSt[current_BobId].DelA++;
  BobPsSt[current_BobId].LastPri = p->Pri;
  BobPsSt[current_BobId].Del++;
  ARCH_LCKON(dst[dsi]);
  Val = --(dst[dsi].NbNode);
  dst[dsi].LastPri=p->Pri;
  ARCH_LCKOFF(dst[dsi]);
#if LBTP==3 || LBTP==4 || LBTP==5
  DSI_Set(dsi,current_BobId,p->Pri,Val);
#endif
  ARCH_LCKON( busy[current_BobId] );
  busy[current_BobId].Val = 1;
  ARCH_LCKOFF( busy[current_BobId] );
  return p;
}



void Bob_GPQDelG(bks)
int bks;
{
	int current_BobId = setget_bobId(0);
int dsi = BobPsArchSt[current_BobId].InsDS;
int nb;

  nb = Bob_PQDelG(dsi,bks);
  BobPsSt[current_BobId].Delg+=nb;
  ARCH_LCKON(dst[dsi]);
  dst[dsi].NbNode-=nb;
  ARCH_LCKOFF(dst[dsi]);
}


