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
 *  File   : upbd.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: Global Upper/Lower Bound for ARTP=SHARED,MDTP=ASC.
 */


#include <stdio.h>
#include "../../Include/bb.h"
#include "stat.h"

BobTULB BobUpBd;


/*----------------------------------------------------------------------*/
int Bob_ULBGet() {
  return BobUpBd.Val;
}

/*----------------------------------------------------------------------*/
void Bob_ULBInit(ub,bs)
int ub;
BobSolution *bs;
{
int i;

  BobSol = (BobSolution *)SH_MALLOC(sizeof(BobSolution));
  if (BobSol==NULL ) {
     fprintf(stderr,"Not enough memory for BobSolution\n");
     exit(1);
  }
  *BobSol = *bs;
  Bob_DataInit(&BobUpBd);
  BobUpBd.Val = ub;
#if DSTP==FXD || DSTP==MRG
  for(i=0;i<BobCt.NbDS; i++ ) 
    dast[i].UpBd = ub;
#endif

}

/*----------------------------------------------------------------------*/
int Bob_ULBUpd(NewBks,bs) 
int NewBks;
BobSolution *bs;
{
int ret=0;
int current_BobId = setget_bobId(0);

#if DSTP==ONE
int i=0;
#else
int i = BobPsArchSt[current_BobId].InsDS;
#endif

  ARCH_LCKON(BobUpBd);
  if ( Bob_EVALL(BobUpBd.Val,NewBks) ) {
     BobUpBd.Val = NewBks;
     *BobSol = *bs;
     ARCH_LCKOFF(BobUpBd); 
#if DSTP!=ONE
     if ( BobCt.NbDS>1 ) {
       ARCH_LCKON(dst[i]);
       dast[i].UpBd = NewBks;
       ARCH_LCKOFF(dst[i]);
     }
#endif
     Bob_PQDelG(i,NewBks);
     ret= 1;
  } else {
     ARCH_LCKOFF(BobUpBd); 
  }
  return ret;
}

/*----------------------------------------------------------------------*/
int Bob_ULBSup(Pri) 
int Pri;
{
int ret;
int ub;
int current_BobId = setget_bobId(0);

#if DSTP==ONE
int i=0;
#else
int i = BobPsArchSt[current_BobId].InsDS;
#endif

  ARCH_LCKON(BobUpBd);
  ret= (Bob_EVALL(BobUpBd.Val,Pri));
  ub = BobUpBd.Val;
  ARCH_LCKOFF(BobUpBd); 
#if DSTP!=ONE
  if ( BobCt.NbDS>1 ) {
     ARCH_LCKON(dst[i]);
     if ( Bob_EVALL(dast[i].UpBd,ub) ) {
         dast[i].UpBd = ub;
         ARCH_LCKOFF(dst[i]);
         Bob_GPQDelG(ub);
     } else 
         ARCH_LCKOFF(dst[i]);
  }
#endif
  return ret;
}

