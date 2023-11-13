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
 *  File   : define.h
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: shared memory machine types and defines.
 */

typedef int BobTid;

#include "mrk_define.h"

#if defined(__ksr__)
/* #define SPG */
#include "ksr_define.h"

#elif defined(ns32000) 
#include "sequent_define.h"

#else
#include "pthr_define.h"

#endif  

#if PLTP==SEQ
typedef ARCHLock BobPQLock;
#define Bob_PQLINIT(p) 
#define Bob_PQLDSTR(p) 
#define Bob_PQLON(p)  
#define Bob_PQLOFF(p)
#define Bob_EXLON(p)
#define Bob_EXLOFF(p)

#else

typedef ARCHLock BobPQLock;
#define Bob_PQLINIT(p) ARCH_LCKINIT(*(p))
#define Bob_PQLDSTR(p) ARCH_LCKDSTR(*(p))

#if PLTP==EXC
#define Bob_PQLON(p)   
#define Bob_PQLOFF(p) 
#define Bob_EXLON(p)   ARCH_LCKON(*(p))
#define Bob_EXLOFF(p)  ARCH_LCKOFF(*(p))
#else
#define Bob_PQLON(p)   ARCH_LCKON(*(p))
#define Bob_PQLOFF(p)  ARCH_LCKOFF(*(p))
#define Bob_EXLON(p)
#define Bob_EXLOFF(p)
#endif

#endif

#if PLTP==SEQ || PLTP==EXC
typedef char BobNDLock;
#define Bob_NDLINIT(p) 
#define Bob_NDLDSTR(p) 
#define Bob_NDLON(p)  
#define Bob_NDLOFF(p)

#elif PLTP==LCK
typedef ARCHLock BobNDLock;
#define Bob_NDLINIT(p) ARCH_LCKINIT(*(p))
#define Bob_NDLDSTR(p) ARCH_LCKDSTR(*(p))
#define Bob_NDLON(p)   ARCH_LCKON(*(p))
#define Bob_NDLOFF(p)  ARCH_LCKOFF(*(p))
#else
typedef MRKLock BobNDLock;
#define Bob_NDLINIT(p) MRK_LCKINIT(*(p))
#define Bob_NDLDSTR(p) MRK_LCKDSTR(*(p))
#define Bob_NDLON(p)   MRK_LCKON(*(p))
#define Bob_NDLOFF(p)  MRK_LCKOFF(*(p))
#endif

