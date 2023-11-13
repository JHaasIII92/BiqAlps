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
 *  File   : sequent_define.h
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: SEQUENT BALANCE machine defines and types.
 */

#include<parallel/parallel.h>
#include<parallel/microtask.h>

#undef MAXPROCS
#define MAXPROCS 9

/*--- VERROU ETC... -------*/

/*------- Definition de type du Verrou de l'Arbre -----------*/
typedef slock_t ARCHLock;

/*------- Definition des actions du Verrou de l'Arbre -------*/
#define ARCH_LOCKON(p)  S_LOCK(&((p).Lock))
#define ARCH_LOCKOFF(p) S_UNLOCK(&((p).Lock))

/*------- Definition de l'initialisation du verrou de l'arbre -*/
#define ARCH_LOCKINIT(h) S_INIT_LOCK(&((h).Lock))
#define ARCH_LOCKDSTR(h) 


/*--- MEMOIRE PARTAGEE ---*/

#define INIT_SHMEM

#define SH_MALLOC shmalloc

#define SH_FREE shfree
 

/*--- Variable Partagee ---*/

#define VAR_ALIGN 

#define VAR_SHARED shared

#define VAR_PRIVATE  


/*--- Get ID ---*/
#define ARCH_GETID(notask) notask=m_get_myid()


