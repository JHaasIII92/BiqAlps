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
 *  File   : ksr_define.h
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: KSR machine defines and types.
 */


#include<pthread.h>

#define MAXPROCS 40


#ifdef SPG 

/*------- Definition de type du Verrou de l'Arbre -----------*/
typedef char ARCHLock;

/*------- Definition des actions du Verrou de l'Arbre -------*/
#define ARCH_LCKON(p) _gspwt(&(p))
#define ARCH_LCKOFF(p) _rsp(&(p)) 


/*------- Definition de l'initialisation du verrou de l'arbre -*/
#define ARCH_LCKINIT(h) 
/*------- Definition de la suppression du verrou de l'arbre -*/
#define ARCH_LCKDSTR(h) 

#else 
/*--- VERROU ETC... -------*/

/*------- Definition de type du Verrou de l'Arbre -----------*/
typedef pthread_mutex_t ARCHLock;

/*------- Definition des actions du Verrou de l'Arbre -------*/
#define ARCH_LCKON(p) pthread_mutex_lock(&((p).Lock))
#define ARCH_LCKOFF(p) pthread_mutex_unlock(&((p).Lock))

/*------- Definition de l'initialisation du verrou de l'arbre -*/
#define ARCH_LCKINIT(h) pthread_mutex_init(&((h).Lock),NULL)
/*------- Definition de la suppression du verrou de l'arbre -*/
#define ARCH_LCKDSTR(h) pthread_mutex_destroy(&((h).Lock))

#endif


/*--- MEMOIRE PARTAGEE ---*/

#define INIT_SHMEM

#define SH_MALLOC malloc

#define SH_FREE   free



/*--- Variable Partagee ---*/

#define VAR_ALIGN __align128

#define VAR_SHARED __align128 __shared

#define VAR_PRIVATE __private


/*--- Get ID ---*/
#define ARCH_GETID(notask) 


