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
 *  File   : mrk_define.h
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: Booelan-mark macros.
 */

/*-------------------------------------------------------------------

---------------------------------------------------------------------*/


#define MARKED 1
#define UNMARKED 0

/*------- Definition de type du Verrou d'un noeud -----------*/
typedef char MRKLock;

/*------- Definition des actions du verrou d'un noeud -------*/
#define MRK_LCKON(p) { while ( (p).Lock==MARKED); \
                     (p).Lock=MARKED;}
#define MRK_LCKOFF(p) (p).Lock=UNMARKED

/*------- Definition de l'initialisation du verrou de l'arbre -*/
#define MRK_LCKINIT(p) (p).Lock = UNMARKED
#define MRK_LCKDSTR(p) 

