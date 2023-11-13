
/* 
 *  File   : ksr_define.h
 *  Author : B. Le Cun.
 *  Date   : Nov. 24 2014.
 *  Comment: pthread defines and types
 */


#include<pthread.h>

#define MAXPROCS 40


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

/*--- MEMOIRE PARTAGEE ---*/

#define INIT_SHMEM

#define SH_MALLOC malloc

#define SH_FREE   free



/*--- Variable Partagee ---*/

#define VAR_ALIGN 

#define VAR_SHARED 

#define VAR_PRIVATE 


/*--- Get ID ---*/
#define ARCH_GETID(notask) 


