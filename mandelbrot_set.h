/*
 * EDS - Parallel Mandelbrot set generation
 *
 * Author: Krzysztof Voss [shobbo@gmail.com]
 *
 */

#ifndef MANSETH
#define MANSETH

#include <pthread.h>	

#define DEBUG
#ifdef TESTED
#undef DEBUG
#endif

/*
 * SHARED VARIABLES
 */

extern int freeProc;	/* threads are numbered from 0 on, -1 means invalid or none */
extern int frees;	/* the number of workers waiting for a job */
extern int boxing;	/* number of thread currently splitting its box */

extern pthread_cond_t cond;	/* budzenie managera; worker wakes up manager */
extern pthread_cond_t mfree;	/* manager moze przyjac nowy watek; manager can serve calling thread */
extern pthread_cond_t mdone;	/* manager skonczyl zmieniac dane watku; manager has finished working on current thread */
extern pthread_cond_t newBox;	/* watek zaczyna nowy box; new box is going to be processed */

extern pthread_mutex_t mutex;	/* one needs manager */
extern pthread_mutex_t muti;	/* being locked by a worker when it finishes its job */
extern pthread_mutex_t sharingBox;	/* for sharing a job in MagicBox */


/*
 * MandelbrotSet struct
 */
typedef struct {
    double ydiff, xdiff;	/* distances between adjoining pixels */
    char** tab;		/* table with results */

    /* image's parameters */
    double xmin, xmax;	/* x range */
    double ymin, ymax; 	/* y range */
    int resolution;	/* resolution of the picture */

    int maxiter;		/* maximal number of iterations */
    double T;		/* threshold */

    int num_proc;		/* number of threads */
    int use_mb, use_omp;	/* whether to use MagicBox or not */
    int sbs;		/* smallest box size for MagicBox (in square pixels) */

    /* Workers' individual data */
    int yl, yh, xl, xh;	/* assigned work */
    int bl, bh;		/* assigned work - BoxLow BoxHigh */
    int wID;		/* worker's ID */
    pthread_mutex_t* mutt;	/* thread's mutex */
    int status;		/* thread's status (0 - free, 1 - busy, 2 - released) */

} fdata;

typedef struct {
    double xmin, xmax;
    double ymin, ymax;
} cords;

#endif

