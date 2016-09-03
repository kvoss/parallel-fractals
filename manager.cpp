/*
 * EDS - Parallel Mandelbrot set generation
 *
 * Author: Krzysztof Voss [shobbo@gmail.com]
 *
 */

#include <unistd.h>
#include <sys/time.h>
#include <errno.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "mandelbrot_set.h"
#include "mandelbrot_set_omp.h"
#include "mandelbrot_set_sq.h"
#include "manager.h"
#include "worker.h"

/*
 * Global data initialization
 *
 */

int freeProc = -1;	/* threads are numbered from 0 on, -1 means invalid or none */
int frees = 0;		/* the number of workers waiting for a job */
int boxing = -1;		/* number of thread currently splitting its box */

pthread_cond_t cond   = PTHREAD_COND_INITIALIZER;	/* worker wakes up manager */
pthread_cond_t mfree  = PTHREAD_COND_INITIALIZER;	/* manager can serve calling thread */
pthread_cond_t mdone  = PTHREAD_COND_INITIALIZER;	/* manager has finished working on current thread */
pthread_cond_t newBox = PTHREAD_COND_INITIALIZER;	/* new box is going to be processed */

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;	/* one needs manager */
pthread_mutex_t muti  = PTHREAD_MUTEX_INITIALIZER;	/* being locked by a worker when it finishes its job */
pthread_mutex_t sharingBox = PTHREAD_MUTEX_INITIALIZER;	/* for sharing a job in MagicBox */


/*
 * Functions' definitions
 *
 */

static int
copy_fd(const fdata* fd1, fdata* fd2)
{
#ifdef DEBUG
	printf("\t[Manager]->copy_fd\n");
#endif
	memcpy(fd2, fd1, sizeof(fdata));

	return 0;
}

////////////////////////////////////////

static int
init_raport(fdata* raport, const fdata* wzor, int numer_procesu, pthread_mutex_t* mutt)
{
	int dy;	/* delta y */

#ifdef DEBUG
	printf("\t[Manager]->init_raport: %d\n", numer_procesu);
#endif

	copy_fd(wzor, raport);

	dy = (int) floor(wzor->resolution/wzor->num_proc);
	raport->yl = 0 + numer_procesu * dy;
	if ( wzor->num_proc == (numer_procesu + 1) )
		raport->yh = wzor->resolution; // because numeration begins with 0
	else
		raport->yh = 0 + (numer_procesu + 1) * dy;

	raport->xl = 0;
	raport->xh = wzor->resolution;

	raport->wID = numer_procesu;
	raport->mutt = mutt;
	raport->status = 1;	// on the beggining we will have job to do

	return 0;
}

////////////////////////////////////////

static int
init_raport_mb(fdata* raport, const fdata* wzor, int numer_procesu, pthread_mutex_t* mutt)
{

#ifdef DEBUG
	printf("\t[Manager]->init_raport_mb: %d\n", numer_procesu);
#endif

	copy_fd(wzor, raport);

	raport->bl = 0 ;
	if ( numer_procesu == 0 )
	{
		raport->bh = 1; // first thread computes the only box existing so far
		raport->status = 1;
	}
	else
	{
		raport->bh = 0; //others have nothing to do
		raport->status = 0;
	}
	raport->wID = numer_procesu;
	raport->mutt = mutt;

	return 0;
}

/*
 * Funkcja odpowiedzialna za przydzielanie pracy
 *
 * Gdy pojawia sie watek z zadaniem, wlasnie tutaj jest ono rozpatrywane
 */
/*
 * Function responsible for assigning work
 * When a thread appears with a request, this is the place where it's considered
 */
static int
manage(pthread_t* threads, fdata** raporty, pthread_mutex_t** mutexy, int num_proc)
{
	int i;
	int nrProc = num_proc;	/* how many threads exist */

#ifdef DEBUG
	printf("[Manager]->manage\n");
#endif

	/*
	 * Pozostajemy w petli dopoki dziala jakis watek
	 * We stay in the loop until any thread is working
	 */
	while( nrProc > 0 )
	{
#ifdef DEBUG
		printf("\t\t[Manager]->Going to sleep..\n");
#endif
		/* zwalnia mutex a zaraz po obudzeniu zajmuje go ponownie */
		/* releases a mutex and just after waking up, aquires and locks it again */
		pthread_cond_wait(&cond, &mutex);

#ifdef DEBUG
		printf("\t\t[Manager]->MANAGER AWAKENED by worker-%d!!!\n", freeProc);
#endif

		/* blokujemy mutex wolnego procesu by miec dostep do jego danych i kontrole nad jego startem */
		/* blocking mutex of the free thread to gain exclusive access to its data and to have a control over its start */
		if ( pthread_mutex_lock(mutexy[freeProc]))
			perror("Manager zaraz po obudzeniu\n");
#ifdef DEBUG
		else
			printf("\t\t[Manager]->Has locked waking thread\n");
#endif


		pthread_mutex_lock(&muti);

		/*
		 * Blokujemy drugi watek
		 * Blocking another thread
		 */
		i = freeProc;
		frees = 1; /* on the beggining we know about only one free thread - the calling one */
		while (1)
		{
#ifdef DEBUG
			printf("\t\t[Manager]->Looking for a busy thread: nrProc=%d, frees=%d\n", nrProc, frees);
#endif
			i = (++i) % num_proc; // looking for another thread

			/* nie chcemy zajac samych siebie a przelecielismy juz wszystkie inne */
			/* we don't want to lock ourself and we have tried all others */
			if ( i == freeProc )
			{
				/* skoro przeszukalismy wszystkie i nie znalezlismy zadnego zajetego to pora konczyc */
				/* niektore moga czekac jeszcze na zmiennej warunkowej lub liczyc ostatnia linie */
				
				/* as we have looked through all threads and have found none busy, it's time to fiinish */
				/* some threads can be waiting on condition or be counting the last line */
#ifdef DEBUG
				printf("\t\t\t[Manager]->Finishing thread #%d\n", i);
#endif
				raporty[freeProc]->status = 2;
				nrProc--;

				pthread_mutex_unlock(mutexy[freeProc]);
				pthread_mutex_unlock(&muti);
				pthread_cond_signal(&mdone);

				freeProc = -1;
				pthread_cond_signal(&mfree);
				break;
			}

			if (! pthread_mutex_lock(mutexy[i]) )
			{
				/*
				 * skoro sie udalo to musimy sprawdzic
				 * czy proces jest naprawde zajety
				 *
				 * if we managed to block it now we have to check
				 * if the thread is really busy
				 */
				if ( raporty[i]->status == 1 )
				{
					if ( raporty[i]->yh == (raporty[i]->yl + 1))
					{
						pthread_mutex_unlock(mutexy[i]);
						continue;
					}
#ifdef DEBUG
					printf("\t\t\t[Manager]->found a busy thread: i=%d [yl: %d, yh: %d]\n", i, raporty[i]->yl, raporty[i]->yh);
#endif
					/*
					 * Mamy zajety watek o numerze i
					 * We have locked thread #i
					 *			 
					 * teraz rozdzielamy prace pomiedzy te watki w nastepujacy sposob:
					 * now we're dividing work between locked threads as follows:
					 *
					 *	job2->yh = job1->yh
					 *	job1->yh = job1->yl + (job1->yh - job1->yl)/2
					 *	job2->yl = job1->yh
					 *
					 *	job2 = raporty[freeProc] //finished
					 *	job1 = raporty[i]	//stopped
					 */
#ifdef DEBUG
					printf("\t\t\t[Manager]->Changing threads' raports.. freeProc=%d and i=%d\n", freeProc, i);
#endif

					raporty[freeProc]->yh = raporty[i]->yh;
					raporty[i]->yh = raporty[i]->yl + (int)floor((raporty[i]->yh - raporty[i]->yl) / 2);
					raporty[freeProc]->yl = raporty[i]->yh;

					if ( raporty[freeProc]->yl < raporty[freeProc]->yh )
						raporty[freeProc]->status = 1;
					else
						raporty[freeProc]->status = 0;

					if ( raporty[i]->yl < raporty[i]->yh )
						raporty[i]->status = 1;
					else
						raporty[i]->status = 0;

#ifdef DEBUG
					printf("\t\t\t[Manager]->Unlocking threads..\n");
#endif
					/* Najpierw startujemy watek, ktory zaczyna od daleszej czesci */
					/* First we are strarting a thread that beggins with latter part */
					pthread_mutex_unlock(mutexy[i]);	
					pthread_mutex_unlock(mutexy[freeProc]);
					pthread_mutex_unlock(&muti);
					pthread_cond_signal(&mdone);

#ifdef DEBUG
					printf("\t\t\t[Manager]->Unlocked mutexes on: %d and %d\n", freeProc, i);
#endif

					freeProc = -1;
					pthread_cond_signal(&mfree);
					break;
				}
				/*
				 * watek zostal zakonczony 
				 * wiec wykonujemy kolejne okrazenie w petli nie wliczajac go do wolnych
				 *
				 * the thread has been finished
				 * so we are doing another lap in the loop, this time not counting it in as a free thread
				 */
				else if ( raporty[i]->status == 0 )
				/* zablokowany watek tez jest wolny wiec nie ma zadnej pracy zeby z nim dzielic */
				/* the blocked thread is also free so there is no work those two can share */
				{
					frees++;
					pthread_mutex_unlock(mutexy[i]);
				}
				else if ( raporty[i]->status == 2 )
				{
					pthread_mutex_unlock(mutexy[i]);
				} //koniec sprawdzania watkow, the end of checking threads
			} // koniec if ktore sprawdza czy udalo sie zajac mutex, the end of the if checking if we managed to aquire the mutex
		} //koniec while do poszukiwania watkow, the end of the while loop for looking for threads
	/*
	 * zrobione wszystko rowna sie z tym, ze zakonczyla sie praca ostatniego watku
	 * having everything done implies working of the last thread is over
	 */
	} /* end of the while loop */

	return 0;
}
///////////////////////////////////////
static int
manage_mb(pthread_t* threads, fdata** raporty, pthread_mutex_t** mutexy, int num_proc)
{
	int timeout = 1;	/* after timeout seconds of waiting for new box we finish any calling thread */
	int i, mstat;		/* iterator and mutexStat */
	int nrProc = num_proc;	/* how many threads exist */
	int busy;		/* how many busy threads we have */
	struct timespec ts;	/* for pthread_cond_timedwait */

#ifdef DEBUG
	printf("[Manager]->manage_mb\n");
#endif

	/*
	 * Obsluga watku - pozostajemy w petli dopoki dziala jakis watek
	 * Thread's service - staying in the loop as long as any thread is working
	 */
	while( nrProc > 0 )
	{
#ifdef DEBUG
		printf("\t[Manager]->Going to sleep..\n");
#endif
		/* zwalnia mutex a zaraz po obudzeniu zajmuje go ponownie */
		/* releases a mutex and just after waking up, aquires and locks it again */
		pthread_cond_wait(&cond, &mutex);
#ifdef DEBUG
		printf("\t\t[Manager]->MANAGER AWAKENED by worker-%d!!!\n", freeProc);
#endif

		/* czekamy az mutex muti sie zwolni co bedzie znaczylo, ze watek
		 * jest gotowy otrzymac komunikat mdone
		 *
		 * waiting for muti mutex being released what means that the calling thread having it before
		 * is ready for getting mdone message
		 */
		pthread_mutex_lock(&muti);

		/* blokujemy mutex wolnego watku by miec dostep do jego danych i kontrole nad jego startem */
		/* blocking mutex of the free thread to gain exclusive access to its data and to have a control over its start */
		if ( pthread_mutex_lock(mutexy[freeProc]))
			perror("Manager cannot lock a thread after waking up\n");
#ifdef DEBUG
		else
			printf("\t\t[Manager]->Has locked waking thread\n");
#endif
		
		/* waiting for a new box being executed */
#ifdef DEBUG
		printf("\t\t[Manager]->Waiting for a thread splitting box\n");
#endif
		pthread_mutex_lock(&sharingBox);
		clock_gettime(CLOCK_REALTIME, &ts);
		ts.tv_sec += timeout;
		if ( (mstat = pthread_cond_timedwait(&newBox, &sharingBox, &ts)) == ETIMEDOUT )
		{
			pthread_mutex_unlock(&sharingBox);
			/* we were waiting over 2 seconds for new box */
#ifdef DEBUG
			printf("\t\t\t[Manager]->Finishing thread #%d[TIMEDOUT]\n", freeProc);
#endif
			raporty[freeProc]->status = 2;
			nrProc--;

			pthread_mutex_unlock(mutexy[freeProc]);
			pthread_mutex_unlock(&muti);
			pthread_cond_signal(&mdone);

			freeProc = -1;
			pthread_cond_signal(&mfree);
			continue;
		}
		else //if ( mstat == 0 )
		{
			i = boxing;
			//boxing = -1;
			pthread_mutex_unlock(&sharingBox);
#ifdef DEBUG
			printf("\t\t\t[Manager]->Discovered some new boxes from thread #%d\n", i);
#endif
			if( i >= 0 )
			{
				pthread_mutex_lock(mutexy[i]);
#ifdef DEBUG
				printf("\t\t\t\t[Manager]->Having thread #%d locked\n", i);
#endif
			}
			else
			{
				/* zaden worker jeszcze nie podzielil */
				/* none worker has split any box */
#ifdef DEBUG
				printf("[Manager]->BOXING < 0 !!!\n");
#endif
				raporty[freeProc]->status = 0;

				pthread_mutex_unlock(mutexy[freeProc]);
				pthread_mutex_unlock(&muti);
				pthread_cond_signal(&mdone);

				freeProc = -1;
				pthread_cond_signal(&mfree);
				continue;
			}

			/*
			 * skoro sie udalo to musimy sprawdzic
			 * czy proces jest naprawde zajety
			 *
			 * now let's check if the thread is really busy at the moment
			 */
			if ( raporty[i]->status == 1 )
			{
				if ( raporty[i]->bh == (raporty[i]->bl + 1))
				{
					pthread_mutex_unlock(mutexy[i]);

					raporty[freeProc]->status = 0;

					pthread_mutex_unlock(mutexy[freeProc]);
					pthread_mutex_unlock(&muti);
					pthread_cond_signal(&mdone);

					freeProc = -1;
					pthread_cond_signal(&mfree);
					continue;
				}
#ifdef DEBUG
				printf("\t\t\t[Manager]->Thread #%d is splitting a box, [bl: %d, bh: %d]\n", i, raporty[i]->bl, raporty[i]->bh);
#endif
				/*
				 * Mamy zajety watek o numerze i
				 * Thread #i is locked
				 *			 
				 * teraz rozdzielamy prace pomiedzy te watki w nastepujacy sposob:
				 * now we share the work between those two threads as follows:
				 *
				 *	job2->yh = job1->yh
				 *	job1->yh = job1->yl + (job1->yh - job1->yl)/2
				 *	job2->yl = job1->yh
				 *
				 *	job2 = raporty[freeProc] //finished
				 *	job1 = raporty[i]	//stopped
				 */
#ifdef DEBUG
				printf("\t\t\t[Manager]->Changing threads' raports.. freeProc=%d and i=%d\n", freeProc, i);
#endif
				raporty[freeProc]->bh = raporty[i]->bh;
				raporty[i]->bh = raporty[i]->bl + (int)floor((raporty[i]->bh - raporty[i]->bl) / 2);
				raporty[freeProc]->bl = raporty[i]->bh;

				if ( raporty[freeProc]->bl < raporty[freeProc]->bh )
					raporty[freeProc]->status = 1;
				else
					raporty[freeProc]->status = 0;

				if ( raporty[i]->bl < raporty[i]->bh )
					raporty[i]->status = 1;
				else
					raporty[i]->status = 0;

#ifdef DEBUG
				printf("\t\t\t[Manager]->Unlocking threads..\n");
#endif
				/* Najpierw startujemy watek, ktory zaczyna od daleszej czesci */
				/* First we start a thread that owns later part */
				pthread_mutex_unlock(mutexy[i]);	
				pthread_mutex_unlock(mutexy[freeProc]);
				pthread_mutex_unlock(&muti);
				pthread_cond_signal(&mdone);

#ifdef DEBUG
				printf("\t\t\t[Manager]->Unlocked mutexes on: %d and %d\n", freeProc, i);
#endif

				freeProc = -1;
				pthread_cond_signal(&mfree);
				continue;
			}
			else
			{
#ifdef DEBUG
				printf("\t\t\t[Manager]->Status of locked thread #%d: %d\n", i, raporty[i]->status);
#endif			
				pthread_mutex_unlock(mutexy[i]);
				pthread_mutex_unlock(mutexy[freeProc]);
				pthread_mutex_unlock(&muti);
				pthread_cond_signal(&mdone);
				freeProc = -1;
				pthread_cond_signal(&mfree);
				continue;
			}		
		}
	} /* end of the while loop */

	return 0;
}
///////////////////////////////////////
int
manager(const fdata* wzor)
{

	int i;
	fdata* raport;
	fdata** raporty;		/* list of  raports assigned to workers */

	pthread_mutex_t* mutt;
	pthread_mutex_t** mutexy;	/* list of mutexes assigned to workers */

	pthread_t* threads;		/* list of thread_ids */

#ifdef DEBUG
	printf("[Manager]->manager\n");
#endif

	/* Using sequential algorithm */
	if ( wzor->num_proc == 1 )
	{
#ifdef DEBUG
	printf("[Manager]->gen_fractal_sq\n");
#endif
		gen_fractal_sq(wzor);
		return 0;
	}

#ifdef DEBUG
	printf("[Manager]->Number of threads: %d\n", wzor->num_proc);
#endif
	/* Using OpenMP */
	if ( wzor->use_omp )
	{
#ifdef DEBUG
	printf("[Manager]->gen_fractal_omp\n");
#endif
		gen_fractal_omp(wzor);
		return 0;
	}

	/*
	 * Hiring personel
	 */
	threads = (pthread_t*) malloc(wzor->num_proc * sizeof(pthread_t));

		raporty = (fdata**) malloc(wzor->num_proc * sizeof(fdata*));
		mutexy = (pthread_mutex_t**) malloc(wzor->num_proc * sizeof(pthread_mutex_t*));

	/*
	 * Zakladamy mutex poniewaz tworzone procesy zaczynaja upominac sie o przydzialy pracy
	 * Locking the mutex because being created threads start to ask for a work
	 *
	 * Nastepnie tworzymy kolejno wymagane struktury danych, inicjalizujemy je a nastepnie
	 * startujemy watki.
	 * Next we create requiered data structures, initilize them and start the threads
	 */
	pthread_mutex_lock(&mutex);
		for(i = wzor->num_proc - 1; i >= 0; i--)
		{

#ifdef DEBUG
			printf("[Manager]->Creating thread: %d\n", i);
#endif

			/* creation and initialisation of the worker's mutex with default values */
			mutexy[i] = mutt = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t));
			pthread_mutex_init(mutt, NULL);

			raporty[i] = raport = (fdata*) malloc(sizeof(fdata));

			/* init raport */
			if(wzor->use_mb)
				init_raport_mb(raport, wzor, i, mutt);
			else
				init_raport(raport, wzor, i, mutt);

			/* create thread */
			if( pthread_create(&threads[i], NULL, worker, (void*)raport) )
			{
				printf("\t[Manager]->ERROR; return code from pthread_create() is not 0\n");
			}
		}


		/*
		 * Tworzenie watkow zostalo ukonczone, przekazujemy sterowanie do funkcji, ktora nimi zarzadza
		 * The threads have been created, we pass the control to the function that manages them
		 */
		if(wzor->use_mb)
			manage_mb(threads, raporty, mutexy, wzor->num_proc);
		else
			manage(threads, raporty, mutexy, wzor->num_proc); //pamietajmy o mutexie

	/* po skonczonym zarzadzaniu i zamknieciu innych watkow mozemy zwolnic mutex */
	/* after finishing managing and having other threads joined, we can release the mutex*/
	pthread_mutex_unlock(&mutex);

	/*
	 * Oczekiwanie na zakonczenie kolejnych watkow
	 * Waiting for other threads being finished
	 */
	 for(i=0; i < wzor->num_proc; i++)
	 {
#ifdef DEBUG
	 	printf("[Manager]->Joining thread #%d\n", i);
#endif
	 	pthread_join(threads[i], NULL);
	 }


	/*
	 * Zwalnianie zasobow
	 * Releasing resources
	 */
	for(i=0; i < wzor->num_proc; i++)
	{
		pthread_mutex_destroy(mutexy[i]);
		free(mutexy[i]);
		free(raporty[i]);
	}

	free(mutexy);
	free(raporty);
	free(threads);

	return 0;
}

