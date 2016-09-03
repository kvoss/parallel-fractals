/*
 * EDS - Parallel Mandelbrot set generation
 *
 * Author: Krzysztof Voss [shobbo@gmail.com]
 *
 */

#include <cstdio>
#include <pthread.h>
#include <complex>
#include <cmath>
//#include <sched.h>

#include "mandelbrot_set.h"
#include "worker.h"

using namespace std;

/*
 *	Z0 = 0
 *	Zn = Z(n-1) + C
 *
 *	abs(c) <T
 */
static int 
fractal_point(const complex<double> c, fdata *fd)
/* oblicza jak szybko ucieka punkt o wspolrzednych zespolonych */
/* counts how fast the point described with complex coordinates is moving from its origins */
{
	complex<double> Zn;
	int n;

	if ( abs(c) >= fd->T ) //point is already over the range
		return 0;

	for ( n=1, Zn=c;  // we start with 1 because point is not over the range on the very beginning
			abs(Zn)< fd->T && n < fd->maxiter+1;
			 n++)
	{
		Zn = Zn * Zn + c;
	}

	return n-1; //due to the last incrementation in the FOR loop
}

///////////////////////////////////////

static void
get_job(fdata* fd)
{

#ifdef DEBUG
	printf("[Worker-%d]->Locking manager\n", fd->wID);
#endif
	if (! pthread_mutex_lock(&mutex))
	{
#ifdef DEBUG
		// we have locked the manager
		printf("\t[Worker-%d]->Manager locked!!!\n", fd->wID);
#endif
		// setting the number of the free thread
		while(1)
		{
			if ( freeProc < 0 )
			{
#ifdef DEBUG
				printf("\t[Worker-%d]->Manager is ready to serve us!!!\n", fd->wID);
#endif
				freeProc = fd->wID;
				break;
			}
			// a jezeli ktos juz ustawil swoj numer to czekamy az menadzer go zwolni
			// and if somebody has already set its number, we wait until manager changes it back
			else
			{
#ifdef DEBUG
				printf("\t[Worker-%d]->Sleeping until peer has been served\n", fd->wID);
#endif
				pthread_cond_wait(&mfree, &mutex);
				continue;
			}
		}

		/* zdejmuje mutex z szefa przy spelnionym warunku */
		/* unlocks the manage's mutex when condition is satisfied*/
		pthread_mutex_unlock(&mutex);

		pthread_mutex_lock(&muti);

			/* budzimy szefa */
			/* waking manager up */
			pthread_cond_signal(&cond);
#ifdef DEBUG
			printf("\t\t[Worker-%d]->Waking manager up\n", fd->wID);
#endif

			pthread_cond_wait(&mdone, &muti);
#ifdef DEBUG
			//printf("\t\t\t[Worker-%d]->Woke up after pthread_cond_wait\n", fd->wID);
#endif
		pthread_mutex_unlock(&muti);
	}
}

///////////////////////////////////////

static int
gen_fractal(fdata* d)
{
	int yl;
	int xl;
	complex<double> c;

	double ydiff = d->ydiff;
	double xdiff = d->xdiff;

	char** tab = d->tab;

	double xmin = d->xmin;
	double xmax = d->xmax; 
	double ymin = d->ymin;
	double ymax = d->ymax; 

#ifdef DEBUG
	printf("[Worker-%d]->gen_fractal: %d %d\n", d->wID, d->yl, d->yh);
#endif

	pthread_mutex_lock(d->mutt);
	while(1)
	{
		/* sprawdzamy czy mamy co liczyc */
		/* do we have anything to count? */
		yl = d->yl;
		if ( yl < d->yh )
		{
			pthread_mutex_unlock(d->mutt);
		}
		else
		{
			// nie mamy co robic
			// we don't have anything to do
			d->status = 0;
			break;
		}

		for(xl=d->xl; xl < d->xh; xl++)
		{
			/* dodajemy jeszcze pol roznicy bo chcemy liczyc od srodka pierwszego pola */
			/* adding half of the field size because we want to count beggining with the middle of the first field */
			c = complex<double>(xmin+(xl+0.5)*xdiff, ymin+(yl+0.5)*ydiff); 
			tab[yl][xl] = fractal_point(c, d);
		}

		
		// po zrobieniu linii zamykamy klodke, jezeli jest zamknieta to znaczy, ze szef zatrzymuje tutaj watek.
		// after completing the line we lock the mutex, if its already locked it means manager has stopped the thread over here
		// blad: MOZLIWE PRZEOCZENIE LINII - ale nie zmieniamy d->yl wiec spokojnie
		pthread_mutex_lock(d->mutt);

			++yl;
			d->yl = yl;
			if ( !(yl < d->yh) )
			{
				d->status = 0;
				break;
			}
	}
	pthread_mutex_unlock(d->mutt);

	return 0;
}

//////////////////////////////////////
// MagicBox Implementation
//////////////////////////////////////

static void
map(int type, cords* c)
{
	switch (type)
	{
		case 0:
			c->xmin = (c->xmin + c->xmax) / 2;
			c->ymax = (c->ymin + c->ymax) / 2;
			break;
		case 1:
			c->xmax = (c->xmin + c->xmax) / 2;
			c->ymax = (c->ymin + c->ymax) / 2;
			break;
		case 2:
			c->xmax = (c->xmin + c->xmax) / 2;
			c->ymin = (c->ymin + c->ymax) / 2;
			break;
		case 3:
			c->xmin = (c->xmin + c->xmax) / 2;
			c->ymin = (c->ymin + c->ymax) / 2;
			break;
		default:
			printf("Unknown switch command\n");
	}

	return;
}
//////////////////////////////////////
static int 
addr(fdata* fd, int sq, cords* c)
{
	int osq; //other square

	if ( sq == 0 )
	{
		c->xmin = fd->xmin;
		c->xmax = fd->xmax;
		c->ymin = fd->ymin;
		c->ymax = fd->ymax;
		return 0;
	}
	osq = (int)floor(( sq - 1 ) / 4);
	/* RECURRENCE */
	addr(fd, osq, c);
	map( sq%4, c);
}
//////////////////////////////////////
static int
addrPix(fdata* fd, cords* c)
/* we are checking what pixels that would be and we are saving them in our raport */
{
	fd->xl = floor((c->xmin - fd->xmin) / fd->xdiff);
	fd->xh = floor((c->xmax - fd->xmin) / fd->xdiff);
	fd->yl = floor((c->ymin - fd->ymin) / fd->ydiff);
	fd->yh = floor((c->ymax - fd->ymin) / fd->ydiff);

	return 0;
}
//////////////////////////////////////
//////////////////////////////////////
static int
sizeBox(fdata* fd)
{
	int x, y;
	x = fd->xh - fd->xl;
	y = fd->yh - fd->yl;
	return x*y;
}
//////////////////////////////////////
static void
countBox(fdata *fd)
{
	int xl, yl;

	complex<double> c;

	double ydiff = fd->ydiff;
	double xdiff = fd->xdiff;

	char** tab = fd->tab;

	double xmin = fd->xmin;
	double xmax = fd->xmax; 
	double ymin = fd->ymin;
	double ymax = fd->ymax; 

#ifdef DEBUG
//	printf("\t\t[Worker-%d]->countBox\n", fd->wID);
#endif
	for ( yl = fd->yl; yl < fd->yh; yl++ )
		for ( xl = fd->xl; xl < fd->xh; xl++ )
		{
			c = complex<double>(xmin+(xl+0.5)*xdiff, ymin+(yl+0.5)*ydiff); 
			tab[yl][xl] = fractal_point(c, fd);
		}

	return;
}
//////////////////////////////////////
static void
fulfillBox(fdata* fd, int v)
/* fulfills whole box b with value v */
{
	int xl, yl;
	char** tab = fd->tab;
#ifdef DEBUG
//	printf("\t\t[Worker-%d]->fulfillBox v=%d, xl=%d, xh=%d, yl=%d, yh=%d\n", fd->wID, v, fd->xl, fd->xh, fd->yl, fd->yh);
#endif

	for ( yl = fd->yl; yl < fd->yh; yl++ )
		for ( xl = fd->xl; xl < fd->xh; xl++ )
			tab[yl][xl] = v;

	return;
}
//////////////////////////////////////
static int gen_fractal_mb(fdata*);
//////////////////////////////////////
static void
splitBox(fdata* fd, int b)
/* splits box b into 4 smaller boxes and starts processing them */
{
	int bl, bh;

#ifdef DEBUG
//	printf("\t\t[Worker-%d]->splitBox, box:%d\n", fd->wID, b);
#endif
	/* first we save bl and bh on stack and we assign them new values */
	pthread_mutex_lock(fd->mutt);
		bl = fd->bl;
		bh = fd->bh;
		fd->bl = (b * 4) + 1;
		fd->bh = (b * 4) + 5; // +1+4
	pthread_mutex_unlock(fd->mutt);

	/* informing manager about an opportunity to share a job */
	/* tutaj mozemy spotkac sie z tym, ze jeden watek wejdzie do
	 * sharingBox, wpisze swoja wartosc a po nim zapisze wartosc inny
	   ale to w niczym nie przeszkadza bo w rezultacie i tak mamy 
	   numer zajetego watku.
	 */
	pthread_mutex_lock(&sharingBox);
		boxing = fd->wID;
	pthread_mutex_unlock(&sharingBox);
	/* I prefer to call it outside mutex, first less instructions inside
	 * second, after receving signal manager won't wait for mutex
	 */
	pthread_cond_signal(&newBox);

	/* then we call a fucntion that will count new boxes */
	/* RECURRENCE */
	gen_fractal_mb(fd);

	/* after completing smaller box we get back to counting another big one */
	pthread_mutex_lock(fd->mutt);
		fd->bl = bl;
		fd->bh = bh;
	pthread_mutex_unlock(fd->mutt);

	return;	
}
//////////////////////////////////////
static int
processBox(fdata* fd, int b)
/* Processes box b */
{
	complex<double> c0, c1, c2, c3;
	int yl, yh, xl, xh, i, range;
	int p0, p1, p2, p3, p;
	double xmin, ymin, xdiff, ydiff;

#ifdef DEBUG
	//printf("\t[Worker-%d]->processBox\n", fd->wID);
#endif
	/* we need exact coordinates of the box we are about to check */
	cords c;
	addr(fd, b, &c);

	/* we are checking what pixels that would be and we are saving them in our raport */
	addrPix(fd, &c);

	/* finally we are checking if every value on the border is equal
	 * if not we are splitting the box into 4 smaller ones
	 */
	xl = fd->xl, xh = fd->xh;
	yl = fd->yl, yh = fd->yh;
	xmin = fd->xmin, ymin = fd->ymin;
	xdiff = fd->xdiff, ydiff = fd->ydiff;
	range = xh - xl;
	for(i=0; i < range; i++ )
	{
		c0 = complex<double>(xmin+(xl+0.5)*xdiff, ymin+(yl+i+0.5)*ydiff);
		c1 = complex<double>(xmin+(xl+i+0.5)*xdiff, ymin+(yh-0.5)*ydiff);
		c2 = complex<double>(xmin+(xh-0.5)*xdiff, ymin+(yl+i+0.5)*ydiff);
		c3 = complex<double>(xmin+(xl+i+0.5)*xdiff, ymin+(yl+0.5)*ydiff);

		p0 = fractal_point(c0, fd);
		p1 = fractal_point(c1, fd);
		p2 = fractal_point(c2, fd);
		p3 = fractal_point(c3, fd);

		/* initialization of the variable keeping past value */
		if ( i == 0 )
			p = p0;

		/* if values are different */
		if ( (p ^ p0) || (p ^ p1) || (p ^ p2) || (p ^ p3) )
		{
			/* we have to check if box is big enought to consider splitting it, otherwise we count it normally */
			if ( sizeBox(fd) < fd->sbs )
			{
				countBox(fd);
				return 0;
			}
			splitBox(fd, b);
			return 0;
		}
	}
#ifdef DEBUG
	//printf("\t[Worker-%d]->processBox(before fulfillBox, xl=%d, xh=%d, yl=%d, yh=%d\n", fd->wID, fd->xl, fd->xh, fd->yl, fd->yh);
#endif
	fulfillBox(fd, p);

	return 0;
}

//////////////////////////////////////
static int
gen_fractal_mb(fdata* d)
{
	int bl;

#ifdef DEBUG
	printf("[Worker-%d]->gen_fractal_mb: %d %d\n", d->wID, d->bl, d->bh);
#endif

	pthread_mutex_lock(d->mutt);
	while(1)
	{
		/* do we have anything to count? */
		bl = d->bl;

		if ( bl < d->bh )
		{
			d->status = 1;
			pthread_mutex_unlock(d->mutt);
		}
		else
		{
			d->status = 0;
			break;
		}

		/* There is at least one box to process */
		processBox(d, bl);

		/* after finishing the box we lock our mutex, if it's already closed it means manager has locked it */
		pthread_mutex_lock(d->mutt);

			++bl;
			d->bl = bl;
			if ( !(bl < d->bh) )
			{
				d->status = 0;
				break;
			}
	}
	pthread_mutex_unlock(d->mutt);

	return 0;
}

//////////////////////////////////////

void*
worker(void* d) //d jak dane ;]
{
	fdata* fd = (fdata*) d;
	int lstat; // local status

	while(1)
	{
		pthread_mutex_lock(fd->mutt);
		lstat = fd->status; // bc of the weird statuses on exit I made a local copy

#ifdef DEBUG
		printf("[Worker-%d]->started; status: %d\n", fd->wID, lstat);
#endif

		/* sprawdzamy czy mamy co robic albo konczymy
		 *
		 * jezeli przydzielono nam pusta prace to prosimy o kolejny przydzial
		 * przeciez nie mozemy sami sie zwolnic z pracy bo nie mamy wiedzy na temat tego
		 * czy nie bedziemy potrzebni gdzie indziej
		 */
		/* we check if we have anything to do or we finish
		 *
		 * if we have null work assiged we ask for another assignment
		 * we cannot fire ourselves from job because we don't know if we are not needed
		 * somewhere else in a moment
		 */
		if ( lstat == 1 )
		{
			pthread_mutex_unlock(fd->mutt);		
		}
		else if ( lstat == 0 )
		{
			pthread_mutex_unlock(fd->mutt);
			get_job(fd);
			continue;
		}
		else if ( lstat == 2 )
		{
			pthread_mutex_unlock(fd->mutt);
			break;
		}
		else
		{
			printf("\t[Worker-%d]->Unknown status value: %d. Finishing.\n", fd->wID, lstat);
			pthread_mutex_unlock(fd->mutt);
			break;
		}

		/* 
		 * Wykonujemy prace, ktora nam zadano; moze ona zostac przerwana
		 * po skonczeniu ustawia status w raporcie na 0
		 *
		 * We do the work we were given; it can be interrupted
		 * after finishing it we have our status set to 0
		 */
		if (fd->use_mb)
			gen_fractal_mb(fd);
		else
			gen_fractal(fd);

		get_job(fd);
	}
#ifdef DEBUG
	printf("[Worker-%d]->RIP !!!\n", fd->wID);
#endif

	return 0;
}
