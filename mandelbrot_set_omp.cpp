/*
 * EDS - Parallel Mandelbrot set generation
 *
 * Author: Krzysztof Voss [shobbo@gmail.com]
 *
 */

/*
 * Wielowatkowe generowanie obrazu zbioru mandelbrota
 * (ujecie wykorzystujace OpenMP)
 *
 * Author: Krzysztof Voss
 * shobbo@gmail.com
 *
 */

#include <cstdio>
#include <complex>
#include <omp.h>

#include "mandelbrot_set.h"
#include "mandelbrot_set_omp.h"

using namespace std;

static int 
fractal_point(complex<double> c, const fdata *fd) //oblicza jak szybko ucieka punkt o wspolrzednych zespolonych
/*
	Z0 = 0
	Zn = Z(n-1) + C

	abs(c) <T
*/
{
	complex<double> Zn;
	int n;

	if ( abs(c) >= fd->T ) //point is already over the range
		return 0;

	for ( n=1, Zn=c;  //zaczynamy od 1 bo przeciez punkt nie jest od razu poza zasiegiem
			abs(Zn)< fd->T && n < fd->maxiter+1;
			 n++)
	{
		Zn = Zn * Zn + c ;
	}

	return n-1; //due to last incrementation in the FOR loop
}

///////////////////////////////////////


int
gen_fractal_omp(const fdata* d)
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

	//void omp_set_num_threads(int num_threads)
	omp_set_num_threads(d->num_proc);
#pragma omp parallel default(shared) private(xl, yl)
	{
		#pragma omp for schedule(dynamic)
		for(yl = 0 ; yl < d->resolution; yl++)
		{
			for(xl=0; xl < d->resolution; xl++)
			{
				c = complex<double>(xmin+(xl+0.5)*xdiff, ymin+(yl+0.5)*ydiff);
				tab[yl][xl] = fractal_point(c, d);
			}
		}
	}

	return 0;
}
