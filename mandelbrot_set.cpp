/*
 * EDS - Parallel Mandelbrot set generation
 *
 * Author: Krzysztof Voss [shobbo@gmail.com]
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <unistd.h>
#include <ctype.h>

#ifdef TESTED

#include <sys/time.h>
#include <time.h>

#endif

#include "mandelbrot_set.h"
#include "worker.h"
#include "manager.h"
///////////////////////////////////////
char *ofile = NULL;
///////////////////////////////////////

static int
write_ppm(fdata* fd, char* filename)
{
    FILE *fp;
    int x, y;
    char r, g, b;
    char colour;

#ifdef DEBUG
    printf("[Main]->write_ppm\n");
#endif
    if ( filename == NULL )
        fp = fopen("mandelbrot_set.ppm", "wb");
    else
        fp = fopen(filename, "wb");

    /* Inserting the PPM's header */
    fprintf(fp, "P6\n%d %d\n%d\n", fd->resolution, fd->resolution, 255);
    /* Writing down sequentially pixels' values starting with the maximum Y, so with the top */
    for(y=fd->resolution-1; y >= 0; y--) {
        for(x=0; x < fd->resolution; x++) {
            colour = fd->tab[y][x];
            r = g = (9 * colour) % 255;
            /* kolejne kolory teczy. najbardziej rzadkie to najdluzsza fala -> czerwone, najczestsze to krotka fala - fiolet */
            b = (r ^ g) % 255;
            //					r = g = b = (9 * colour) % 255; 
            fprintf(fp, "%c%c%c", r, g, b);
        }
    }
    fprintf(fp, "\n");
    fclose(fp);

    return 0;
}

///////////////////////////////////////

    static int 
clean_table(fdata* fd)
{
    int i;

#ifdef DEBUG
    printf("[Main]->clean_table\n");
#endif

    for(i=0; i < fd->resolution; i++) {
        free(fd->tab[i]);
        fd->tab[i] = NULL;
    }
    free(fd->tab);
    fd->tab = NULL;

    return 0;
}

///////////////////////////////////////

static int
gen_table(fdata* fd)
{
    int i;
    char** tab;

#ifdef DEBUG
    printf("[Main]->gen_table\n");
#endif

    tab = (char**)malloc(fd->resolution * sizeof(char*));
    /* do we need to initialize it with 0s? */
    for(i=0; i < fd->resolution; i++)
        tab[i] = (char*)malloc(fd->resolution * sizeof(char));
    fd->tab = tab;

    return 0;
}

///////////////////////////////////////

static int 
usage(char* appName)
{
    /* displays information on how to run application */
    printf("Usage: %s [options]\n", appName);
    printf("-x\t\tMinimal x value [default: -2.0]\n");
    printf("-X\t\tMaximal x value [default: 2.0]\n");
    printf("-y\t\tMinimal y value [default: -2.0]\n");
    printf("-Y\t\tMaximal y value [default: 2.0]\n");
    printf("-r\t\tResolution of a picture [default: 1024]\n");
    printf("\n");
    printf("-i\t\tMaximal number of iterations [default: 200]\n");
    printf("-t\t\tThreshold [default: 2]\n");
    printf("-n\t\tNumber of simultanously running threads [default: 1 (runs sequentially)]\n");
    printf("\n");
    printf("-m\t\tImplies using MagicBox and POSIX Threads [default: not set]\n");
    printf("-o\t\tImplies using OpenMP (turns off MagicBox) [default: not set]\n");
    printf("-p\t\tImplies using POSIX Threads [default: set]\n");
    printf("-s\t\tSmallest box size (when using MagicBox maximal number of times the rectangle is divided) [default: 4]\n");
    printf("-f\t\tOutput filename [default: mandelbrot_set.ppm]\n");
    printf("-h\t\tPrints this help\n");

    return 0;
}

///////////////////////////////////////

static int
verify(fdata* fd)
{
    if ( fd->xmax <= fd->xmin || fd->ymax <= fd->ymin ) {
        printf("Error: Wrong rectangle parameters were given: ");
        printf("%f %f %f %f\n", fd->xmin, fd->xmax, fd->ymin, fd->ymax);
        return 1;
    }
    if (fd->resolution < 1) {
        printf("Error: Wrong resolution was given\n");
        return 1;
    }
    if (fd->use_mb && fd->sbs <= 0) {
        printf("Error: Wrong smallestBoxSize was given\n");
        return 1;
    }
    if (fd->num_proc < 1) {
        printf("Error: Too few processes were set\n");
        return 1;
    }

    return 0;
}


///////////////////////////////////////

static int
init_fd(fdata* fd, int argc, char* argv[])
{
    int sBox;

#ifdef DEBUG
    printf("[Main]->init_fd\n");
#endif
    int c;

    memset(fd, 0, sizeof(fdata));
    fd->xmin = -2.0;
    fd->xmax = +2.0;
    fd->ymin = -2.0;
    fd->ymax = +2.0;

    fd->resolution = 1024;

    fd->maxiter = 200;
    fd->T = 2;

    fd->num_proc = 1;
    fd->use_mb = 0;
    fd->use_omp = 0;
    sBox = 4;
    fd->sbs = 0;

    opterr = 0;

    while ((c = getopt (argc, argv, "x:X:y:Y:r:i:t:n:f:mophs:")) != -1)
        switch (c) {
            case 'x':
                fd->xmin = atof(optarg);
                break;
            case 'X':
                fd->xmax = atof(optarg);
                break;
            case 'y':
                fd->ymin = atof(optarg);
                break;
            case 'Y':
                fd->ymax = atof(optarg);
                break;
            case 'r':
                fd->resolution = atoi(optarg);
                break;
            case 'i':
                fd->maxiter = atoi(optarg);
                break;
            case 't':
                fd->T = atof(optarg);
                break;
            case 'n':
                fd->num_proc = atoi(optarg);
                break;
            case 'f':
                ofile = optarg;
                break;
            case 'm':
                fd->use_mb = 1;
                fd->use_omp = 0;
                break;
            case 'o':
                fd->use_mb = 0;
                fd->use_omp = 1;
                break;
            case 'p':
                fd->use_omp = 0;
                break;
            case 'h':
                usage(argv[0]);
                return 1;
            case 's':
                fd->use_mb = 1;
                fd->use_omp = 0;
                sBox = atoi(optarg);
                break;

            case '?':
                if ( strchr("xXyYritnfs", optopt) && optopt != 0 )
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n",	optopt);
                return 1;
            default:
                return 1;
        }

    if ( sBox > 0 )
        fd->sbs = (int) pow( (fd->resolution / pow(2,sBox)), 2);

    if ( verify(fd) ) {
        usage(argv[0]);
        return 1;
    }

    fd->wID = -1;
    fd->mutt = NULL;

    fd->xdiff = (fd->xmax - fd->xmin) / fd->resolution;
    fd->ydiff = (fd->ymax - fd->ymin) / fd->resolution;

    return 0;
}

///////////////////////////////////////
#ifdef TESTED
double my_wtime()
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1E-6;
}
#endif

///////////////////////////////////////
///////////////////////////////////////

int 
main (int argc, char* argv[])
{
    fdata* fd;
#ifdef TESTED
    double etime;
    int sManager; // manager exit status
#endif

    printf("=======  EDS - Eve's Distribution System  =======\n");
    printf("Author:  Krzysztof Voss [shobbo@gmail.com]\n\n");

    if ( argc == 1 ) {
        usage(argv[0]);
        return 1;
    }

    fd = (fdata*)malloc(sizeof(fdata));

    if ( init_fd(fd, argc, argv) ) {
        free(fd);
        return 1;
    }
    gen_table(fd);

#ifdef TESTED
    etime = - my_wtime();
    sManager = manager(fd);
    etime += my_wtime();
    printf("Elapsed time: %.3f\n", etime);

    //if ( ! sManager ) write_ppm(fd, ofile);
#else
    if ( ! manager(fd))
        write_ppm(fd, ofile);
#endif

    clean_table(fd);
    free(fd);

    return 0;
}

