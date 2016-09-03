CPP=g++
LFLAGS=-pthread -lrt
#CXXFLAGS=-march=core2 -msse4.1 -O2 -pipe
CXXFLAGS=-O2 -pipe -fopenmp


OBJS = mandelbrot_set.o worker.o manager.o mandelbrot_set_omp.o mandelbrot_set_sq.o

mandelbrot_set: $(OBJS)
#		@ echo "Compiling $<..."
		$(CPP) $(CXXFLAGS) $(LFLAGS) $^ -o $@

mandelbrot_set.o: mandelbrot_set.cpp mandelbrot_set.h
mandelbrot_set_sq.o: mandelbrot_set_sq.cpp mandelbrot_set_sq.h
mandelbrot_set_omp.o: mandelbrot_set_omp.cpp mandelbrot_set_omp.h
manager.o: manager.cpp manager.h
worker.o: worker.cpp worker.h

clean:
	-rm -f *.o *.ppm test_procedure-output mandelbrot_set
