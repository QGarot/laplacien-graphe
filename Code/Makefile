FC=gfortran
OPT=-g

#---------executables------------------
main : meth_puiss.o meth_defl.o
	$(FC) $(OPT) meth_puiss.o meth_defl.o -o main

#---------objects------------------
meth_puiss.o : meth_puiss.f
	$(FC) $(OPT) meth_puiss.f -c
	
meth_defl.o : meth_defl.f
	$(FC) $(OPT) meth_defl.f -c

clean:
	rm -f *.o main
