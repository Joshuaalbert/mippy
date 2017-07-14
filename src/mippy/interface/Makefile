CFLAGS=-I../../../local/include 
#LDFLAGS=-L../../../local/lib
all: 
	gcc --std=c99 lpsolver_interface.c $(CFLAGS) ../../../local/lib/libglpk.a -lm  -o lpsolver_interface
#$(LDFLAGS) 
#-lglpk -lm -o lpsolver_interface
#	cp ../../../local/lib/libglpk.a .
#	gcc lpsolver_interface.o libglpk.a 
