CC=gfortran
all :   
	$(CC) -c -O3 *.f95 
	$(CC) *.o -O3 -o ale2d



