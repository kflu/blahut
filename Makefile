blahut:
	gcc -g -c blahut.c main.c
	gcc blahut.o main.o -o blahut -lm `gsl-config --libs`

lib:
	gcc -g -c -fPIC blahut.c
	gcc -shared blahut.o -o blahut.so -lm `gsl-config --libs`

clean:
	rm -f *.o
	rm -f *.so
	rm -f blahut
	rm -f blahut.exe
	rm -f blahut
	rm -f *~
