CC = gcc
DEBUG = -g
COPT = -c -Wall $(DEBUG) -I. -I..
LINKOPT = -shared
LIBS = -lm `gsl-config --libs`
SHELL = /bin/sh

LibSource = blahut.c
LibObj = blahut.o
LibName = blahut.so
ExampleDir = examples
ExampleSources = $(ExampleDir)/bsc.c \
	         $(ExampleDir)/example1.c \
	         $(ExampleDir)/example2.c
ExamplesObj = $(ExampleSources:.c=.o)
ExamplesTarget = $(ExampleSources:.c=.out)

# gcc -g -c -fPIC blahut.c
# gcc -shared blahut.o -o blahut.so -lm `gsl-config --libs`

%.o: %.c
	$(CC) $(COPT) $< -o $@

lib: $(LibObj)
	$(CC) $(LINKOPT) $? -o $(LibName) $(LIBS)

$(LibObj): $(LibSource)
	$(CC) $(COPT) -fPIC $? -o $(LibObj)

examples: $(ExamplesTarget)

$(ExamplesTarget): $(LibObj) $(ExamplesObj) 
	for i in $(ExamplesObj); \
	    do \
	      Target=$${i/.o/.out}; \
	      $(CC) $$i $(LibObj) -o $$Target $(LIBS); \
	done;

clean:
	rm -f *.o
	rm -f *.so
	rm -f blahut
	rm -f blahut.exe
	rm -f blahut
	rm -f *~
	rm -f $(ExampleDir)/*.o
	rm -f $(ExampleDir)/*.out
	rm -f $(ExampleDir)/*.txt
	rm -f $(ExampleDir)/*~
