paraloc.exe: paraloc.c
	gcc -Wall -Werror paraloc.c -lgsl -lgslcblas -fopenmp -lm -o paraloc

hello_gsl: hello_gsl.cpp
	gcc hello_gsl.cpp -lgsl -lgslcblas -lm -o hello_gsl
	./hello_gsl

run:paraloc.exe
	./paraloc.exe GRAPH_5CYT.csv
