paraloc.exe: paraloc.c va13-1.0.0/src/va13d.f va13-1.0.0/src/ddeps.f main_va13.cpp
	gcc -Wall -Werror paraloc.c -fopenmp -lm  -o paraloc

hello_gsl: hello_gsl.cpp
	gcc hello_gsl.cpp -lgsl -lgslcblas -lm -o hello_gsl
	./hello_gsl

run:paraloc.exe
	./paraloc.exe GRAPH_5CYT.csv
