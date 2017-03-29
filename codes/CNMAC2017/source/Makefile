buildupopt: buildupopt.c
	gcc -Wall -Werror buildupopt.c -lgsl -lgslcblas -fopenmp -lm -o buildupopt

run:buildupopt
	./buildupopt 1 1 ../instances/GRAPH_1PLW_N00075.csv
