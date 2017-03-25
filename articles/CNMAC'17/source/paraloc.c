#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <omp.h>
#include <math.h>

typedef struct {
	int nnodes;
	int nedges;
	int *i;
	int *p;
	double *d;
	double *e;
	int *degree;
	int max_degree;
} Graph;

typedef struct {
	int key;
	double value;
} KeyValue;

static int keyvalue_compare(const void *a, const void *b) {
	KeyValue *A = (KeyValue*)a;
	KeyValue *B = (KeyValue*)b;
	if (A->value > B->value) return 1;
	if (A->value < B->value) return -1;
	return 0;
}

void print_int(char *name, int *d, int n) {
	int i;
	if (n == 1) {
		printf("%s = %d\n", name, d[0]);
	}
	else {
		printf("%s = [ ", name);
		for (i = 0; i < n; i++) {
			printf(" %d ", d[i]);
		}
		printf("]\n");
	}
}

void print_double(char *name, double *d, int n) {
	int i;
	if (n == 1) {
		printf("%s = %.8g\n", name, d[0]);
	}
	else {
		printf("%s = [ ", name);
		for (i = 0; i < n; i++) {
			printf(" %.8g ", d[i]);
		}
		printf("]\n");
	}
}


double vec_dot(double *u, double *v) {
	return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

double vec_norm(double *u) {
	return sqrt(vec_dot(u, u));
}

void vec_cross(double *u, double *v, double *x) {
	x[0] = u[1] * v[2] - u[2] * v[1];
	x[1] = u[2] * v[0] - u[0] * v[2];
	x[2] = u[0] * v[1] - u[1] * v[0];
	return;
}

void vec_diff(double *x, double *y, double *z) {
	z[0] = x[0] - y[0];
	z[1] = x[1] - y[1];
	z[2] = x[2] - y[2];
	return;
}

void vec_normalize(double *x) {
	double nrm = vec_norm(x);
	x[0] /= nrm;
	x[1] /= nrm;
	x[2] /= nrm;
	return;
}

void vec_copy(double *x, double *y) {
	y[0] = x[0];
	y[1] = x[1];
	y[2] = x[2];
	return;
}

double vec_dist(double *x, double *y) {
	double dx, dy, dz;
	dx = x[0] - y[0];
	dy = x[1] - y[1];
	dz = x[2] - y[2];
	double dist = sqrt(dx * dx + dy * dy + dz * dz);
	/*print_double("x", x, 3);
	print_double("y", y, 3);
	printf("|x-y| = %g\n", dist);*/
	return dist;
}

void vec_print(char *name, double *x) {
	int i, n = 3;
	printf("%s = [ ", name);
	for (i = 0; i < n; i++) {
		printf("%g ", x[i]);
	}
	printf("]\n");
}

void graph_read_edges(char *filename, int **ei, int **ej, double **ed, int *nnodes, int *nedges) {
	char * line = NULL;
	size_t len = 0, nlines = 0;
	int read;
	FILE *fid = fopen(filename, "r");
	double dij;
	int i, j;

	if (fid == NULL) {
		printf("The file %s could not be opened.\n", filename);
		exit(EXIT_FAILURE);
	}

	// memory allocation (counting the file number of lines)
	while ((read = getline(&line, &len, fid)) != -1) {
		nlines++;
	}
	*nedges = nlines - 1;

	*ei = (int*)malloc(sizeof(int) * (*nedges));
	*ej = (int*)malloc(sizeof(int) * (*nedges));
	*ed = (double*)malloc(sizeof(double) * (*nedges));
	if (*ei == NULL || *ej == NULL || *ed == NULL) {
		printf("Memory could not be allocated\n");
	}

	printf("\nReading graph file: %s\n", filename);
	rewind(fid);
	// skipping file header
	read = getline(&line, &len, fid);
	nlines = 0;
	while ((read = getline(&line, &len, fid)) != -1) {
		if (sscanf(line, "%d,%d,%lf\n", &i, &j, &dij) != 3) {
			printf("   Line[%zu] %s could not be read.\n", nlines + 1, line);
			exit(EXIT_FAILURE);
		}

		// update edges
		(*ei)[nlines] = i - 1; // convert to zero-based
		(*ej)[nlines] = j - 1; // convert to zero-based
		(*ed)[nlines] = dij;
		nlines++;

		// update nnodes
		if (i < j) i = j;
		if ((*nnodes) < i) (*nnodes) = i;
	}

	fclose(fid);
}

void graph_read(char *filename, Graph *G) {
	int i, j, k, nedges, nnodes, *degree;
	int *ei, *ej;
	double *ed;

	graph_read_edges(filename, &ei, &ej, &ed, &nnodes, &nedges);

	// set degrees
	degree = (int*)calloc(nnodes, sizeof(int));
	for (k = 0; k < nedges; k++) {
		degree[ei[k]]++;
		degree[ej[k]]++;
	}

	// compressing edges
	int    *gp = (int*)malloc((nnodes + 1) * sizeof(int));
	int    *gi = (int*)malloc(2 * nedges * sizeof(int));
	double *gd = (double*)malloc(2 * nedges * sizeof(double));
	double *ge = (double*)malloc(2 * nedges * sizeof(double));
	int    *gn = (int*)calloc(nnodes, sizeof(int));

	if (gp == NULL || gi == NULL || gd == NULL) {
		printf("Memory could not be allocated.\n");
		exit(EXIT_FAILURE);
	}

	int max_degree = 0;
	gp[0] = 0;
	for (k = 0; k < nnodes; k++) {
		gp[k + 1] = gp[k] + degree[k];
		if (max_degree < degree[k]) {
			max_degree = degree[k];
		}
	}

	for (k = 0; k < nedges; k++) {
		i = ei[k];
		j = ej[k];
		gi[gp[i] + gn[i]] = j;
		gd[gp[i] + gn[i]] = ed[k];
		gn[i]++;

		gi[gp[j] + gn[j]] = i;
		gd[gp[j] + gn[j]] = ed[k];
		gn[j]++;
	}

	// set graph
	G->i = gi;
	G->d = gd;
	G->p = gp;
	G->e = ge;
	G->nnodes = nnodes;
	G->nedges = nedges;
	G->degree = degree;
	G->max_degree = max_degree;

	printf("       #Edges : %d\n", nedges);
	printf("       #Nodes : %d\n", nnodes);
	printf("   max(degree): %d\n", max_degree);

	free(ei);
	free(ej);
	free(ed);
	free(gn);
}

void graph_get_neighs(Graph *G, int i, int **neighs, int *n) {
	*neighs = &(G->i[G->p[i]]);
	*n = G->p[i + 1] - G->p[i];
}

int graph_edge(Graph *G, int i, int j, double *dij, double *eij) {
	int k;
	for (k = G->p[i]; k < G->p[i + 1]; k++) {
		if (G->i[k] == j) {
			if (dij != NULL) *dij = G->d[k];
			if (eij != NULL) *eij = G->e[k];
			return 1;
		}
	}
	//printf("G(%d,%d): %g\n", i, j, dij);
	return 0;
}

double graph_dij(Graph *G, int i, int j) {
	double dij;
	graph_edge(G, i, j, &dij, NULL);
	return dij;
}

double graph_eij(Graph *G, int i, int j) {
	double eij;
	graph_edge(G, i, j, NULL, &eij);
	return eij;
}

void graph_update_error_edge(Graph *G, int i, int j, double eij) {
	int k;
	
	//printf("Update[G(%d,%d)]: eij = %g\n", i, j, eij);

	// loc edge ij
	for (k = G->p[i]; k < G->p[i + 1]; k++) {
		if (G->i[k] == j) {
			G->e[k] = eij;
			break;
		}
	}

	// swap i,j to reuse the same code to update edge ji
	k = i; i = j; j = k;
	for (k = G->p[i]; k < G->p[i + 1]; k++) {
		if (G->i[k] == j) {
			G->e[k] = eij;
			break;
		}
	}

	//print_double("e", G->e, 2 * G->nedges);
}

void graph_update_error_neighs(Graph *G, int k, double *x, int *isloc) {
	int i, j, n, *neighs;
	double eij, dij;
	graph_get_neighs(G, k, &neighs, &n);
	//printf("N[%d] = ", k); print_int("", neighs, n);
	for (i = 0; i < n; i++)
	{
		j = neighs[i];
		if (isloc[j]) {
			dij = graph_dij(G, k, j);
			eij = fabs(vec_dist(&x[3 * k], &x[3 * j]) - dij) / dij;
			graph_update_error_edge(G, k, j, eij);
		}
	}
}

double graph_error_node(Graph *G, int k, int *isloc) {
	double sum_eij = 0, eij;
	int *neighs, n, i, j;
	if (isloc[k]) {
		graph_get_neighs(G, k, &neighs, &n);
		for (i = 0; i < n; i++) {
			j = neighs[i];
			if (isloc[j]) {
				graph_edge(G, k, j, NULL, &eij);
				sum_eij += eij;
			}
		}
		return sum_eij;
	}
	else {
		return 0;
	}
}

void graph_print_error_node(Graph *G, int k, int *isloc) {
	if (isloc[k]) {
		printf("e[%d] = %g\n", k, graph_error_node(G, k, isloc));
	}
	else {
		printf("e[%d] = 0 (isloc[%d] == 0)\n", k, k);
	}
}

void graph_print_error_stats(Graph *G, int *isloc) {
	double max_e = 0, mean_e = 0, eij;
	int i;
	for (i = 0; i < G->nnodes; i++) {
		eij = graph_error_node(G, i, isloc);
		if (eij > max_e) max_e = eij;
		mean_e += eij;
		if (G->nnodes < 30) printf("e[%d] = %g\n", i, eij);
	}
	mean_e /= G->nnodes;
	printf("   Mean(RelError) = %e\n", mean_e);
	printf("   Max (RelError) = %e\n", max_e);
}

double mat_det(double *A) {
	return A[0] * (A[4] * A[8] - A[5] * A[7])
		- A[1] * (A[3] * A[8] - A[5] * A[6])
		+ A[2] * (A[3] * A[7] - A[4] * A[6]);
}

void mat_print(char *name, double *A, int m, int n) {
	int i, j;
	printf("%s =\n", name);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			printf("%3.2g   ", A[i * n + j]);
		}
		printf("\n");
	}
}

void mat_solve(double *A, double *b, double *x) {
	double detA;
	double detC;
	double c[3];
	int i;

	detA = mat_det(A);
	//printf("detA = %g\n", detA);
	for (i = 0; i < 3; i++) {
		// save the i-th columns
		c[0] = A[0 + i];
		c[1] = A[3 + i];
		c[2] = A[6 + i];

		// replace the i-th column by b
		A[0 + i] = b[0];
		A[3 + i] = b[1];
		A[6 + i] = b[2];

		detC = mat_det(A);
		x[i] = detC / detA;

		//if (~VERBOSE) printf("detC[%d] = %g\n", i, detC);

		// restore the i-th column
		A[0 + i] = c[0];
		A[3 + i] = c[1];
		A[6 + i] = c[2];
	}
	/*print_double("   A", A, 9);
	print_double("   b", b, 3);
	printf("   detA= %g\n", detA);
	print_double("x", x, 3);*/
}


void intersect3Spheres(double *a, double ra, double *b, double rb, double *c, double rc, double *x_pos, double *x_neg) {
	double p[3], y[3], A[9], *u, *w, *n;
	double a2, b2, c2, ra2, rb2, rc2, dpa, rp;

	if (x_pos == NULL && x_neg == NULL) {
		return;
	}

	u = A + 0;
	w = A + 3;
	n = A + 6;

	// calculates the normal
	// u = b - a
	// w = c - a
	// n = normalize(u x w)
	vec_diff(b, a, u);
	vec_diff(c, a, w);
	vec_cross(u, w, n);
	vec_normalize(n);

	// set y = rhs
	a2 = vec_dot(a, a);
	b2 = vec_dot(b, b);
	c2 = vec_dot(c, c);
	ra2 = ra * ra;
	rb2 = rb * rb;
	rc2 = rc * rc;
	y[0] = (ra2 - rb2 + b2 - a2) / 2.0;
	y[1] = (ra2 - rc2 + c2 - a2) / 2.0;
	y[2] = vec_dot(n, a);

	// Ap = y
	mat_solve(A, y, p);

	dpa = vec_dist(p, a);
	rp = sqrt(ra2 - dpa*dpa);

	// x_pos = p + rp * n
	if (x_pos != NULL) {
		x_pos[0] = p[0] + rp * n[0];
		x_pos[1] = p[1] + rp * n[1];
		x_pos[2] = p[2] + rp * n[2];
	}

	// x_neg = p - rp * n
	if (x_neg != NULL) {
		x_neg[0] = p[0] - rp * n[0];
		x_neg[1] = p[1] - rp * n[1];
		x_neg[2] = p[2] - rp * n[2];
	}
}

void intersect4Spheres(double *a, double ra, double *b, double rb, double *c, double rc, double *d, double rd, double *x) {
	double y[3], A[9], *u, *v, *w;
	double a2, b2, c2, d2, ra2, rb2, rc2, rd2;

	/*print_double("   a", a, 3);
	print_double("   b", b, 3);
	print_double("   c", c, 3);
	print_double("   d", d, 3);
	printf("   ra = %.8g\n", ra);
	printf("   rb = %.8g\n", rb);
	printf("   rc = %.8g\n", rc);
	printf("   rd = %.8g\n", rd)*/;

	// set pointers
	u = &A[0];
	v = &A[3];
	w = &A[6];

	a2 = vec_dot(a, a);
	b2 = vec_dot(b, b);
	c2 = vec_dot(c, c);
	d2 = vec_dot(d, d);
	ra2 = ra * ra;
	rb2 = rb * rb;
	rc2 = rc * rc;
	rd2 = rd * rd;

	// set matrix coeffs
	vec_diff(b, a, u);
	vec_diff(c, a, v);
	vec_diff(d, a, w);

	// set rhs 
	y[0] = (ra2 - rb2 + b2 - a2) / 2.0;
	y[1] = (ra2 - rc2 + c2 - a2) / 2.0;
	y[2] = (ra2 - rd2 + d2 - a2) / 2.0;

	mat_solve(A, y, x);

	//printf("error: %g %g %g %g\n",
	//	fabs(vec_dist(x, a) - ra) / ra,
	//	fabs(vec_dist(x, b) - rb) / rb,
	//	fabs(vec_dist(x, c) - rc) / rc,
	//	fabs(vec_dist(x, d) - rd) / rd);
}


int base_isclique(Graph *G, int *B, int n) {
	int i, j, k, Bi, Bk, count, *neighs, nneighs;

	for (i = 0; i < n; i++) {
		Bi = B[i];
		graph_get_neighs(G, Bi, &neighs, &nneighs);
		//printf("[%d]: ", B[i]);
		//neighs_print(neighs, nneighs);
		count = 0;
		for (j = 0; j < nneighs; j++) {
			for (k = 0; k < n; k++) {
				Bk = B[(i + k) % n];
				//printf("Bk: %d\n", Bk);
				if (neighs[j] == Bk) {
					count++;
				}
			}
		}
		if (count < (n - 1)) return 0;
	}
	return 1;
}

int base_anchor(Graph *G, int k, int *isloc, int *b) {
	int *neighs, i, j, n, bn = 0;
	KeyValue *kv = (KeyValue*)malloc(sizeof(KeyValue) * G->max_degree);

	// get neighs
	graph_get_neighs(G, k, &neighs, &n);

	// set kv
	for (i = 0; i < n; i++) {
		j = neighs[i];
		if (isloc[j]) {
			kv[bn].key = j;
			kv[bn].value = graph_error_node(G, j, isloc);
			// discount the error of edge (k,j)
			/*if (isloc[k]) {
				double eij = graph_eij(G, k, j);
				kv[bn].value -= eij;
			}*/
			bn++;
		}
	}

	if (bn > 3) {
		// select the neighs with smallest error
		qsort(kv, bn, sizeof(KeyValue), keyvalue_compare);
		for (i = 0; i < 4; i++) {
			b[i] = kv[i].key;
		}
		free(kv);
		return 1;
	}
	else {
		return 0;
	}
}

void base_update_x(Graph *G, int *b, int k, double *x) {
	double *A = &x[3 * b[0]];
	double *B = &x[3 * b[1]];
	double *C = &x[3 * b[2]];
	double *D = &x[3 * b[3]];
	double ra = graph_dij(G, k, b[0]);
	double rb = graph_dij(G, k, b[1]);
	double rc = graph_dij(G, k, b[2]);
	double rd = graph_dij(G, k, b[3]);
	intersect4Spheres(A, ra, B, rb, C, rc, D, rd, &x[3 * k]);
}

void base_init(Graph *G, int *B, double *x) {
	int i;

	//set x[B[0]]
	i = 3 * B[0];
	x[i + 0] = 0.0;
	x[i + 1] = 0.0;
	x[i + 2] = 0.0;
	// vec_print("x1", x + i);

	//set x[B[1]]
	double d01 = graph_dij(G, B[0], B[1]);
	i = 3 * B[1];
	x[i + 0] = d01;
	x[i + 1] = 0.0;
	x[i + 2] = 0.0;
	// vec_print("x2", x + i);

	//set x[B[2]]
	double d02 = graph_dij(G, B[0], B[2]);
	double d12 = graph_dij(G, B[1], B[2]);
	double cos_theta = (d02*d02 + d01*d01 - d12*d12) / (2 * d01*d02);
	double sin_theta = sqrt(1 - cos_theta * cos_theta);
	i = 3 * B[2];
	x[i + 0] = d02 * cos_theta;
	x[i + 1] = d02 * sin_theta;
	x[i + 2] = 0.0;
	// vec_print("x3", x + i);

	//set x[B[3]]
	double d03 = graph_dij(G, B[0], B[3]);
	double d13 = graph_dij(G, B[1], B[3]);
	double d23 = graph_dij(G, B[2], B[3]);
	intersect3Spheres(x + 3 * B[0], d03, x + 3 * B[1], d13, x + 3 * B[2], d23, x + 3 * B[3], NULL);
	// vec_print("x4", x + 3 * B[3]);
}

void base_locx(Graph *G, int *B, double *x) {
	int i, j, k, sn, sn_new, nneighs, nfixed = 0;
	int *s, *s_new, *swap, *neighs, b[4];
	int *scores, *isloc;
	int n = 4;

	// memory allocation
	isloc = (int*)malloc(sizeof(int) * G->nnodes);
	s = (int*)malloc(sizeof(int) * G->nnodes);
	s_new = (int*)malloc(sizeof(int) * G->nnodes);
	scores = (int*)malloc(sizeof(int) * G->nnodes);

	// check allocation
	if (s == NULL || s_new == NULL || isloc == NULL || scores == NULL) {
		printf("Memory could not be allocated.\n");
		exit(EXIT_FAILURE);
	}

	// set isloc = scores = zeros
	for (i = 0; i < G->nnodes; i++) {
		isloc[i] = 0;
		scores[i] = 0;
	}

	// init isloc
	for (i = 0; i < n; i++) isloc[B[i]] = 1;

	// init scores
	for (i = 0; i < n; i++) scores[B[i]] = n;

	// init s
	sn = n;
	for (i = 0; i < n; i++) s[i] = B[i];

	// init erros
	for (i = 0; i < n; i++) graph_update_error_neighs(G, B[i], x, isloc);

	do {
		nfixed += sn;
		//if (VERBOSE) print_int("s", s, sn);
		sn_new = 0;
		// set location and update the list of nodes to be fixed
		for (i = 0; i < sn; i++) {
			graph_get_neighs(G, s[i], &neighs, &nneighs);
			//if (VERBOSE)printf("N[%d]", s[i]); print_int("", neighs, nneighs);
			for (j = 0; j < nneighs; j++) {
				k = neighs[j];
				// update neighs scores
				scores[k]++;
				//if (VERBOSE)  printf("score[%d] = %d\n", k, scores[k]);
				//if (VERBOSE) printf("isloc[%d] = %d\n", k, isloc[k]);
				// check if k coords can be calculated
				if (scores[k] >= n && isloc[k] == 0) {
					// k can be fixed
					s_new[sn_new++] = k;
					// set x[k] location
					isloc[k] = base_anchor(G, k, isloc, b);
					//printf("Loc %d with ", k); print_int("b", b, 4);
					//if (VERBOSE) print_int("s_new", s_new, sn_new);
					if (isloc[k]) {
						base_update_x(G, b, k, x);
						graph_update_error_neighs(G, k, x, isloc);
					}
					else {
						printf("The node %d could not be loc.\n", k);
						exit(EXIT_FAILURE);
					}
				}
			}
		}

		// update s
		swap = s;
		s = s_new;
		sn = sn_new;
		s_new = swap;
	} while (sn > 0);

	// print outputs
	printf("   Number of fixed points: %d/%d\n", nfixed, G->nnodes);
	if (nfixed < G->nnodes && nfixed < 30) {
		printf("     NotFixed: ");
		for (i = 0; i < G->nnodes; i++) {
			if (isloc[i] == 0) {
				printf("%d ", i);
				x[3 * i + 0] = NAN;
				x[3 * i + 1] = NAN;
				x[3 * i + 2] = NAN;
			}
		}
	}
	printf("\n");
	graph_print_error_stats(G, isloc);

	free(s);
	free(s_new);
	free(scores);
	free(isloc);
}

void nodes_save(char *filename, double *x, int n) {
	int i;
	FILE *fid = fopen(filename, "w");
	printf("Saving file: %s\n", filename);
	fprintf(fid, "x,y,z\n");
	for (i = 0; i < n; i++) {
		fprintf(fid, "%g,%g,%g\n", x[3 * i], x[3 * i + 1], x[3 * i + 2]);
	}
	fclose(fid);
}

int main(int argc, char **argv) {
	int num_procs, num_threads;
	Graph G;
	char isclique;
	double *x;
	char filename[256];

	// set base
	int B[4] = { 0,1,2,3 }, n = 4;

	num_procs = omp_get_num_procs();
	omp_set_num_threads(num_procs - 1);
	num_threads = omp_get_max_threads();
	printf("Number of procs: %d\n", num_procs);
	printf("Running with %d threads\n", num_threads);

	//// READ GRAPH ////////////////////////////////////////////////////////
	sprintf(filename, "%s.csv", argv[1]);
	double tic = omp_get_wtime();
	graph_read(filename, &G);
	tic = omp_get_wtime() - tic;
	printf("   ETA: %3.2g seconds\n", tic);

	//// INIT X ////////////////////////////////////////////////////////////		
	x = (double*)malloc(sizeof(double)* (3 * G.nnodes));
	if (x == NULL) {
		printf("Memory could not be allocated\n");
		exit(EXIT_FAILURE);
	}

	////////////////////////////////////////////////////////////////////////	
	printf("Init base\n");
	tic = omp_get_wtime();
	base_init(&G, B, x);
	tic = omp_get_wtime() - tic;
	printf("   ETA: %3.2g seconds\n", tic);

	///// SET X COORDS ///////////////////////////////////////////////////////
	printf("Set x coords\n");
	print_int("   B", B, 4);
	tic = omp_get_wtime();
	isclique = base_isclique(&G, B, n);
	printf("   isclique: %s\n", isclique ? "yes" : "no");
	if (isclique) {
		base_locx(&G, B, x);
	}
	tic = omp_get_wtime() - tic;
	printf("   ETA: %3.2g seconds\n", tic);

	// Save solution
	sprintf(filename, "%s_xsol.csv", argv[1]);
	nodes_save(filename, x, G.nnodes);
	return EXIT_SUCCESS;
}