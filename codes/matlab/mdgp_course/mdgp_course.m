%% Class 1.A: Creating an instance from a known solution ==================
x = rand(5,3);
problem = create_mdgp_instance(x);
plot_dgp_problem(problem);

%% Class 1.B: Creating a random instance ==================================
problem = create_mdgp_instance_rand(5);
plot_dgp_problem(problem);

%% Class 1.C: Calculating the protein angles ==============================
x = [0 0 0;
     1 0 0;
     1 1 0;
     1 1 1];
problem = create_mdgp_instance(x);
plot_dgp_problem(problem);
view_molecule(problem.x);
calc_protein_angles(problem.x);

%% Class 1.C: Checking the quality of a MDGP solution =====================
problem  = create_mdgp_instance_rand(10);
check_dgp_solution(problem, problem. x);

x = problem.x + 0.01 * rand(problem.natoms, 3);
check_dgp_solution(problem, x);

%% Class 1.D: Visualizing a protein (solution) ============================
problem = create_mdgp_instance_rand(10);
% plot_molecule(problem.x);
view_molecule(problem.x);

%% Class 1.E: Solving a DGP instance using nonlinear equations ============
problem = create_mdgp_instance_rand(30);
x = solver_sys(problem);
check_dgp_solution(problem, x);

%% Class 2.A: Creating and solnving a small problem =======================
x = [0 0 0;
    1 0 0;
    1 1 0;
    1 1 -1;
    1 0 -1];
% i = [1 2 3 4, 1, 2, 3, 1];
% j = [2 3 4 5, 3, 4, 5, 5];
problem = create_mdgp_instance(x);
plot_dgp_problem(problem);
y = solver_sys(problem);
check_dgp_solution(problem, y);

%% Class 2.B: RMSD - Check if the solution is different ===================
% plot y and x
figure
view_molecule(x);
hold on
view_molecule(y);

% check if x and y represent the same structure
[rmsd,x_new,y_new] = calc_rmsd(x,y);
fprintf('rmsd(x,y) = %3.2e\n', rmsd);

% plot y and x
figure
view_molecule(x_new);
hold on
view_molecule(y_new);

%% Class 3.A: Solving the optimization formulation of DGP =================
problem = create_mdgp_instance_rand(20);
tic
[x,f] = solver_opt(problem);
fprintf('Elapsed time %3.2f seconds\n', toc);
check_dgp_solution(problem, x);

%% Class 3.B: Solving a DGP instance with all constraints =================
problem = create_mdgp_instance_rand(100, inf);
tic
x = solver_eig(problem);
fprintf('Elapsed time %3.2f seconds\n', toc);
check_dgp_solution(problem, x);

%% Class 3.C: Solving a system with three quadratic equations =============
% creates a random problem
a = rand(1,3);
b = rand(1,3);
c = rand(1,3);
x = rand(1,3);
da = norm(a - x);
db = norm(b - x);
dc = norm(c - x);
% get the intersection
y = intersect3spheres(a,da,b,db,c,dc);
clc
plot_vector('x ',x);
plot_vector('y1',y(1,:));
plot_vector('y2',y(2,:));

%% Class 3.D: Solving a system with four quadratic equations ==============
a = rand(1,3);
b = rand(1,3);
c = rand(1,3);
d = rand(1,3);
x = rand(1,3);
da = norm(a - x);
db = norm(b - x);
dc = norm(c - x);
dd = norm(d - x);
y = intersect4spheres(a,da,b,db,c,dc,d,dd);
clc
plot_vector('x',x);
plot_vector('y',y);

%% Class 5.D: Find a solution to the following instances
x_sol = [0 0 0; 1 0 0; 1 1 0; 1 1 1; 1 0 1; 0 0 1];

%% 5.DA: All information
i = [1 1 1 1 1 2 2 2 2 3 3 3 4 4 5];
j = [2 3 4 5 6 3 4 5 6 4 5 6 5 6 6];
problem = create_mdgp_instance(x_sol,i,j);
plot_dgp_problem(problem)

%% 5.DB: Data is incomplete, but it is exact
i = [1 1 1 1 1 2 2 2 3 3 3 4 4 5];
j = [2 3 4 5 6 3 4 5 4 5 6 5 6 6];
problem = create_mdgp_instance(x_sol,i,j);
plot_dgp_problem(problem)

%% 5.DC: data is incomplete, but it is exact (fewer edges)
i = [1 1 1 1 2 3 3 3 4 4 5];
j = [2 3 4 6 3 4 5 6 6 5 6];
problem = create_mdgp_instance(x_sol,i,j);
plot_dgp_problem(problem)

%% 5.DD: data is incomplete and imprecise
i   = [1 1   1   1 2 2   2   2 3 3   3 4 4 5];
j   = [2 3   4   6 3 4   5   6 4 5   6 5 6 6];
eps = [0 0 0.1 0.1 0 0 0.1 0.1 0 0 0.1 0 0 0];
problem = create_mdgp_instance(x_sol,i,j,eps);
plot_dgp_problem(problem);
y = solver_opt(problem);
check_dgp_solution(problem,y);
problem = create_mdgp_instance(x_sol,i,j,eps);
plot_dgp_problem(problem);

%% Solver: BP
% create instance (y = solution)
i = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 6];
j = [2 3 4 5 6 3 4 5 4 5 6 5 6 7 6 7 7];
y = [0,0,0;1,0,0;1,1,0;1,1,1;0,1,1;0,0,1;0,0,2];
problem = create_mdgp_instance(y,i,j);
close all
subplot(1,2,1)
view_molecule(y)
subplot(1,2,2)
plot_dgp_problem(problem);

% call solver
[X,B] = solver_bp(problem);

%% Solver: IBP
% create instance (y = solution)
i = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 6];
j = [2 3 4 5 6 3 4 5 4 5 6 5 6 7 6 7 7];
y = [0,0,0;1,0,0;1,1,0;1,1,1;0,1,1;0,0,1;0,0,2];
dij_eps = zeros(size(i));
for k = 1:length(i)
if abs(i(k)-j(k)) > 2
    dij_eps(k) = 0.1;
end
end
problem = create_mdgp_instance(y,i,j,dij_eps);
close all
subplot(1,2,1)
view_molecule(y)
subplot(1,2,2)
plot_dgp_problem(problem);

% call solver
[x,p,b] = solver_ibp(problem,3);

%% Solver: IBPR

% create instance (y = solution)
i = [1 1 1 1 2 2 2 3 3 4];
j = [2 3 4 5 3 4 5 4 5 5];
y = [0,0,0;
    -1,0,0;
    -3/2,sqrt(3)/2,0;
    -1.311,1.552,0.702;
    -0.728,2.337,0.491];
problem = create_mdgp_instance(y,i,j);
problem.l(3) = problem.d(3); problem.u(3) = problem.d(3);
problem.l(4) = 2.45; problem.u(4) = 2.55;
problem.l(7) = 2.20; problem.u(7) = 2.60;

% close all
% subplot(1,2,1)
% view_molecule(y)
% subplot(1,2,2)
plot_dgp_problem(problem);

[x,p,b] = solver_ibpr(problem, 3);