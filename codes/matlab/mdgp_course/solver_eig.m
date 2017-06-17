function x = solver_eig(problem)
natoms = max(max(problem.i), max(problem.j));

% create distance matrix
D = sparse(problem.i, problem.j, problem.d, natoms, natoms);

% create distance squared matrix
M = zeros(natoms);
for i = 1:natoms
    for j = i:natoms
        M(i,j) = (D(1,i)^2 + D(1,j)^2 - D(i,j)^2) / 2;
    end
end
M = M + M' - diag(diag(M));

% get eigenvalues
[x,lambda] = eigs(M,3);

% set solution
x = x * sqrt(lambda);
end