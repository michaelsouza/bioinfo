function [x,f] = solver_opt(problem)
natoms = max(max(problem.i), max(problem.j));

% starting point (vector)
x = rand(3 * natoms, 1);
% x = problem.x';
% x = x(:);

% set objective function
fobj = @(x)dgp_fobj(problem, x);

% solve optimization formulation
[x,f] = fminunc(fobj, x);

% convert output to matrix
x = reshape(x, 3, natoms)';
end