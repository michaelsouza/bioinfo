function x = solver_sys(problem)
natoms = problem.natoms;
x = rand(3 * natoms, 1);
% x = problem.x';x = x(:);
fsys = @(x) dgp_sys(problem, x);
x = fsolve(fsys, x);
x = reshape(x,3,[])';
end