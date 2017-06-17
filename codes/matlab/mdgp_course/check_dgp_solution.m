function relative_residues = check_dgp_solution(problem, x)
nedges   = length(problem.d);
relative_residues = zeros(nedges, 1);
for k = 1:nedges
    xi = problem.i(k);
    xj = problem.j(k);
    dk = norm(x(xi,:) - x(xj,:));
    lk = problem.l(k);
    uk = problem.u(k);
    relative_residues(k) = max([(lk - dk) / lk, (dk - uk) / uk, 0]);
end
natoms = problem.natoms;
R = sparse(problem.i, problem.j, relative_residues, natoms, natoms);
figure
image(R,'CDataMapping','scaled')
xlabel('atom index')
ylabel('atom index')
title('Relative Residue');
colormap jet;
colorbar

fprintf('Checking MDGP solution:\n');
fprintf('  Min  Relative Residue: %3.2e\n', min(relative_residues));
fprintf('  Max  Relative Residue: %3.2e\n', max(relative_residues));
fprintf('  Mean Relative Residue: %3.2e\n', mean(relative_residues));
end