function D = plot_dgp_problem(problem)
fprintf('natoms: %d\n', problem.natoms);
fprintf('nedges: %d\n', problem.nedges);
fprintf('id\t i\t j\t  d\t\t  l\t\t  u\n');
fprintf('%c%c%c %c%c%c %c%c%c %c%c%c%c%c\t%c%c%c%c%c\t%c%c%c%c%c\n',175 * ones(24,1))

D = sparse(problem.i,problem.j,problem.d,problem.natoms,problem.natoms);
L = sparse(problem.i,problem.j,problem.l,problem.natoms,problem.natoms);
U = sparse(problem.i,problem.j,problem.u,problem.natoms,problem.natoms);
[i,j,~] = find(D');
[j,ind] = sort(j);
i = i(ind);
for k = 1:problem.nedges
    ik = i(k);
    jk = j(k);
    fprintf('%3d\t%3d\t%3d\t%4.3f\t%4.3f\t%4.3f\n', ...
    k, jk, ik, full(D(jk,ik)), full(L(jk,ik)),full(U(jk,ik)));
end
spy(D);
hold on
plot(1:problem.natoms,1:problem.natoms);
xlabel('Atom ID');
ylabel('Atom ID');
grid on
title('Connectivity Matrix')
hold off
end