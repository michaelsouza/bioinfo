function f = dgp_sys(problem, x)
x = reshape(x, 3, [])';
f = zeros(problem.nedges, 1);
for k = 1:problem.nedges
    i = problem.i(k);
    j = problem.j(k);
    d = problem.d(k);
    f(k) = norm(x(i,:) - x(j,:))^2 - d^2;
end
end