function f = dgp_fobj(problem, x)
x = reshape(x, 3, [])';
f = 0;
nedges = problem.nedges;
for k = 1:nedges
    xi = problem.i(k);
    xj = problem.j(k);
    dk = norm(x(xi,:) - x(xj,:))^2;
    lk = problem.l(k)^2;
    uk = problem.u(k)^2;
    rk = max([(lk - dk)/lk,(dk - uk)/uk,0]);
    f = f + rk;
end
end