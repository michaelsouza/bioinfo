function [done,msg] = sph_stop(G, x, dij_tol, maxit, nit_global, df, dx)
    msg = '';
    if is_solution(G,x,dij_tol)
        msg = 'A solution has been found';
    elseif nit_global >= maxit
        msg = 'The maximum number of iterations has been reached';
    elseif dx < 1E-8
        msg = 'Solution stagnation';
    elseif abs(df) < 1E-8
        msg = 'Function stagnation';
    end
    done = ~isempty(msg);
end

function bval = is_solution(G,x,dij_tol)
bval = true;
for k = 1:G.nedges
    i = G.i(k);
    j = G.j(k);
    lij = (1-dij_tol) * G.l(k);
    uij = (1+dij_tol) * G.u(k);
    xi = x(:,i);
    xj = x(:,j);
    dij = norm(xi - xj);
    f_lij = max(lij - dij, 0);
    f_uij = max(dij - uij, 0);
    eij = max([f_lij, f_uij]);
    if eij > 0
        bval = false;
        return
    end
end
end