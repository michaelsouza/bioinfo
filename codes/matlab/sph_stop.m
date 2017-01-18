function [done,msg] = sph_stop(maxit,nit_global,fx,df,dx)
    msg = '';
    if fx < 1E-8
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