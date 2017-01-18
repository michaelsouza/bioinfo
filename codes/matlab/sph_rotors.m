function x = sph_rotors(G,xref,options)
    options.plot = ~strcmp(options.plot, 'off');
    tstart = tic;
	nedges = G.nedges;
	D     = (G.l + G.u) / 2.0;
	d     = sort(D);
	tau   = d(ceil(options.tau_level * nedges));
	alpha = options.alpha;
	rho   = options.rho;
	maxit = options.maxit;
	dtol  = options.dtol;
    acc   = options.acc;
	fx    = mdgp_fobj(G,xref,dtol);
    nit_local  = 0;
	nit_global = 0;
    
    if is_xinit_valid(G,xref) == false
        error('Invalid initial point.')
    end
   
    fprintf('SPH_ROT\n')
	fprintf('   alpha = %g\n', alpha);
	fprintf('   tau   = %g\n', tau);
	fprintf('   rho   = %g\n', rho);
	fprintf('   maxit = %g\n', maxit);
	fprintf('   dtol  = %g\n', dtol);
	fprintf('   fx    = %g\n', fx);
    fprintf('   acc   = %g\n\n', acc);
	    
    if options.plot
        close all
        FIG  = figure;
        AXES = axes('Parent',FIG);
        view(AXES,[32 16]);
        box  on;
        grid on;
    end

    x = xref;
    s = 4:G.nnodes;
    w = zeros(size(s));
    done = false;
    while ~done
        nit_global = nit_global + 1;
        fx_old = fx;
        w_old  = w;
        x_old  = x;
        % defining and solving unconstrained problem
        f = @(w) mdgp_fobj_smooth_rotors(alpha, tau, G, xref, s, w);
        tic
        [w,~,~,output] = fminunc(f, w, options.optim);
        telapsed = toc;
        x  = rotors_apply(xref,s,w);
        fx = mdgp_fobj(G,x,dtol);
        nit_local = nit_local + output.iterations;
        
        dw = norm(w - w_old) / max(norm(w_old), 1);
        df = (fx - fx_old) / fx_old;
        dx = norm(x - x_old) / max(norm(x_old), 1);
        speed = ceil(1 + max(0, -acc * log(dx)));
        % check if the points are aligned
        if options.plot
            view_coords(xref);
            drawnow
        end
        
        if mod(nit_global,20) == 1
            fprintf(' iter    fx        df       dx      dw      nit    rho     speed   time\n')
        end
        if mod(nit_global,1) > -1
            fprintf('%5d %5.2E % 5.2E %5.2E %5.2E %4d  %5.2E %5d %6.2f\n', ...
                nit_global, fx, df, dx, dw, output.iterations, tau, speed, telapsed);
        end
        % stop criteria
        [done,msg] = sph_stop(maxit,nit_global,fx,df,dx);
        tau = tau * rho^speed;
    end
    fprintf('%5d %5.2E % 5.2E %5.2E %5.2E %4d  %5.2E %5d %6.2f\n', ...
        nit_global, fx, df, dx, dw, output.iterations, tau, speed, telapsed);
    fprintf('\nFinal results\n')
	fprintf('   nit global = %g\n', nit_global);
	fprintf('   nit local  = %g\n', nit_local);
	fprintf('   max_error  = %g\n', mdgp_calculate_erros(G,xref));
    fprintf('   tElapsed   = %3.2fs\n', toc(tstart));
    fprintf('   message    = %s\n', msg);
    close 
end

function b = is_xinit_valid(G,x)
for k = 1:G.nedges
    i = G.i(k);
    j = G.j(k);
    if j - i > 2
        continue
    end
    lij = G.l(k);
    uij = G.u(k);
    xi = x(:,i);
    xj = x(:,j);
    dij = norm(xi-xj);
    b = lij - dij < 1E-8 && dij - uij < 1E-8;
    if b == false
        return
    end
end
end