function x = sph(G,x,options)
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
	fprintf('SPH\n')
	fprintf('   alpha = %g\n', alpha);
	fprintf('   tau   = %g\n', tau);
	fprintf('   rho   = %g\n', rho);
	fprintf('   maxit = %g\n', maxit);
	fprintf('   dtol  = %g\n', dtol);
   
	fx = mdgp_fobj(G,x,dtol);
	fprintf('   fx    = %g\n\n', fx);
	% check_xsol(G,x, detailed=True)
	nit_local  = 0;
	nit_global = 0;
    
    if options.plot
        close all
        FIG  = figure;
        AXES = axes('Parent',FIG);
        view(AXES,[32 16]);
        box  on;
        grid on;
    end

    done = false;
    while ~done
        nit_global = nit_global + 1;
        fx_old = fx;
        x_old  = x;
        % defining and solving unconstrained problem
        f = @(y) mdgp_fobj_smooth(alpha, tau, G, y);
        x = x(:); % matrix to vector
        tic
        [x,~,~,output] = fminunc(f, x, options.optim);
        telapsed = toc;
        x  = reshape(x,3,[]); % vector to matrix
        fx = mdgp_fobj(G,x,dtol);
        nit_local = nit_local + output.iterations;
        
        df = (fx - fx_old) / fx_old;
        dx = norm(x - x_old) / max(norm(x), 1);
        
        % check if the points are aligned
        if options.plot
            view_coords(x);
            drawnow
        end
        
        if mod(nit_global,20) == 1
            fprintf(' iter    fx        df       dx      nit   rho     time\n')
        end
        if mod(nit_global,1) > -1
            fprintf('%5d %5.2E % 5.2E %5.2E %4d  %5.2E %3.2f\n', ...
                nit_global, fx, df, dx, output.iterations, tau, telapsed);
        end
        % stop criteria
        [done,msg] = sph_stop(maxit,nit_global,fx,df,dx);
        tau = tau * rho^(1 + max(0, -log(dx)));
        % tau = tau * rho;
    end
	fprintf('%5d %5.2E % 5.2E %5.2E %4d  %5.2E %3.2f\n', ...
		nit_global, fx, df, dx, output.iterations, tau, toc);
	fprintf('\nFinal results\n')
	fprintf('   nit global = %g\n', nit_global);
	fprintf('   nit local  = %g\n', nit_local);
	fprintf('   max_error  = %g\n', mdgp_calculate_erros(G,x));
    fprintf('   tElapsed   = %3.2fs\n', toc(tstart));
    fprintf('   message    = %s\n', msg);
    close 
end