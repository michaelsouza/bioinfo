function options = sph_options(varargin)
    options.optim  = optimoptions('fminunc',...
        'GradObj','on',...
        'DerivativeCheck','off',...
        'Display','off',...
        'Algorithm','quasi-newton');
    options.plot      = 'off';
    options.tau_level = 0.5;
	options.alpha     = 0.5;
	options.rho       = 0.99;
	options.maxit     = 100;
	options.dtol      = 1E-2;
    options.acc       = 20;
end