function [x,G,xsol] = mdgp()
fname = '../../instances/mdpg_1GPV_N008_EPS0.16';
[G, xsol] = mdgp_load_problem(fname);
x = sph(G);
[value, Q, x, y] = rmsd(x,xsol,[],'Coords', {'x', 'xsol'});
fprintf('rmsd = %g\n', rmsd(x, xsol)); 

%% Check functions
% check_theta()
% check_phi()
% check_fobj_smooth(G,xsol)
% check_xinit(G)
% check_rotate()
% check_rotors_apply()
% check_rotors_diff()
%    check_fobj_smooth_rot(G,xsol)
%    sph(G)
end

function disp_vec(name, x)
    fprintf('%s =', name)
    for i = 1:length(x)
        fprintf(' % g', x(i))
    end
    fprintf('\n');
end

function disp_mat(name, x)
	fprintf('%s =\n', name)
	[nrow,ncol] = size(x);
	for i = 1:nrow
		for j = 1:ncol
			fprintf(' % 8.5g', x(i, j))
		end
	fprintf('\n');    
	end
end

function x = rotate(vec_p,vec_d,theta,vec_x)
	% Rotates the vec_x around the axis vec_d anchored at the point vec_p.
	%References:
	%[1] http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisrotate/

	% unpacking
	vec_d = vec_d / norm(vec_d); % normalizing
	a=vec_p(1); b=vec_p(2); c=vec_p(3); % base point
	u=vec_d(1); v=vec_d(2); w=vec_d(3); % line direction
	x=vec_x(1); y=vec_x(2); z=vec_x(3); % point to be rotated

	cos_theta = cos(theta);
	sin_theta = sin(theta);
	vec_d_dot_vec_x = u*x+v*y+w*z;
	omt = 1 - cos_theta;
	C1 = (a*(v^2+w^2) - u*(b*v+c*w-vec_d_dot_vec_x));
	C2 = (b*(u^2+w^2) - v*(a*u+c*w-vec_d_dot_vec_x));
	C3 = (c*(u^2+v^2) - w*(a*u+b*v-vec_d_dot_vec_x));
	S1 = (-c*v + b*w - w*y + v*z);
	S2 = ( c*u - a*w + w*x - u*z);
	S3 = (-b*u + a*v - v*x + u*y);
	x = [C1*omt + x*cos_theta + S1*sin_theta,...
		C2*omt + y*cos_theta + S2*sin_theta,...
		C3*omt + z*cos_theta + S3*sin_theta];
end

function g = rotate_diff(vec_p,vec_d,theta,vec_x)
	% unpacking
	vec_d = vec_d / norm(vec_d); % normalizing
	a=vec_p(1); b=vec_p(2); c=vec_p(3); % base point
	u=vec_d(1); v=vec_d(2); w=vec_d(3); % line direction
	x=vec_x(1); y=vec_x(2); z=vec_x(3); % point to be rotated

	cos_theta = cos(theta);
	sin_theta = sin(theta);
	vec_d_dot_vec_x = u*x+v*y+w*z;
	C1 = (a*(v^2+w^2) - u*(b*v+c*w-vec_d_dot_vec_x));
	C2 = (b*(u^2+w^2) - v*(a*u+c*w-vec_d_dot_vec_x));
	C3 = (c*(u^2+v^2) - w*(a*u+b*v-vec_d_dot_vec_x));
	S1 = (-c*v + b*w - w*y + v*z);
	S2 = ( c*u - a*w + w*x - u*z);
	S3 = (-b*u + a*v - v*x + u*y);
	g = [(C1-x)*sin_theta+S1*cos_theta,...
		(C2-y)*sin_theta+S2*cos_theta,...
		(C3-z)*sin_theta+S3*cos_theta];
end

function x = rotors_apply(x,wi,wx)
	x = array2matrix(x);
	% applying rotates
	k_old=-1;
	for j = 1:length(wi)
		k = wi(j);
		a = wx(j);
		if k <= k_old
			error('The wi array must be sorted in ascending order.')
		end
		for i = k:length(x)
			x(i,:) = rotate(x(k-1,:),x(k-1,:)-x(k-2,:),a,x(i,:));
		end
	end
end

function g = rotors_diff(x,wi,wx)
	x = array2matrix(x);
	% diff rotate
	g = zeros(length(x),3,length(wi));
	for j = 1:length(wi)
		k = wi(j);
		a = wx(j);
		for i = k:length(x)
			g(i,:,j) = rotate_diff(x(k-1,:),x(k-1,:)-x(k-2,:),a,x(i,:));
		end
	end
end

function x = xinit_jjmore(G)
	nnodes = G.nnodes;
	nedges = G.nedges;

	% D : sparse matrix (CSR)
	Di = zeros(nnodes^2,1);
	Dj = Di;
	Dx = Di;
	nz = 0;
	for k = 1:nedges
		lij = G.l(k);
		uij = G.u(k);
		nz = nz + 1;
		Di(nz) = G.i(k);
		Dj(nz) = G.j(k);
		Dx(nz) = (lij + uij) / 2.0;
	end
	D = sparse(Di(1:nz),Dj(1:nz),Dx(1:nz),nnodes,nnodes);
	x = zeros(nnodes,3);
	for i = 1:nnodes
		xi = x(i,:);
		[~,J] = find(D(i,:));
		for j = J
			% random spherical coords centered on x(i)
			dij   = D(i,j);
			phi   = rand * pi;
			theta = rand * pi;
			x(j,:) = xi + dij * [sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)];
		end
	end
end

function x = xinit_rotors(G)
	nnodes = G.nnodes;
	d = sparse(G.i, G.j, (G.l + G.u) / 2, nnodes, nnodes);
	% planar angles
	planar_angle_ijk = @(i,j,k) acos((d(i,j)^2 + d(i,k)^2 - d(j,k)^2)/(2*d(i,j)*d(j,k)));
	theta = zeros(nnodes, 1);
	for i = 3:nnodes
		theta(i) = planar_angle_ijk(i-2,i-1,i);
	end
	% random dihedral angles
	w = rand(nnodes, 1) * 2 * pi;
	
	% set random positions
	x = zeros(nnodes, 3);
	x(2,:) = [-d(1,2),0,0];
	x(3,:) = [-d(1,2) + d(2,3) * cos(theta(3)), dBC * sin(theta(3)) , 0];
	for i = 4:nnodes
		A = x(i-3,:);
		B = x(i-2,:);
		C = x(i-1,:);
		n = cross(A-C, B-C); % normal to the plane ABC
		dBC = d(i-2, i-1);
		dCD = d(i-1, i);
		D = C + dDC * (C - B) / dBC;
		% rotation on the plane ABC
		D = rotate(C, n, theta(i), D); 
		% rotation around BC
		x(i,:) = rotate(B, C - B, w(i), D);
	end
end

function f = fobj(G,x,dij_tol)
	if nargin < 2
		dij_tol = 0.0;
	end
	f = 0.0;
	for k = 1:G.nedges
		lij = (1-dij_tol) * G.l(k);
		uij = (1+dij_tol) * G.u(k);
		xi = x(G.i(k),:);
		xj = x(G.j(k),:);
		dij = norm(xi - xj);
		f_lij = max(lij - dij, 0);
		f_uij = max(dij - uij, 0);
		f = f + f_lij + f_uij;
	end
end

function [f,g] = theta(tau, x)
	f = sqrt(tau^2 + norm(x)^2);
	if nargout > 1
		g = x / f;
	end
end

function [f,g] = phi(alpha, tau, x)
	f = alpha * x + sqrt((alpha * x)^2 + tau^2);
	if nargout > 1
		g = alpha * f / (f - alpha * x);
	end
end

function [f,g] = fobj_smooth(alpha, tau, G, y)
	x = array2matrix(y);
	% initializations
	f = 0.0;
	g = zeros(length(x), 3);
	for k = 1:G.nedges
		lij = G.l(k);
		uij = G.u(k);
		i = G.i(k);
		j = G.j(k);
		xi = x(i,:);
		xj = x(j,:);
		[theta_ij, g_theta]  = theta(tau, xi - xj);
		[phi_lij, g_phi_lij] = phi(alpha, tau, lij - theta_ij);
		[phi_uij, g_phi_uij] = phi(alpha, tau, theta_ij - uij);
		f =  f + phi_lij + phi_uij;
		if nargout > 1
			g(i,:) = g(i,:) + (g_phi_uij - g_phi_lij) * g_theta;
			g(j,:) = g(j,:) - g(i,:);
		end
	end
	if nargout > 1
		g = matrix2array(g);
	end
end

function [f,g] = fobj_smooth_rot(alpha, tau, G, x, wi, wx)
	x = array2matrix(x);

	% applying rotates
	x = rotors_apply(x, wi, wx)

	% eval smooth function
	f = 0.0;
	for k = 1:G.nedges
		lij = G.l(k);
		uij = G.u(k);
		i = G.i(k);
		j = G.j(k);
		xi = x(i,:);
		xj = x(j,:);
		theta_ij, g_theta  = theta(tau, xi - xj);
		phi_lij, g_phi_lij = phi(alpha, tau, lij - theta_ij);
		phi_uij, g_phi_uij = phi(alpha, tau, theta_ij - uij);
		f = f + phi_lij + phi_uij;
	end
end

function x = sph(G)
	nedges = G.nedges;
	D = (G.l + G.u) / 2.0;
	d = sort(D);
	tau   = d(ceil(0.75 * nedges));
	alpha = 0.5;
	rho   = 0.99;
	maxit = 1000;
	ftol  = 1E-8;
	dtol  = 1E-2;
	fprintf('SPH\n')
	fprintf('   alpha = %g\n', alpha);
	fprintf('   tau   = %g\n', tau);
	fprintf('   rho   = %g\n', rho);
	fprintf('   maxit = %g\n', maxit);
	fprintf('   ftol  = %g\n', ftol);
	fprintf('   dtol  = %g\n', dtol);
	fopts = optimoptions('fminunc','GradObj','on','Display','off','Algorithm','quasi-newton' );
	x  = xinit_jjmore(G);
	fx = fobj(G,x,dtol);
	fprintf('   fx    = %g\n\n', fx);
	% check_xsol(G,x, detailed=True)
	nit_local = 0;
	for nit_global = 1:maxit
		fx_old = fx;
		x_old = x;
		% defining and solving unconstrained problem
		f = @(y) fobj_smooth(alpha, tau, G, y);
		x = matrix2array(x); % matrix to vector
		tic
		[x,~,~,output] = fminunc(f, x, fopts);
		telapsed = toc;
		x = array2matrix(x);
		fx = fobj(G,x,dtol);
		nit_local = nit_local + output.iterations;
		
		df = (fx - fx_old) / fx_old;
		dx = max(norm(x - x_old)) / max(norm(x));
		if mod(nit_global,100) == 1
			fprintf(' iter    fx        df       dx      nit   rho     time\n')
		end
		if mod(nit_global,20) == 1
			fprintf('%5d %5.2E % 5.2E %5.2E %4d  %5.2E %3.2f\n', ...
				nit_global, fx, df, dx, output.iterations, tau, telapsed);
		end
		% stop criteria
		if fx < ftol
			break
		end
		tau = tau * rho;
	end
	fprintf('%5d %5.2E % 5.2E %5.2E %4d  %5.2E %3.2f\n', ...
		nit_global, fx, df, dx, output.iterations, tau, toc);
	fprintf('\nFinal results\n')
	fprintf('   nit global = %g\n', nit_global);
	fprintf('   nit local  = %g\n', nit_local);
	fprintf('   max_error  = %g\n', mdgp_calculate_erros(G,x));
end

% function check_theta()
% fprintf('Checking function theta\n')
% tau = 0.001;
% xij = [1,2,3];
% dij = norm(xij);
% fprintf('dij = %g\n', dij);
% [f, g] = theta(tau, xij);
% fprintf('f = %g\n', f)
% disp_vec('g   ', g)
% f = @(y) theta(tau, y);
% g_num = numdiff(f, xij);
% disp_vec('gnum', g_num);
% end

% function check_phi()
% fprintf('Checking function phi\n')
% alpha = 0.5;
% tau   = 0.001;
% y     = pi/3;
% fun   = max(y, 0);
% [f,g] = phi(alpha, tau, y);
% fprintf('fun = %g\n', fun)
% fprintf('phi = %g\n', f)
% f = @(y) phi(alpha, tau, y);
% fprintf('g     = %g\n', g);
% fprintf('g_num = %g\n', numdiff(f,y));
% end

% function check_fobj_smooth(G,xsol)
% dtol  = 0.01;
% fun   = fobj(G,xsol,dtol);
% tau   = 0.001;
% alpha = 0.5;
% [f,g] = fobj_smooth(alpha, tau, G, xsol);
% fprintf('fun = %g\n',fun)
% fprintf('f   = %g\n',f)
% f = @(y) fobj_smooth(alpha, tau, G, y);
% disp_vec('g    ', g);
% x = matrix2array(xsol);
% disp_vec('g_num', numdiff(f,x));
% end

% function check_fobj_smooth_rot(G,xsol)
% fprintf('Checking fobj_smooth_rot function')
%    y = rotate(xsol[])
%    fun = fobj(G,xsol)
%    tau   = 0.001
%    alpha = 0.5
%    wi = [4,5]
%    wx = [pi/8, pi/10]
%    f,g = fobj_smooth_rot(alpha, tau, G, xsol, wi, wx, grad=True)
%    fprintf('fun=',fun)
%    fprintf('f  =',f)
%    f = lambda wx: fobj_smooth_rot(alpha, tau, G, xsol, wi, wx, grad=False)
%    x = xsol.reshape((3 * len(xsol),))
%    fprintf('g    =', g)
%    fprintf('g_num=', numdiff(f,x))
% end


% function check_xinit(G)
% fprintf('Generating solution using xinit_random\n')
% x = xinit_random(G);
% check_xsol(G, x, true);
% fprintf('Generating solution using xinit_jjmore\n')
% x = xinit_jjmore(G);
% check_xsol(G, x, true);
% end

% function check_rotate()
% fprintf('Checking rotate function\n');
% [~,x] = load_problem('../../instances/mdpg_1GPV_N008_EPS0.16');
% x = array2matrix(x);
% p = x(3,:);
% d = x(3,:) - x(2,:);
% theta = pi / 3;
% disp_vec('y',rotate(p,d,theta,x(4,:)));
% f = @(t) rotate(p,d,t,x(4,:));
% disp_vec('g    ',rotate_diff(p,d,theta,x(4,:)));
% disp_vec('g_num',numdiff(f,theta));
% end

% function check_rotors_apply()
% fprintf('Checking rotors apply')
% [~,x] = load_problem('../../instances/mdpg_1GPV_N008_EPS0.16');
% wi = [4, 5];
% wx = [pi / 8, pi / 10];
% disp_mat('x', rotors_apply(x,wi,wx))
% end

% function check_rotors_diff()
	% fprintf('Checking rotors apply')
	% [~,x] = load_problem('../../instances/mdpg_1GPV_N008_EPS0.16');
	% wi = [4,5];
	% wx = [pi / 8, pi / 10];
	% x = matrix2array(x);
	% f = @(wx) rotors_apply(x,wi,wx);
	% fprintf('g=\n')
	% disp(rotors_diff(x,wi,wx));
	% fprintf('g_num=\n')
	% disp(numdiff(f,wx))
% end