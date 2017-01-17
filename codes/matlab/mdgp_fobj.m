function f = mdgp_fobj(G,x,dij_tol)
	if nargin < 2
		dij_tol = 0.0;
	end
	f = 0.0;
	for k = 1:G.nedges
		lij = (1-dij_tol) * G.l(k);
		uij = (1+dij_tol) * G.u(k);
		xi = x(:,G.i(k));
		xj = x(:,G.j(k),:);
		dij = norm(xi - xj);
		f_lij = max(lij - dij, 0);
		f_uij = max(dij - uij, 0);
		f = f + f_lij + f_uij;
	end
end