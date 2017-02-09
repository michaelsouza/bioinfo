function f = mdgp_fobj(G,x)
	f = 0.0;
	for k = 1:G.nedges
		lij = G.l(k);
		uij = G.u(k);
		xi = x(:,G.i(k));
		xj = x(:,G.j(k),:);
		dij = norm(xi - xj);
		f_lij = max(lij - dij, 0);
		f_uij = max(dij - uij, 0);
		f = f + f_lij + f_uij;
	end
end