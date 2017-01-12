function max_eij = mdgp_calculate_erros(G, x, verbose)
	if nargin < 3
		verbose = false;
	end
	max_eij = 0.0;
	for k = 1:G.nedges
		i = G.i(k);
		j = G.j(k);
		lij = G.l(k);
		uij = G.u(k);
		xi = x(i,:);
		xj = x(j,:);
		dij = norm(xi - xj);
		f_lij = max(lij - dij, 0);
		f_uij = max(dij - uij, 0);
		eij = max([f_lij, f_uij]);
		if eij > max_eij
			max_eij = eij;
		end
		if verbose
			fprintf('(%d,%d) [%3.2f,%3.2f] = %6.2f (%6.3f)\n', i, j, lij, uij, dij, eij);
		end
	end
	if verbose
		fprintf('MaxError: %g\n', max_eij)
	end
end