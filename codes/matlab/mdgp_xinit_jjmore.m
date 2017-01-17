function x = mdgp_xinit_jjmore(G)
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
	x = zeros(3,nnodes);
    for i = 1:nnodes
        xi = x(:,i);
        [~,J] = find(D(i,:));
        for j = J
            % random spherical coords centered on x(i)
            dij   = D(i,j);
            phi   = rand * pi;
            theta = rand * pi;
            x(:,j) = xi + dij * [sin(theta) * cos(phi); sin(theta) * sin(phi); cos(theta)];
        end
    end
    
    % centering in the origin
    c = mean(x, 2);
    for i = 1:3
        x(i,:) = x(i,:) - c(i);
    end
end