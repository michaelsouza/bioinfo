function [G,xsol] = mdgp_load_problem(fname)
	fid_problem = sprintf('%s%s', fname, '.csv');
	fprintf('Reading problem : %s\n', fid_problem)
	G = readtable(fid_problem);
	nedges = length(G.I);
	for k = 1:nedges
		i = G.I(k);
		j = G.J(k);
		lij = G.L(k);
		uij = G.U(k);
		if i > j
			fprintf('   Changing edge index order from (%d,%d) to (%d,%d)\n',i,j,j,i)
			G.I(k) = j;
			G.J(k) = i;
		end
		if lij > uij
			error('Inconsistent data (lij > uij): lij = %g and uij %g', lij, uij)
		end
	end
	nnodes = max(max(G.I),max(G.J));
	fid_xsol = sprintf('%s%s',fname, '_xsol.csv');
	fprintf('Reading solution: %s\n', fid_xsol)
	xsol = readtable(fid_xsol);
	xsol = [xsol.x';xsol.y';xsol.z'];
	G = struct('fname',fname,'nnodes',nnodes,'nedges',nedges,'i',G.I,'j',G.J,'l',G.L,'u',G.U);
	disp_problem(G)
end

function disp_problem(G)
	fprintf('Problem\n')
	disp(G)
end

