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
