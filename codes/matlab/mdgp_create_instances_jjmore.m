function mdgp_create_instances_jjmore(nnodes)
    %% get pdb struct
    disp('Retrieving data from PDB');
    pdb = getpdb('1GPV');

    %% get coords
    x = [pdb.Model.Atom.X];
    y = [pdb.Model.Atom.Y];
    z = [pdb.Model.Atom.Z];
    x = [x',y',z'];

    %% restriction to chain A
    chainID = [pdb.Model.Atom.chainID] == 'A';
    resSeq = [pdb.Model.Atom.resSeq];
    resSeq = resSeq(chainID);
    x = x(chainID,:);

    %% calculate distances
    disp('Creating distance array');
    xi = zeros(nnodes^2,1);
    xj = zeros(nnodes^2,1);
    xx = zeros(nnodes^2,1);
    nedges = 0;
    for i = 1:nnodes
        ri = resSeq(i);
        for j = (i+1):nnodes
            rj = resSeq(j);
            if abs(ri - rj) <= 1
                nedges = nedges + 1;
                xi(nedges) = i;
                xj(nedges) = j;
                xx(nedges) = norm(x(i,:) - x(j,:));
            end
        end
    end
    xi = xi(1:nedges);
    xj = xj(1:nedges);
    xx = xx(1:nedges);

    %% creating instances
    eps_xx = [0.0,0.04,0.08,0.12,0.16];
    for k = 1:length(eps_xx)
        eps_xx_k = eps_xx(k);
        L = (1 - eps_xx_k) * xx;
        U = (1 + eps_xx_k) * xx;
        fname_prefix = sprintf('mdpg_1GPV_N%03d_EPS%3.2f', nnodes, eps_xx_k);
        fprintf('Saving instance: %s\n', fname_prefix);
        disp('   Saving instance solution');
        fid = fopen([fname_prefix '_xsol.csv'], 'w');
        fprintf(fid,'x,y,z\n');
        for i = 1:nnodes
            fprintf(fid,'%g,%g,%g\n',x(i,:));
        end
        fclose(fid);
        disp('   Saving instance edges');
        fid = fopen([fname_prefix '.csv'], 'w');
        fprintf(fid,'I,J,L,U\n');
        for i = 1:nedges
            fprintf(fid,'%g,%g,%g,%g\n',xi(i),xj(i),L(i),U(i));
        end
        fclose(fid);
    end
end