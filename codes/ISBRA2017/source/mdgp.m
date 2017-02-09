function [x_sph,x_rot,G,xsol] = mdgp(fname, calc_rmsd)
    if nargin < 2
        calc_rmsd = false;
    end
    
    fid = fopen([fname,'_log.csv'], 'w');
    fprintf(fid,'method,fobj,nit_global,nit_local,total_time,msg\n');
    
    [G, xsol] = mdgp_load_problem(fname);
    x_ini     = mdgp_xinit_rotors(G);
    
    options = sph_options();
    
    options.optim.GradObj = 'on';
    fprintf('\n\nOptions\n');
    fprintf('   GradObj = %s\n', options.optim.GradObj);
    fprintf('   Plot    = %s\n\n', options.plot);
    [x_sph,fx,nit_global,nit_local,total_time_elapsed,msg] = sph(G,x_ini,options);
    fprintf(fid,'sph,%g,%g,%g,%g,"%s"\n',fx,nit_global,nit_local,total_time_elapsed,msg);
   
    options.optim.GradObj = 'on';
    fprintf('\n\nOptions\n');
    fprintf('   GradObj = %s\n', options.optim.GradObj);
    fprintf('   Plot    = %s\n', options.plot);
    [x_rot,fx,nit_global,nit_local,total_time_elapsed,msg] = sph_rotors(G,x_ini,options);
    fprintf(fid,'rot,%g,%g,%g,%g,"%s"\n',fx,nit_global,nit_local,total_time_elapsed,msg);
    
    fclose(fid);
    if calc_rmsd
        rmsd(x_sph,xsol,[],'SPH:');
        rmsd(x_rot,xsol,[],'SPH-ROTORS:');
    end
end