%% sets the random seed in order to reproduce the results
rng(1);

%% gets all instance names
files  = dir('../instances/*N8*_ins.csv');
nfiles = length(files);
[~,i]  = sort([files.bytes]);
files  = files(i); % sort by size
for i = 1:nfiles
    files(i).name = ['../instances/' files(i).name];
end

%% execute all experiments
for i = 1:nfiles
    fname = strrep(files(i).name,'_ins.csv','');
    disp(fname);
    mdgp(fname);
end

%% compiling the results
fid = fopen('table_tests.csv','w');
fprintf(fid,'n,dmax,nedges,method,time_in_sec,fobj,nit_global,nit_local,solved,msg,uid\n');
for i = 1:length(files)
    fname  = strrep(files(i).name,'_ins.csv','_log.csv');
    if ~exist(fname, 'file')
        fprintf('skip: %s\n', fname);
        continue
    end
    table  = readtable(files(i).name);
    nedges = length(table.I);
    table  = readtable(fname); 
    tokens = strsplit(fname,'_');
    n           = tokens{2};
    dmax        = tokens{3};
    method      = table.method;
    fobj        = table.fobj;
    time_in_sec = table.total_time;
    nit_global  = table.nit_global;
    nit_local   = table.nit_local;
    msg         = table.msg;
    uid         = tokens{4};
    for j = 1:length(table.nit_global)
        fprintf(fid,'%s,%s,%d,%s,%g,%g,%g,%g,%d,"%s",%s\n',...
            n,dmax,nedges,method{j},time_in_sec(j),fobj(j),nit_global(j),...
            nit_local(j),strcmp(msg(j),'A solution has been found'),msg{j},uid);
    end
end
fclose(fid);

%% analysing results
table    = readtable('table_tests.csv');
query    = struct('n',{'N8'},'method',{'sph'});
selected = zeros(length(table.n),2);
for i = 1:length(table.n)
    selected(i,1) = strcmp(table.n(i),query.n);
    selected(i,2) = strcmp(table.method(i),query.method);
end
index = logical(prod(selected,2));