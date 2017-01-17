x = sym('x',[3,1],'real');
y = sym('y',[3,1],'real');
tau = sym('tau','real');
alpha = sym('alpha','real');
lij = sym('lij','real');
uij = sym('uij','real');
dij = sqrt(tau^2 + sum((x - y).^2));
elij = lij - dij;
euij = dij - uij;
phi = @(x) alpha * x + sqrt((alpha * x).^2 + tau.^2);
fij = phi(elij) + phi(euij);

%% create files
fprintf('Creating files\n');
g = gradient(fij,[x;y]);
matlabFunction([fij;g], 'File', 'fij.m');