function [X,B] = solver_bp(problem, b, x, i)
xtol = 1E-5;

% set distance (bound) matrix
D = sparse(problem.i, problem.j, problem.d, problem.natoms, problem.natoms);
D = D + D';
L = sparse(problem.i, problem.j, problem.l, problem.natoms, problem.natoms);
L = L + L';
U = sparse(problem.i, problem.j, problem.u, problem.natoms, problem.natoms);
U = U + U';

% set initial point
x = x_start(D);

% start branch vector
if nargin > 2 || isempty(b)
    b = ones(1,problem.natoms);
    b(1:3) = 2;
else
end

% solution and branch cells
X = cell(1);
B = cell(1);
n = 0; % number of sols

i = 4;
while i > 3
    % get next x
    x(i,:) = x_next(D,b(i),x,i);
%   fprintf('x(%d) = [%3.2f,%3.2f,%3.2f]\n',i,x(i,:))
    % check partial solution
    errx = check_solution(L,U,x,i,xtol);
    prune = errx > xtol;
    if ~prune && i == problem.natoms
        n = n + 1;
        X{n} = x;
        B{n} = b;
        fprintf('b = '); fprintf('%g ', b(1:i)); fprintf('\n');
        fprintf('BP:Solution (%3d) found, natoms: %d\n', n, problem.natoms);
        % only the first solution is returned
        if nargin > 1
            X = x;
            B = b;
            return;
        end
    end
    if prune
        fprintf('b = '); fprintf('%g ', b(1:i)); fprintf('\n');
    end

    prune = prune || i == problem.natoms;    
    % backward
    if prune && b(i) == 2
        while i > 3 && b(i) == 2
            i = i - 1;
        end
    end
    
    % reset subtree
    if prune && b(i) == 1
        b(i) = 2;
        for j = (i+1):problem.natoms
            b(j) = 1;
        end
    end
   
    if ~prune
        i = i + 1;
    end
end

% no solution
if n == 0
    X = [];
    B = [];
    fprintf('BP:No solution\n');
end
end

function x = x_start(D)
x = zeros(length(D), 3);
x(2,1) = D(1,2);
theta  = acos((D(1,2)^2 + D(1,3)^2 - D(2,3)^2) / (2 * D(1,2) * D(1,3)));
x(3,:) = rotate(theta,[0,0,0],[0,0,1],[D(1,3),0,0])';
end

function xi = x_next(D,b,x,i)
A = x(i-3,:);
B = x(i-2,:);
C = x(i-1,:);
dA = D(i-3,i);
dB = D(i-2,i);
dC = D(i-1,i);
xs = intersect3spheres(A,dA,B,dB,C,dC);
xi = xs(b,:);
end

function errx = check_solution(L, U, x, i, xtol)
% considering only the precedent atoms
[~,~,l] = find(L(i,1:i));
[~,j,u] = find(U(i,1:i));
xi  = x(i,:);
errx = 0;
for k = 1:length(j)
    dk = norm(x(j(k),:) - xi);
    lk = l(k);
    uk = u(k);
    if lk > dk
        errx = errx + (lk - dk) / lk;
    end
    if uk < dk
        errx = errx + (dk - uk) / uk;
    end
    if errx > xtol
        return
    end
end
end