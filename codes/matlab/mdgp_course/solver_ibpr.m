function [x,p,b] = solver_ibpr(problem, np)
% set distance (bound) matrix
D = sparse(problem.i, problem.j, problem.d, problem.natoms, problem.natoms);
D = D + D';
L = sparse(problem.i, problem.j, problem.l, problem.natoms, problem.natoms);
L = L + L';
U = sparse(problem.i, problem.j, problem.u, problem.natoms, problem.natoms);
U = U + U';

% set initial point
x = x_start(D);

% start branches
p = ones(1, problem.natoms);
p(1:3) = np;
b = ones(1, problem.natoms);
b(1:3) = 2;
i_start = 4;

i = 4;
while i > 3
    [x,bi] = bpr(x,i_start,p,np,b(1:i),D,L,U);
    prune  = isempty(x);
    
    if ~prune
        if i == problem.natoms
            fprintf('p = '); fprintf('%g ', p(1:i)); fprintf('\n');
            fprintf('IBP:Solution found\n');
            return
        else
            % update start branch
            b(1:i)  = bi;
            i_start = i + 1;
        end
    else
        % reset ibp
        i_start  = 4;
        b(4:end) = 1;
        fprintf('IBP:prune\n')
    end
    
    prune = prune || i == problem.natoms;
    % backward
    if prune && p(i) == np
        while i > 3 && p(i) == np
            i = i - 1;
        end
    end
    
    % reset subtree
    if prune && p(i) < np
        p(i) = p(i) + 1;
        for j = (i+1):problem.natoms
            p(j) = 1;
        end
    end
    
    if ~prune
        i = i + 1;
    end
end
if isempty(x)
    b = [];
    fprintf('IBPR:No solution\n');
end
end

function lambda = get_lambda(u,v,p,y)
w = v - u;
w = w / norm(w);
k = dot(w, u - p);
z = u(1)*(w(2)^2+w(3)^2)-w(1)*(k-u(1)*w(1));
A = p(1) - z;
B = -u(3)*w(2)+u(2)*w(3)-w(3)*p(2)+w(2)*p(3);
C = z-y(1);
D = sqrt((A*B)^2+B^4-(B*C)^2);
E = A^2+B^2;

if abs(A) >= 1E-8 && abs(B) <= 1E-8
    lambda = acos(-C/A);
end
if abs(A) < 1E-8 && abs(B) >= 1E-8
    lambda = asin(-C/B);
end
if abs(A) >= 1E-8 && abs(B) >= 1E-8
    lambda1 = atan2((-C+(A*(A*C+D))/E)/B,(-A*C-D)/E);
    lambda2 = atan2((-C+(A*(A*C-D))/E)/B,(-A*C+D)/E);
    y1 = rotate(lambda1,u,v,p)';
    y2 = rotate(lambda2,u,v,p)';
    if norm(y1 - y) < norm(y2 - y)
        lambda = lambda1;
    else
        lambda = lambda2;
    end
    if abs(A) < 1E-8 && abs(B) < 1E-8
        error('IBPR::Lambda could not be calculated\n');
    end
end
end

function [x,b] = bpr(x,i,p,np,b,D,L,U)
xtol = 1E-5;
n = length(b);
while i > 3
    % reference point (lambda = 0)
    D(i-3,i) = L(i-3,i);
    P = x_next(D,b(i),x,i,i-3);
    
    % arc
    [phi,flip] = get_phi(x,i,b,P,D,L,U);
    
    prune = phi(1) > phi(2); % empty intersection
    if ~prune
        lambda = phi(1) + (phi(2)-phi(1)) / np * p(i); % rotation angle
        x(i,:) = rotate((-1)^flip * lambda,x(i-2,:),x(i-1,:),P);
        errx   = check_solution(L,U,x,i,xtol);
        prune  = errx > xtol;
    end
    
    if ~prune && i == n
        fprintf('BPR:Solution found, natoms: %d\n',n);
        return
    end
    
    prune = prune || i == n;
    % backward
    if prune && b(i) == 2
        while i > 3 && b(i) == 2
            i = i - 1;
        end
    end
    
    % reset subtree
    if prune && b(i) == 1
        b(i) = 2;
        for j = (i+1):n
            b(j) = 1;
        end
    end
    
    if ~prune
        i = i + 1;
    end
end
x = [];
b = [];
fprintf('BPR:No solution\n');
end

function [phi,flip] = get_phi(x,i,b,P,D,L,U)
[~,j,~] = find(L(i,1:i-3));
u = x(i-2,:); v = x(i-1,:);
lambda = zeros(length(j),2);
for k = 1:length(j)
    jk = j(k);
    D(jk,i) = L(jk,i);
    A = x_next(D,b(i),x,i,jk);
    lambda(k,1) = get_lambda(u,v,P,A);
    D(jk,i) = U(jk,i);
    B = x_next(D,b(i),x,i,jk);
    lambda(k,2) = get_lambda(u,v,P,B);
end
flip = lambda(end,1) > lambda(end,2);
if flip
    lambda = -lambda;
end
phi = [max(lambda(:,1)), min(lambda(:,2))];
end

function x = x_start(D)
x = zeros(length(D), 3);
x(2,1) = -D(1,2);
theta  = acos((D(1,2)^2 + D(1,3)^2 - D(2,3)^2) / (2 * D(1,2) * D(1,3)));
x(3,:) = rotate(-theta,[0,0,0],[0,0,1],[-D(1,3),0,0])';
end

function xi = x_next(D,b,x,i,j)
A = x(j,:);
B = x(i-2,:);
C = x(i-1,:);
dA = D(j,i);
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