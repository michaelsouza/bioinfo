function [x,p,b] = solver_ibp(problem,np)
% np: number of partitions

% start branches
p = ones(1, problem.natoms);
p(1:3) = np;
b = ones(1, problem.natoms);
b(1:3) = 2;

i = 4;
while i > 3
    b(i) = 1;
    % set and solve subproblem
    subproblem = set_subproblem(problem, np, p, i);
    [x,bi] = solver_bp(subproblem, b(1:i));
    
    prune = isempty(x);
    if ~prune
        if i == problem.natoms
            fprintf('p = '); fprintf('%g ', p(1:i)); fprintf('\n');
            fprintf('IBP:Solution found\n');
            return
        else
            % update start branch
            b(1:i) = bi;
        end
    else
        % reset branch
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
    fprintf('IBP:No solution\n');
end
end

function subproblem = set_subproblem(problem, np, p, i)
index = problem.i <= i & problem.j <= i;
subproblem.i = problem.i(index);
subproblem.j = problem.j(index);
subproblem.d = problem.d(index);
subproblem.l = problem.l(index);
subproblem.u = problem.u(index);
subproblem.natoms = i;
subproblem.nedges = length(subproblem.i);
for k = 1:subproblem.nedges
    ik = subproblem.i(k);
    jk = subproblem.j(k); 
    if abs(ik - jk) == 3 
        if ik > jk
            jk = ik;
        end
    end
    lk = subproblem.l(k);
    uk = subproblem.u(k);
    subproblem.d(k) = lk + (uk - lk) * p(jk) / np;
end
end