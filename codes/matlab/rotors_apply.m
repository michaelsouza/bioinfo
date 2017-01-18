function [x,g] = rotors_apply(x,S,W)
% x: coordinates' matrices
% s: index set
% w: angles
    if ~issorted(S)
        error('The index set s must be sorted in ascending order.')
    end
    % convert to matrix
    [nrows, ncols] = size(x);
    if nrows == 3 && ncols > nrows
        x = reshape(x,3,[]);
    end
    nnodes = length(x);
    % jacobian matrix
    J = zeros(3,nnodes,10,length(S));
	
    % applying rotations
    for k = 1:length(S)
        j = S(k);
        w = W(k);
        if j < 4
            continue;
        end
        % update x
        for q = j:nnodes
            [x(:,q),J(:,q,:,k)] = rotate(w,x(:,q),x(:,j-1),x(:,j-2));
        end
    end
    
    % calculating jacobian
    if nargout > 1
        g = zeros(3,nnodes,length(S));
        for k = 1:length(S)
            i = S(k);
            for q = i:length(x)
                % largest rotation index that affects Xq
                j = max(S(S<=q));
                g(:,q,k) = diff_x(q,j,i,S,J);
            end
        end
        
        % convert 3D matrix to jacobian
        % g = d x_i / d w_j
        g = reshape(g,3*nnodes,[]);
    end
end

 % d Xjq / d Wi
 function g = diff_x(q,j,i,S,J)
 j = min(j, max(S(S<=q)));
 [~,k] = ismember(j,S);
 if q < i || j < i || j > q
     g = [0,0,0];
 elseif i == j
     Jqk = reshape(J(:,q,:,k), [3, 10]);
     u = [1;0;0;0;0;0;0;0;0;0];
     g = (Jqk * u)';
 else
     Jqk = reshape(J(:,q,:,k),[3,10]);
     jnew = max(S(S<j));
     dp = diff_x(  q,jnew,i,S,J);
     dy = diff_x(j-1,jnew,i,S,J);
     dz = diff_x(j-2,jnew,i,S,J);
     u = [0,dp,dy,dz];
     g = (Jqk * u')';
 end
 end
 
 function [b,k] = ismember(j,S)
 if isempty(j)
     b = false;
     k = nan;
 else
     for k = 1:length(S)
         if S(k) == j
             b = true;
             return
         end
     end
 end
 end