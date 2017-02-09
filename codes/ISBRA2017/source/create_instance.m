function create_instance(nins, np, dmax)
% Parameters
% nins : number of instances
% np   : number of points. Instance size
% dmax : max distance. All distances dij = ||xi-xj|| > dmax for |i-j| > 3 are discarded 
for k = 1:nins
    fprintf('Creating instance (%3d/%3d)\n', k, nins);
    eps_dij   = 0.1;
    eps_omega = 0.05;
    x = generate_random_backbone(np, eps_omega);
    fprintf('   Calculating bounds\n');
    I = zeros(np^2,1);
    J = zeros(np^2,1);
    L = zeros(np^2,1);
    U = zeros(np^2,1);
    nedges = 0;
    for i = 1:np
        for j = (i+1):np
            dij = norm(x(i,:) - x(j,:));
            if dij < dmax && (j - i) >= 3
                % imprecise distance
                nedges = nedges + 1;
                I(nedges) = i;
                J(nedges) = j;
                L(nedges) = (1 - eps_dij) * dij;
                U(nedges) = (1 + eps_dij) * dij;
            elseif (j - i) < 3
                % exact distance
                nedges = nedges + 1;
                I(nedges) = i;
                J(nedges) = j;
                L(nedges) = dij;
                U(nedges) = dij;
            end
        end
    end
    % saving file
    import java.util.UUID;
    uid   = char(UUID.randomUUID());
    fname = sprintf('../instances/mdgp_N%d_D%3.2f_%s', np, dmax, uid);
    fprintf('   Saving files: %s\n', fname);
    table = struct('x',x(:,1),'y',x(:,2),'z',x(:,3));
    writetable(struct2table(table), [fname '_sol.csv']); % saving solution
    table = struct('I',I(1:nedges),'J',J(1:nedges),'L',L(1:nedges),'U',U(1:nedges));
    writetable(struct2table(table), [fname '_ins.csv']); % saving instance
end
end

function x = generate_random_backbone(np, eps_omega)
    fprintf('   Generating random backbone\n');
    % length of covalent bonds
    d = 1.526 * ones(np,1);
    d(1) = 0.0;
    % plane angles
    theta = 1.91 * ones(np,1);
    % dihedral angles with perturbation of 5% at most
    angles = [pi/3;pi/2;5*pi/3]; % typical angles
    omega = angles(randi(3,np,1)) .* random('Uniform',1-eps_omega,1+eps_omega,np,1);
    x = set_coords(d,theta,omega);
end

function x = set_coords(d,theta,omega)
    fprintf('      Set coords\n');
    % Calculates Euclidean coords from distances, planar and dihedral angles 
    % using transformation matrices;
    %   Input parameters;
    %   d    : distances       (d(i) = ||x(i,:)-x(i-1,:)||;            
    %   theta: planar angles   (theta(i) = angle(x(i,:),x(i-1,:),x(i-2,:));        
    %   omega: dihedral angles (omega(i) = angle(x(i,:),x(i-1,:),x(i-2,:),x(i-3,:));

    np = length(d);
    x = zeros(np,3);
    ct = cos(theta(3));
    st = sin(theta(3));
    x(1:3,:) = [0,0,0;-d(2),0,0;-d(2)+d(3)*ct,d(3)*st,0];

    di = d(3);
    B1 = eye(4);
    B2 = [-1,0,0,-d(2);0,1,0,0;0,0,-1,0;0,0,0,1];
    B3 = [-ct,-st,0,-di*ct;st,-ct,0,di*st;0,0,1,0;0,0,0,1];

    % set remaining coordinates
    B = B1 * B2 * B3;
    for i = 4:np
        ct = cos(theta(i)); st = sin(theta(i));
        cw = cos(omega(i)); sw = sin(omega(i));
        di = d(i);
        % concatenate transformation matrices
        B = B * [-ct,-st,0,-di*ct;st*cw,-ct*cw,-sw,di*st*cw;st*sw,-ct*sw,cw,di*st*sw;0,0,0,1];
        % apply transformation at origin
        v = B  * [0;0;0;1];
        % update coordinate
        x(i,:) = v(1:3);
    end
end