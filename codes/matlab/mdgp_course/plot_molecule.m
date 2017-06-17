function ax = plot_molecule(x)

% cylinder radius
rc = 0.2; 

% resolution of spheres and cylinders
NS   = 50; % spheres, more looks smoother
NB   = 50; % more looks smoother

natoms = size(x,1); % number of atoms in molecule

% set atoms labels
symbols    = 'NCC'; % backbone pattern
lab        = repmat('H',1,natoms);
lab(2:end) = symbols(mod(0:(natoms-2),3)+1);

% prepare axis
ax = newplot([]);
axes(ax);
set(ax,'Visible','off','DataAspectRatio',[1 1 1]) % fix aspect ratio

light('Position',[1 1 2]); % add some light

% draw spheres
for k = 1:natoms
    colk = col(lab(k));
    radk = rad(lab(k));
    
    % basic sphere
    [sx,sy,sz] = sphere(NS);
    
    % draw sphere
    surface('XData',x(k,1) + radk*sx,'YData',x(k,2) + radk*sy, ...
        'ZData',x(k,3) + radk*sz,'FaceColor',colk, ...
        'EdgeColor','none','FaceLighting','gouraud')
end

% draw cylinders for each bond
for k = 1:(natoms-1)
    xi = x(k  ,:); % coordinates atom i
    xj = x(k+1,:); % coordinates atom 2
    
    % bond angles in spherical coordinates
    v     = (xj - xi)/norm(xj - xi);
    phi   = atan2d(v(2),v(1));
    theta = -asind(v(3));
    
    % bond distance minus sphere radii
    ds   = norm(xi - xj);
    cyli = ds; % length full bond cylinder
    
    % prototype cylinders for bond
    [ciz,ciy,cix] = cylinder(rc,NB); % full bond cylinder
    cix(2,:)      = cix(2,:) * cyli; % adjust length
    
    % rotate cylinders to match bond vector v
    for kk = 1:numel(cix)
        vr = [cix(kk); ciy(kk); ciz(kk);];
        vr = rotz(phi)*roty(theta)*vr;
        cix(kk) = vr(1);
        ciy(kk) = vr(2);
        ciz(kk) = vr(3);
    end
    
    % get colors of both atoms
    coli = col(lab(k));
    
    % full bond color 1
    surface('XData',xi(1) + cix,'YData',xi(2) + ciy,...
        'ZData',xi(3) + ciz,'FaceColor',coli,...
        'EdgeColor','none','FaceLighting','gouraud')
end

%---------------------------------------------
% element specific CPK colors
function c = col(s)
switch s
    case  'H', c = [1 1 1];
    case  'C', c = [0.2 0.2 0.2];
    case  'N', c = [0.5 0.5 0.6];
end

%---------------------------------------------
% element specific radii
function r = rad(s)
switch s
    case  'H', r = 0.4;
    case  'C', r = 0.5;
    case  'N', r = 0.5;
 end