%
% Examples for the usage of molecule3D.m
%
%   Version: 1.2
%   Author:  André Ludwig (aludwig@phys.ethz.ch)

%%
%
% Example 1: Simple call of the function with positions xyz and labels
% 

% geometry of ethylene C2H4 in Angström
xyz = [0.000    0.000    0.000;  % C
       1.340    0.000    0.000;  % C
      -0.545    0.944    0.000;  % H
      -0.545   -0.944    0.000;  % H
       1.885    0.944    0.000;  % H
       1.885   -0.944    0.000]; % H
   
% labels\elements of the xyz columns
labels = {'C' 'C' 'H' 'H' 'H' 'H'};

figure
set(gcf,'Color','w')
molecule3D(xyz,labels)
cameratoolbar

%%
%
% Example 2: Comparison of styles
% 

% geometry of ethylene C2H4 in Angström
xyz = [0.000    0.000    0.000;  % C
       1.340    0.000    0.000;  % C
      -0.545    0.944    0.000;  % H
      -0.545   -0.944    0.000;  % H
       1.885    0.944    0.000;  % H
       1.885   -0.944    0.000]; % H
   
% labels\elements of the xyz columns
labels = {'C' 'C' 'H' 'H' 'H' 'H'};

styles = {'ballstick','licorice','large','superlarge'};
figure
set(gcf,'Color','w')
for k = 1:numel(styles)
    subplot(2,2,k)
    molecule3D(xyz,labels,styles{k})
    text(0.67,1.5+k/4,0,styles{k},'HorizontalAlignment','center')
end

%%
%
% Example 3: Simple call of the function for water
% 

% geometry of water H2O in Angström
xyz = [ 0.000    0.000    0.000;  % O
        0.762   -0.584    0.000;  % H
       -0.762   -0.584    0.000]; % H
   
% labels\elements of the xyz columns
labels = {'O' 'H' 'H'};

figure
set(gcf,'Color','w')
molecule3D(xyz,labels)

%%
%
% Example 4: Simple call of the function for iodo acetylene
% 

% geometry of iodo acetylene ICCH in Angström
xyz = [-1.590    0.000    0.000;  % I
        0.000    0.000    0.000;  % C
        1.200    0.000    0.000;  % C
        2.260    0.000    0.000]; % H
    
% labels\elements of the xyz columns
labels = {'I' 'C' 'C' 'H'};

figure
set(gcf,'Color','w')
molecule3D(xyz,labels)