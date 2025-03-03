% plotMesh Visualizes a 2D mesh grid using the patch function.
%
% This function takes in nodal coordinates and element connectivity 
% to create a visual representation of the mesh. The mesh is plotted 
% as a collection of polygons defined by the vertices and faces.
%
% Syntax:
%   plotMesh(nodes, elements)
%
% Inputs:
%   nodes    - An Nx2 matrix containing the coordinates of the nodes,
%              where N is the number of nodes. Each row represents
%              a node's (x, y) coordinates.
%   elements  - An MxK matrix defining the connectivity of the elements,
%               where M is the number of elements and K is the number
%               of nodes per element. Each row corresponds to an element
%               defined by the node indices.
%
% Outputs:
%   None. The function produces a plot of the mesh grid.

function plotMesh(nodes, elements)
% Create a figure for the mesh plot
figure(1);

% Plot the mesh using the patch function
patch('Faces', elements, 'Vertices', nodes, 'FaceColor', 'white');

% Set equal scaling for the axes
axis equal;

% Add title and axis labels
title('2D Mesh Grid');
xlabel('X Coordinate (m)');
ylabel('Y Coordinate (m)');

end
