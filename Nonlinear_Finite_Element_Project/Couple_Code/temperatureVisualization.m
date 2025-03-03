% temperatureVisualization Visualizes the 2D temperature distribution using a triangular surface plot.
%
% This function takes the coordinates of the nodes, the element connectivity, 
% and the global temperature values to create a surface plot representing 
% the temperature distribution over a 2D mesh.
%
% Syntax:
%   temperatureVisualization(x, y, elementsConnectivity, TGlobal)
%
% Inputs:
%   x                - A vector containing the x-coordinates of the nodes.
%   y                - A vector containing the y-coordinates of the nodes.
%   elementsConnectivity - An Mx3 matrix defining the connectivity of triangular elements,
%                         where M is the number of elements. Each row corresponds to
%                         an element defined by the node indices.
%   TGlobal          - A vector containing the temperature values at each node.
%
% Outputs:
%   None. The function produces a 3D surface plot of the temperature distribution.

function temperatureVisualization(x, y, elementsConnectivity, TGlobal)

% Create a triangular surface plot of the temperature distribution
trisurf(elementsConnectivity, x, y, TGlobal, 'EdgeColor', 'k');

% Set the view to 2D
view(2);

% Create a custom colormap with a gradient from blue to cyan to yellow
colormap turbo;

% Add a colorbar to indicate temperature scale
colorbar;

% Add labels and title to the plot
title('2D Temperature Distribution');
xlabel('X Coordinate (m)');
ylabel('Y Coordinate(m)');
zlabel('Temperature (°C)');

% Apply shading interpolation for a smoother appearance
shading interp;

end
