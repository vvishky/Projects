% gaussPointsWeights computes the Gauss quadrature points and weights.
%
% This function returns the Gauss quadrature weights and points for 
% a two-dimensional quadrilateral element. The quadrature points are 
% used for numerical integration in finite element analysis, particularly 
% for evaluating integrals over the element's domain.
%
% Syntax:
%   [gaussWeights, gaussPoints] = gaussPointsWeights()
%
% Outputs:
%   gaussWeights   - A scalar value representing the weight for each 
%                    Gauss quadrature point. In this case, the value 
%                    is set to 4, which corresponds to the weight used 
%                    in 2D quadrature for a unit square.
%
%   gaussPoints    - A 4x2 matrix containing the coordinates of the 
%                    Gauss quadrature points in the reference element. 
%                    These points are given in the local coordinates 
%                    of the element and are used to evaluate integrals.

function [gaussWeights, gaussPoints] = gaussPointsWeights()

% Set Gauss quadrature weight for the points
gaussWeights = 4; % Weight for each Gauss point in the quadrature rule

% Define the Gauss quadrature points in 2D
gaussPoints = [-sqrt(1/3), -sqrt(1/3);  % Point 1
               sqrt(1/3), -sqrt(1/3);  % Point 2
               sqrt(1/3),  sqrt(1/3);  % Point 3
              -sqrt(1/3),  sqrt(1/3)]; % Point 4
end
