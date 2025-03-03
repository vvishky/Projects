% elementMechThermfem computes the mechanical and thermal stiffness matrix 
% and the force vector for a 2D finite element analysis.
%
% This function evaluates the element stiffness matrix and the residual 
% vector considering both mechanical and thermal effects. The calculation 
% is performed using Gauss quadrature for numerical integration. The 
% material properties may depend on the temperature.
%
% Syntax:
%   [K_e_mt, R_e_mt] = elementMechThermfem(nodeCoordinates, Te, E, nu, alpha, temperatureDependent)
%
% Inputs:
%   nodeCoordinates     - Coordinates of the nodes in the element 
%                         (8x2 matrix).
%
%   Te                  - Reference temperature (scalar).
%
%   E                   - Young's modulus of the material (scalar).
%
%   nu                  - Poisson's ratio of the material (scalar).
%
%   alpha               - Coefficient of thermal expansion (scalar).
%
%   temperatureDependent - A string indicating if the material properties 
%                         are temperature-dependent ('yes' or 'no').
%
% Outputs:
%   K_e_mt             - The element stiffness matrix considering mechanical 
%                        and thermal effects (8x4 matrix).
%
%   R_e_mt             - The element residual vector (8x1 matrix).

function [K_e_mt, R_e_mt] = elementMechThermfem(nodeCoordinates, Te, E, nu, alpha, temperatureDependent)

% Initialize the element stiffness matrix and residual vector
K_e_mt = zeros(8, 4);
R_e_mt = zeros(8, 1);

% Get Gauss weights and points for integration
[gaussWeights, gaussPoints] = gaussPointsWeights();

% Loop over each Gauss point for numerical integration
for i = 1:gaussWeights
    xi = gaussPoints(i, 1);
    eta = gaussPoints(i, 2);

    % Evaluate shape functions and their derivatives at the Gauss point
    [shape, shapeDerivates] = shapefunctions(xi, eta);

    % Compute the Jacobian and its determinant
    jac = shapeDerivates * nodeCoordinates;
    detJ = det(jac);

    % Calculate the derivatives of shape functions with respect to global coordinates
    dN_dxbar = Bmatrix(jac, shapeDerivates);

    % Construct the B matrix using vectorization
    B = [dN_dxbar(1, 1), 0,        dN_dxbar(1, 2), 0,        dN_dxbar(1, 3), 0,        dN_dxbar(1, 4), 0;
         0,        dN_dxbar(2, 1), 0,        dN_dxbar(2, 2), 0,        dN_dxbar(2, 3), 0,        dN_dxbar(2, 4);
         dN_dxbar(2, 1), dN_dxbar(1, 1), dN_dxbar(2, 2), dN_dxbar(1, 2), dN_dxbar(2, 3), dN_dxbar(1, 3), dN_dxbar(2, 4), dN_dxbar(1, 4)];
    
    % Get material properties for the current Gauss point
    [Dtemp, Dderiv, alphaTemp, alphaDeri] = materialMechfem(E, nu, shape, Te, temperatureDependent);
    
    % Use provided alpha if material properties are not temperature-dependent
    if strcmp(temperatureDependent, 'no')
        alphaTemp = alpha;
    end
    
    % Compute the element stiffness matrix and residual vector
    K_e_mt = K_e_mt + (B' * Dtemp * [1; 1; 0] * alphaTemp * shape + ...
                        B' * Dtemp * [1; 1; 0] * alphaDeri * shape * Te * shape + ...
                        B' * Dderiv * [1; 1; 0] * alphaTemp * shape * Te * shape) * detJ;

    R_e_mt = R_e_mt + B' * Dtemp * [1; 1; 0] * alphaTemp * shape * Te * detJ;

end

end
