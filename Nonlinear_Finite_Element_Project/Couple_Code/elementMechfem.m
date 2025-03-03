% elementMechfem computes the element stiffness matrix and the
% residual force vector for a 2D mechanical finite element analysis.
%
% This function evaluates the stiffness matrix and the residual force
% vector for a 2D quadrilateral element using the finite element method.
% It utilizes Gauss quadrature for numerical integration and accounts
% for temperature-dependent material properties.
%
% Syntax:
%   [K_e_m, R_e_m] = elementMechfem(nodeCoordinates, U_e, Te, E, nu, temperatureDependent)
%
% Inputs:
%   nodeCoordinates     - A 4x2 matrix containing the coordinates of
%                         the element nodes in global coordinates. Each
%                         row corresponds to a node, with the first column
%                         representing the x-coordinate and the second column
%                         representing the y-coordinate.
%
%   U_e                 - A 8x1 vector representing the nodal displacements
%                         of the element in the global coordinate system.
%
%   Te                  - A scalar value representing the reference temperature
%                         for the material.
%
%   E                   - Young's modulus of the material.
%
%   nu                  - Poisson's ratio of the material.
%
%   temperatureDependent - A boolean value indicating whether the material
%                         properties are temperature-dependent (true) or
%                         constant (false).
%
% Outputs:
%   K_e_m               - The 8x8 element stiffness matrix for the
%                         mechanical response of the element.
%
%   R_e_m               - The 8x1 residual force vector for the element.
%
% This code is adapted from the work of Srinivas Kokula.
% The Plane Stress Problem from Siva Srinivas Kolukula. https://de.mathworks.com/matlabcentral/fileexchange/31788-the-plane-stress-problem

function [K_e_m, R_e_m] = elementMechfem(nodeCoordinates, U_e, Te, E, nu, temperatureDependent)
% Initialize the element stiffness matrix and residual vector
K_e_m = zeros(8, 8);
R_e_m = zeros(8, 1);

% Get Gauss quadrature weights and points for integration
[gaussWeights, gaussPoints] = gaussPointsWeights();

% Loop over each Gauss point for numerical integration
for i = 1:gaussWeights
    xi = gaussPoints(i, 1);  % Gauss point coordinate in xi direction
    eta = gaussPoints(i, 2);  % Gauss point coordinate in eta direction

    % Evaluate shape functions and their derivatives at the Gauss point
    [shape, shapeDerivates] = shapefunctions(xi, eta);

    % Compute the Jacobian matrix and its determinant
    jac = shapeDerivates * nodeCoordinates;  % Jacobian matrix
    detJ = det(jac);  % Determinant of the Jacobian

    % Compute the derivatives of shape functions in the global coordinates
    dN_dxbar = Bmatrix(jac, shapeDerivates);

    % Construct the B matrix for strain-displacement relations using vectorization
    B = [dN_dxbar(1, 1), 0,        dN_dxbar(1, 2), 0,        dN_dxbar(1, 3), 0,        dN_dxbar(1, 4), 0;
        0,        dN_dxbar(2, 1), 0,        dN_dxbar(2, 2), 0,        dN_dxbar(2, 3), 0,        dN_dxbar(2, 4);
        dN_dxbar(2, 1), dN_dxbar(1, 1), dN_dxbar(2, 2), dN_dxbar(1, 2), dN_dxbar(2, 3), dN_dxbar(1, 3), dN_dxbar(2, 4), dN_dxbar(1, 4)];

    % Evaluate material stiffness matrix D for the current Gauss point
    [Dtemp, ~, ~, ~] = materialMechfem(E, nu, shape, Te, temperatureDependent);

    % Add stiffness matrix on gauss points
    K_e_m = K_e_m + B' * Dtemp * B * detJ;

    % Add the residual force vector on gauss points
    R_e_m = R_e_m + B' * Dtemp * B * U_e * detJ;
end
end
