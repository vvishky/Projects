% numericalTesting Performs numerical tests on shape functions and their derivatives.
%
% This function carries out several numerical tests to verify the properties of 
% shape functions and their derivatives in the finite element method. It checks:
%   - The unity test for shape functions (sum of shape functions should equal 1).
%   - The sum of derivatives of shape functions with respect to local variables.
%   - The sum of derivatives of shape functions with respect to physical coordinates.
%
% Syntax:
%   numericalTesting(N, deriN, B, xi, eta)
%
% Inputs:
%   N      - A vector containing the shape function values at a given Gauss point.
%   deriN  - A matrix containing the derivatives of the shape functions 
%            with respect to the local variables (ξ, η).
%   B      - The B-matrix, which contains the derivatives of the shape functions 
%            with respect to the physical coordinates (x, y).
%   xi     - The Gauss point coordinate in the ξ direction.
%   eta    - The Gauss point coordinate in the η direction.
%
% Outputs:
%   None. The function prints results of the numerical tests to the command window.

% Define nodal temperatures T1, T2, T3, T4
% T = [0; 48; 60; 12];
% 
% %% 1. Analytical Calculation
% % Based on the assumed linear temperature field T(x, y) = a0 + a1*x + a2*y
% % The linear system of equations is:
% % T1 = a0 + a1*x1 + a2*y1
% % T2 = a0 + a1*x2 + a2*y2
% % T3 = a0 + a1*x3 + a2*y3
% % T4 = a0 + a1*x4 + a2*y4
% 
% % Set up the system of equations
% A = [1 nodalCoords(1,1) nodalCoords(1,2);
%     1 nodalCoords(2,1) nodalCoords(2,2);
%     1 nodalCoords(3,1) nodalCoords(3,2);
%     1 nodalCoords(4,1) nodalCoords(4,2)];
% 
% % Solve for the coefficients [a0; a1; a2]
%  coefficients = A\T;
% 
% % Coefficients a0, a1, a2 represent the temperature field: T(x,y) = a0 + a1*x + a2*y
% a1 = coefficients(2);
% a2 = coefficients(3);
% a1 = 48
% a2 = 12

function numericalTesting(N, deriN, B, xi, eta)


% Sum of shape functions
sumN = sum(N); 

% Sum of derivatives of shape functions with respect to local variables
sumDerivN = sum(deriN, 2);

% Sum of derivatives of shape functions with respect to physical coordinates
sumB = sum(B, 2);

% Print the sum of shape functions at the current Gauss point
fprintf('Sum of shape functions at (xi, eta) = (%f, %f): %f\n', xi, eta, sumN);

% Unity test for shape functions (should sum to 1)
if abs(sumN - 1) < 1e-12
    disp('Unity test passed.');
else
    disp('Unity test failed.');
end

% Test for the sum of derivatives with respect to local variables (should sum to 0)
if (sumDerivN(1,1) == 0) && (sumDerivN(2,1) == 0)
    disp('Sum of derivatives of the shape functions with respect to local variable test passed.');
else
    disp('Sum of derivatives of the shape functions with respect to local variable test failed.');
end

% Test for the sum of derivatives with respect to physical coordinates (should sum to 0)
if (sumB(1,1) == 0) && (sumB(2,1) == 0)
    disp('Sum of derivatives of the shape functions with respect to physical coordinates test passed.');
else
    disp('Sum of derivatives of the shape functions with respect to physical coordinates test failed.');
end

end

