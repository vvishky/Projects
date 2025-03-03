% shapefunctions computes the shape functions and their derivatives
% for a 2D quadrilateral element in finite element analysis.
%
% This function evaluates the shape functions based on the local
% coordinates (xi, eta) of the element. The shape functions are
% used to interpolate the field variables over the element domain.
%
% Syntax:
%   [shape, shapeDerivates] = shapefunctions(xi, eta)
%
% Inputs:
%   xi             - Local coordinate in the xi direction, typically
%                    ranging from -1 to 1 for a 2D quadrilateral element.
%
%   eta            - Local coordinate in the eta direction, typically
%                    ranging from -1 to 1 for a 2D quadrilateral element.
%
% Outputs:
%   shape          - A 1x4 vector containing the values of the shape
%                    functions at the specified (xi, eta) coordinates.
%                    For a quadrilateral element, the shape functions
%                    are defined as:
%                    N1 = 0.25 * (1 - xi) * (1 - eta)
%                    N2 = 0.25 * (1 + xi) * (1 - eta)
%                    N3 = 0.25 * (1 + xi) * (1 + eta)
%                    N4 = 0.25 * (1 - xi) * (1 + eta)
%
%   shapeDerivates - A 2x4 matrix containing the derivatives of the
%                    shape functions with respect to xi and eta.
%                    The first row corresponds to derivatives with
%                    respect to xi, and the second row corresponds
%                    to derivatives with respect to eta.
%

function [shape, shapeDerivates] = shapefunctions(xi, eta)
% Compute shape functions for the 2D quadrilateral element
shape(1) = 0.25 * (1 - xi) * (1 - eta);  % Shape function N1
shape(2) = 0.25 * (1 + xi) * (1 - eta);  % Shape function N2
shape(3) = 0.25 * (1 + xi) * (1 + eta);  % Shape function N3
shape(4) = 0.25 * (1 - xi) * (1 + eta);  % Shape function N4

% Initialize the shape derivatives matrix
shapeDerivates = zeros(2, 4);

% Compute the derivatives of the shape functions with respect to xi
shapeDerivates(1, 1) = 0.25 * (-1 + eta);   % dN1/dxi
shapeDerivates(1, 2) = 0.25 * (1 - eta);    % dN2/dxi
shapeDerivates(1, 3) = 0.25 * (eta + 1);    % dN3/dxi
shapeDerivates(1, 4) = 0.25 * (-eta - 1);   % dN4/dxi

% Compute the derivatives of the shape functions with respect to eta
shapeDerivates(2, 1) = 0.25 * (xi - 1);      % dN1/deta
shapeDerivates(2, 2) = 0.25 * (-xi - 1);     % dN2/deta
shapeDerivates(2, 3) = 0.25 * (xi + 1);      % dN3/deta
shapeDerivates(2, 4) = 0.25 * (1 - xi);      % dN4/deta
end
