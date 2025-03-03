% Bmatrix computes the B matrix for finite element analysis.
%
% This function calculates the B matrix, which is used to relate the 
% strains in an element to the displacements. The B matrix is computed 
% as the product of the inverse of the Jacobian matrix (jac) and the 
% shape function derivatives (shapeDerivates).
%
% Syntax:
%   B = Bmatrix(jac, shapeDerivates)
%
% Inputs:
%   jac            - A 2x2 Jacobian matrix, which represents the 
%                   transformation from the reference (local) 
%                   coordinates to the global coordinates. The 
%                   Jacobian should be computed based on the element 
%                   geometry and the mapping of local coordinates to 
%                   global coordinates.
%
%   shapeDerivates - A matrix of shape function derivatives with respect 
%                   to local coordinates. The size of this matrix should 
%                   be consistent with the number of shape functions 
%                   and the number of dimensions of the element.
%
% Outputs:
%   B              - The computed B matrix, which relates the 
%                   strain-displacement relationship in the finite 
%                   element method. The size of this matrix depends 
%                   on the number of shape functions and the dimension 
%                   of the problem.
%

function [B] = Bmatrix(jac, shapeDerivates)
% Calculate the B matrix using the Jacobian and shape derivatives
B = jac \ shapeDerivates; % Use matrix division to solve for B
end
