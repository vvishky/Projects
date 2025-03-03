% materialMechfem computes the material properties matrix for a
% 2D mechanical finite element analysis.
%
% This function evaluates the constitutive matrix for plane stress
% conditions, considering both temperature-dependent and constant
% material properties. It computes the elastic modulus, its
% derivative with respect to temperature, and the thermal expansion
% coefficient and its derivative.
%
% The temperature dependent equations are taken from below references in the field of
% finite element analysis.
% ZhihuiLiu, ZhihuiLi , Qiang Ma. ”Nonlinear finite element algorithm for
% solving fully coupled thermomechanical problems under strong
% aerothermodynamic environment” 2023
%
% Syntax:
%   [Dtemp, Dderiv, alphaTemp, alphaDeri] = materialMechfem(E, nu, shape, Te, temperatureDependent)
%
% Inputs:
%   E                   - Young's modulus of the material (scalar).
%
%   nu                  - Poisson's ratio of the material (scalar).
%
%   shape               - A vector of shape functions evaluated at the
%                         reference temperature (1xN, where N is the
%                         number of shape functions).
%
%   Te                  - Reference temperature (scalar).
%
%   temperatureDependent - A string indicating if the material properties
%                         are temperature-dependent ('yes' or 'no').
%
% Outputs:
%   Dtemp               - The constitutive matrix for plane stress (3x3).
%
%   Dderiv              - The derivative of the constitutive matrix with
%                         respect to temperature (3x3).
%
%   alphaTemp           - The thermal expansion coefficient (scalar).
%
%   alphaDeri           - The derivative of the thermal expansion
%                         coefficient with respect to temperature (scalar).
%


function [Dtemp, Dderiv, alphaTemp, alphaDeri] = materialMechfem(E, nu, shape, Te, temperatureDependent)

% Check if the material properties are temperature-dependent
if strcmp(temperatureDependent, 'yes')

    % Calculate the temperature based on the shape functions and reference temperature
    temp = shape * Te;

    % Update Young's modulus based on temperature
    E = 2e7 - (2e4) * temp;

    % Derivative of Young's modulus with respect to temperature
    Ederiv = -2e4;

    % Constitutive matrix for plane stress conditions
    Dtemp = E / (1 - nu^2) * [1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];

    % Derivative of the constitutive matrix with respect to temperature
    Dderiv = Ederiv / (1 - nu^2) * [1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];

    % Calculate thermal expansion coefficient based on temperature
    alphaTemp = 0.001 - 0.000001 * temp;

    % Derivative of thermal expansion coefficient with respect to temperature
    alphaDeri = -0.000001;

else
    % Constitutive matrix for plane stress conditions with constant properties
    Dtemp = E / (1 - nu^2) * [1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];

    % Derivative of the constitutive matrix with respect to temperature
    Dderiv = 0;

    % Constant thermal expansion coefficient
    alphaTemp = 0;

    % Derivative of thermal expansion coefficient with respect to temperature
    alphaDeri = 0;
end

end
