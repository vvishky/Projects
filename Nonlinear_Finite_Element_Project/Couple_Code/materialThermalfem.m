% materialThermalfem computes the thermal properties of the material based on
% the shape functions and previous time step temperature.
%
% This function calculates the thermal conductivity, density, and specific heat 
% of a material, along with their derivatives, using linear relationships with 
% respect to temperature. The properties are defined based on the local shape 
% functions and the temperature vector from the previous time step.
%
% The temperature dependent equations are taken from below references in the field of
% finite element analysis.
% ZhihuiLiu, ZhihuiLi , Qiang Ma. ”Nonlinear finite element algorithm for
% solving fully coupled thermomechanical problems under strong
% aerothermodynamic environment” 2023
%
% Syntax:
%   [thermalConductivity, thermalConductivityDeriv, rho, specificHeat, rhoderiv, specificHeatderiv] = materialThermalfem(shape, Te)
%
% Inputs:
%   shape    - Shape function values for the element.
%   Te       - Previous time step temperature vector.
%
% Outputs:
%   thermalConductivity      - Thermal conductivity matrix (2x2).
%   thermalConductivityDeriv - Derivative of thermal conductivity matrix (2x2).
%   rho                      - Density of the material.
%   specificHeat             - Specific heat of the material.
%   rhoderiv                 - Derivative of density.
%   specificHeatderiv        - Derivative of specific heat.

function [thermalConductivity, thermalConductivityDeriv, rho, specificHeat, rhoderiv, specificHeatderiv] = materialThermalfem(shape, Te)

% Calculate the temperature at the Gauss point using shape functions
temp = shape * Te;

% Define thermal conductivity as a function of temperature
k = 400 + 0.04 * temp; % Linear relation with temperature
thermalConductivity = [k, 0; 0, k]; % Construct 2x2 thermal conductivity matrix

% Define the derivative of thermal conductivity
kder = 0.04; % Derivative is constant
thermalConductivityDeriv = [kder, 0; 0, kder]; % Construct 2x2 derivative matrix

% Define density as a function of temperature
rho = 9000 + 0.009 * temp; % Linear relation with temperature

% Define specific heat as a function of temperature
specificHeat = 300 + 0.03 * temp; % Linear relation with temperature

% Define the derivatives of density and specific heat
rhoderiv = 0.009; % Derivative of density
specificHeatderiv = 0.03; % Derivative of specific heat

end
