% recoveryStressStrain Recovers stress and strain values at Gauss points for a 2D element.
%
% This function computes the total strain, thermal strain, elastic strain, and 
% thermal stress at Gauss points of a finite element using the provided nodal 
% coordinates, displacement vector, and temperature values. It also accounts 
% for temperature-dependent material properties.
%
% Syntax:
%   [totalStrain, thermalStrain, elasticStrain, thermalStress] = recoveryStressStrain(nodeCoordinates, Ue, Te, E, nu, alpha, temperatureDependent)
%
% Inputs:
%   nodeCoordinates     - Coordinates of the element nodes (4 x 2 matrix for a quadrilateral element).
%   Ue                  - Nodal displacement vector (8 x 1 for 4 nodes, with 2 DOF per node).
%   Te                  - Nodal temperature values (4 x 1 vector for temperature at each node).
%   E                   - Young's modulus.
%   nu                  - Poisson's ratio.
%   alpha               - Coefficient of thermal expansion.
%   temperatureDependent - Flag ('yes'/'no') for whether material properties are temperature-dependent.
%
% Outputs:
%   totalStrain         - Total strain at Gauss points (3 x 4 matrix, each column for one Gauss point).
%   thermalStrain       - Thermal strain at Gauss points (3 x 4 matrix).
%   elasticStrain       - Elastic strain at Gauss points (3 x 4 matrix).
%   thermalStress       - Thermal stress at Gauss points (3 x 4 matrix).

function [totalStrain,thermalStrain,elasticStrain,thermalStress] = recoveryStressStrain(nodeCoordinates,Ue,Te,E,nu,alpha,temperatureDependent)

% Initialize output matrices
totalStrain = zeros(3,4);
thermalStrain = zeros(3,4);
elasticStrain = zeros(3,4);
thermalStress = zeros(3,4);

% Get Gauss weights and points for numerical integration
[gaussWeights,gaussPoints] = gaussPointsWeights();

% Loop through each Gauss point
for i = 1:gaussWeights
    xi = gaussPoints(i,1);
    eta = gaussPoints(i,2);

    [shape,shapeDerivates]=shapefunctions(xi,eta);
    
    jac = shapeDerivates*nodeCoordinates;
    
    dN_dxbar  = Bmatrix(jac,shapeDerivates);

    % Construct the B matrix using vectorization
    B = [dN_dxbar(1,1), 0,        dN_dxbar(1,2), 0,        dN_dxbar(1,3), 0,        dN_dxbar(1,4), 0;
        0,        dN_dxbar(2,1),  0,        dN_dxbar(2,2), 0,        dN_dxbar(2,3), 0,        dN_dxbar(2,4);
        dN_dxbar(2,1), dN_dxbar(1,1),  dN_dxbar(2,2), dN_dxbar(1,2), dN_dxbar(2,3), dN_dxbar(1,3), dN_dxbar(2,4), dN_dxbar(1,4)];

    [Dtemp,~,alphaTemp,~] = materialMechfem(E,nu,shape,Te,temperatureDependent);
    
    if strcmp(temperatureDependent,'no')
        alphaTemp = alpha;
    end
    
    % Calculate total strain at Gauss point
    totalStrain(:,i) =  B*Ue;

    % Calculate thermal strain at Gauss point
    thermalStrain(:,i) =  [1;1;0]*alphaTemp*shape*Te;
    
    % Calculate elastic strain (total strain - thermal strain)
    elasticStrain(:,i) = totalStrain(:,i) - thermalStrain(:,i);
    
     % Calculate thermal stress using the elastic strain and material stiffness matrix
    thermalStress(:,i) = Dtemp*(totalStrain(:,i) - thermalStrain(:,i));

end