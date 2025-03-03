% elementThermalfem computes the element stiffness matrix and residual vector
% for thermal finite element analysis.
%
% This function calculates the stiffness matrix (K_e) and residual vector (R_e)
% based on the given nodal IDs, boundary conditions, material properties,
% and temperature distributions. It handles both implicit time integration
% and various boundary conditions including convective, heat flux, and 
% radiation.
%
% Syntax:
%   [K_e, R_e] = elementThermalfem(nodelIds, leftBoundaryNodeId, rightBoundaryNodeId, topBoundaryNodeId, nodeCoordinates, thermalConductivity, Te, testCase, ...
%                                   heatFlux, heatCoeff, ambientTemp, rho, specificHeat, deltaT, tempUpdate, temperatureDependent, emissivity, stefanBoltzmann, numericalTest)
%
% Inputs:
%   nodelIds             - IDs of the nodes associated with the element.
%   leftBoundaryNodeId   - Node IDs on the left boundary.
%   rightBoundaryNodeId  - Node IDs on the right boundary.
%   topBoundaryNodeId    - Node IDs on the top boundary.
%   nodeCoordinates       - Nx2 matrix of node coordinates (x,y).
%   thermalConductivity   - Thermal conductivity of the material.
%   Te                    - Previous time step temperature vector.
%   testCase             - Test case identifier for boundary conditions.
%   heatFlux             - Heat flux applied on the left boundary.
%   heatCoeff            - Convective heat transfer coefficient.
%   ambientTemp          - Ambient temperature for convection.
%   rho                  - Density of the material.
%   specificHeat         - Specific heat of the material.
%   deltaT               - Time step size for implicit integration.
%   tempUpdate           - Updated temperature from the previous step.
%   temperatureDependent  - Flag for temperature-dependent properties.
%   emissivity           - Emissivity for radiation boundary conditions.
%   stefanBoltzmann      - Stefan-Boltzmann constant for radiation.
%   numericalTest        - Flag to enable numerical testing outputs.
%
% Outputs:
%   K_e                  - Element stiffness matrix.
%   R_e                  - Residual vector.


function [K_e,R_e] = elementThermalfem(nodelIds,leftBoundaryNodeId,rightBoundaryNodeId,topBoundaryNodeId,nodeCoordinates,thermalConductivity,Te,testCase,...
    heatFlux,heatCoeff,ambientTemp,rho,specificHeat,deltaT,tempUpdate,temperatureDependent,emissivity,stefanBoltzmann,numericalTest)

% Initialize element stiffness matrix (K_e) and residual vector (R_e)
K_e = zeros(4, 4);
K_H = zeros(4, 4);   % Matrix for convective boundary condition
K_R = zeros(4, 4);   % Matrix for radiation boundary condition
qFLux = zeros(4, 1); % Heat flux vector
qAmbient = zeros(4, 1); % Ambient heat vector
qAmbientR = zeros(4, 1); % Ambient radiation vector
Q_int_e = zeros(4, 1); % Internal heat vector for elements
Q_ext_e = zeros(4, 1); % External heat vector for elements
Q_int_r = zeros(4, 1); % Internal heat vector for radiation

% Data for numerical testing
area = 0; % Area of the element
tempGradient = [0; 0]; % Temperature gradient initialization
T = [0; 48; 60; 12];

% Get Gauss points and weights for numerical integration
[gaussWeights, gaussPoints] = gaussPointsWeights();

% Loop over Gauss points for integration
for i = 1:gaussWeights

    xi = gaussPoints(i, 1); % xi coordinate
    eta = gaussPoints(i, 2); % eta coordinate

    % Calculate shape functions and their derivatives
    [shape, shapeDerivates] = shapefunctions(xi, eta);

    % Compute the Jacobian and its determinant
    jac = shapeDerivates * nodeCoordinates;
    detJ = det(jac);

    % Calculate the B matrix
    [B] = Bmatrix(jac, shapeDerivates);
    
    % Perform numerical testing if enabled
    if strcmp(numericalTest,'yes')

        numericalTesting(shape,shapeDerivates,B,xi, eta)

        % Accumulate area and temperature gradient
        area = area + detJ;
        tempGradient = tempGradient + B*T*detJ;
        continue;% Skip further calculations for this iteration
    end
    
    % Update thermal properties if they are temperature dependent
    if strcmp(temperatureDependent,'yes')
        [thermalConductivity,thermalConductivityDeriv,rho,specificHeat,rhoderiv,specificHeatderiv] = materialThermalfem(shape,Te);
    else
        thermalConductivityDeriv = [0,0;0,0];
        rhoderiv = 0;
        specificHeatderiv = 0;
    end

    %Implicit Time integration scheme.
    K_e = K_e + ((rho*specificHeat*(shape'*shape))/deltaT + ((rhoderiv*specificHeat*(shape'*shape))*(tempUpdate - Te)*shape)/deltaT + ...
        ((rho*specificHeatderiv*(shape'*shape))*(tempUpdate - Te)*shape)/deltaT+ B'*thermalConductivity*B+ B'*thermalConductivityDeriv*B*Te*shape)*detJ;
    
    % Accumulate internal heat vector
    Q_int_e = Q_int_e+ ((rho*specificHeat*(shape'*shape)*(tempUpdate - Te))/deltaT + B'*thermalConductivity*B*tempUpdate)*detJ;

end

% Check for convective boundary condition on right side and heat flux on left side
if testCase == 4 || testCase == 6

    mergeConstrainNodeIds = [leftBoundaryNodeId,rightBoundaryNodeId];

    % Check if element node IDs are present in boundary node IDs
    [isPresent, ~] = ismember(nodelIds, mergeConstrainNodeIds);

    if any(isPresent)
        gaussPoints2 = [sqrt(1/3),-sqrt(1/3)]; % Gauss points for edge integration

        for j = 1:2

            N = 1/2*[1-gaussPoints2(j), 1+gaussPoints2(j)]; % shape function for edge
            derivN = [-1/2,1/2];
            
            % Check presence of nodes on left boundary and right boundary
            isPresentLeft = ismember(nodelIds, leftBoundaryNodeId);
            isPresentRight = ismember(nodelIds, rightBoundaryNodeId);

            if sum(isPresentLeft) == 2
                % Replace non-matching elements with 0 for left boundary
                shapeEdge =  nodelIds;
                shapeEdge(~isPresentLeft) =  0;
                shapeEdge(isPresentLeft) = N;
                
                JEdge = derivN*nodeCoordinates(isPresentLeft,:);
                
                % Determinant of Jacobian for edge
                detJEdge = norm(JEdge);
                
                % Accumulate heat flux
                qFLux = qFLux+ shapeEdge'*heatFlux*detJEdge;

            elseif sum(isPresentRight) == 2


                % Replace non-matching elements with 0 for right boundary
                shapeEdge =  nodelIds;
                shapeEdge(~isPresentRight) =  0;
                shapeEdge(isPresentRight) = N;
                
                JEdge = derivN*nodeCoordinates(isPresentRight,:);
                
                % Determinant of Jacobian for edge
                detJEdge = norm(JEdge);
                
                % Accumulate convective matrix
                K_H = K_H + (shapeEdge'*shapeEdge)*heatCoeff*detJEdge;

                qAmbient = qAmbient+ shapeEdge'*heatCoeff*ambientTemp*detJEdge;


            end
        end
    end
end

% Check for radiation boundary condition on the top boundary and right boundary
if testCase == 7

    mergeConstrainNodeIds = [topBoundaryNodeId,rightBoundaryNodeId];

    % Check if element node IDs are present in boundary node IDs
    [isPresent, ~] = ismember(nodelIds, mergeConstrainNodeIds);

    if sum(isPresent) == 2
        gaussPoints2 = [sqrt(1/3),-sqrt(1/3)];% Gauss points for edge integration

        for j = 1:2

            N = 1/2*[1-gaussPoints2(j), 1+gaussPoints2(j)]; % shape function for edge
            derivN = [-1/2,1/2];

            % Replace non-matching elements with 0
            shapeEdge =  nodelIds;
            shapeEdge(~isPresent) =  0;
            shapeEdge(isPresent) = N;

            JEdge = derivN*nodeCoordinates(isPresent,:);
            
            % Determinant of Jacobian for edge
            detJEdge = norm(JEdge);

            temp = shapeEdge*tempUpdate;
            
            % Accumulate radiation matrix and residual vector
            K_R = K_R + 4*emissivity*stefanBoltzmann*(shapeEdge'*shapeEdge)*temp^3*detJEdge;
            qAmbientR = qAmbientR + shapeEdge'*emissivity*stefanBoltzmann*ambientTemp^4*detJEdge;
            Q_int_r = Q_int_r + emissivity*stefanBoltzmann*(shapeEdge'*shapeEdge)*tempUpdate.^4*detJEdge;

        end
    end

end

% Final stiffness matrix calculation
K_e = K_e + K_H + K_R;
Q_int_e= Q_int_e + K_H*tempUpdate + Q_int_r;
Q_ext_e= Q_ext_e + qAmbient + qFLux + qAmbientR;

% Final residual vector calculation
R_e = Q_int_e-Q_ext_e;


%Numerical Testing code.
if strcmp(numericalTest,'yes')
    fprintf("Area of unit sqaure using determinant of Jacobian: %d \n",area);
    
    disp("Fem B matrix temperature Gradient:");
    disp(tempGradient);
    disp("Fem Jacobian matrix:");
    disp(jac);
end

end
