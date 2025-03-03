%--------------------------------------------------------------------------
% Title: Nonlinear Finite element implementation of thermoelastic Material Model.
% Name:  Vishal Soni 
% Matrikulation Numbber: 66158
% Purpose: Personal Programming Project
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Inputs :1) readFile
%         2) Loadcase
%         3) testCase
%         4) eigenValueTest
%         5) numericalTesting
%         6) temperatureDependent
%         7) steadyState
%--------------------------------------------------------------------------
% Outputs :
%           1) 2D contour plots of all displacements and temperature.

clear,clc;

% Open the file for reading
readFile = fopen('square100Element.inp', 'r');

%Calling extract_data_inp function to extract all data for meshes from inp file
%generated from GMSH
[nodeCoordinates,elementsConnectivity,leftBoundaryNodeId,rightBoundaryNodeId,topBoundaryNodeId,bottomBoundaryNodeId] = extractDataInpFile(readFile);

% Close the file
fclose(readFile);

%Tolerance for newton raphson loop
tolerance = 1e-6;

%Initialize variable for calcultaions
Totaldofs = size(nodeCoordinates,1)*3; %3 Degree of freedom for Temperature and displacements.
TGlobal = zeros(size(nodeCoordinates,1),1); %global Temperature values
UGlobal = zeros(size(nodeCoordinates,1)*2,1); %global displacement values
WGlobal = zeros(Totaldofs,1); %W global has displacement and temperature values together total 3 degree of freedom.
nElement = size(elementsConnectivity,1); %Number of elements


totalDisplacementdof= size(nodeCoordinates,1)*2; %Total Degree of freedom for Displacement.
totalTemperaturedof = size(nodeCoordinates,1); %Total Degree of freedom for Temperature.

%Dummy values for Convection Boundary condition
heatCoeff = 0;
heatFlux = 0;

%Dummy values for Radiation Boundary condition
ambientTemp = 0;
emissivity = 0;
stefanBoltzmann = 0;

% Uncomment the required type of test case.
testCase = 1; %Must be equal to zero if running different load cases
% testCase = 1; % Exact solution test for temperature distribution
% testCase = 2; % Exact solution test for temperature distribution
% testCase = 3; %patch test for heat transfer element
% testCase = 4; % Exact solution test for convective boundary conditon

% Uncomment the required type of load case.
%LoadCase = 0; %Must be equal to zero to run the testCases
%LoadCase = 5; %Temperature Boundary condition
%LoadCase = 6; %Heat flux and convective boundary condition
LoadCase = 0;  %Radiation Boundary condition

%Write 'yes' or 'no' to perform Zero eigen value and numerical testing.
eigenValueTest = 'yes';
numericalTesting = 'no';

%Write 'yes' or 'no' to consider temperature dependent material properties or not.
temperatureDependent = 'no';
%Write 'yes' or 'no' to perform steady state or transient analysis.
steadyState = 'yes';

% Check if the problem is steady-state
% Material properties are referred from ZhihuiLiu, ZhihuiLi , Qiang Ma. ”Nonlinear finite element algorithm for
% solving fully coupled thermomechanical problems under strong
% aerothermodynamic environment” 2023
if strcmp(steadyState, 'yes')
    
    % For steady-state, set the density (rho) and specific heat to zero
    rho = 0;
    specificHeat = 0;
    
    % In steady-state, the time step and start/end times are arbitrary
    deltaTime = 1;
    tStart = 1;
    tEnd = 1;

else
    
    % For transient analysis, assign physical values for density and specific heat
    rho = 9000;          % Material density (kg/m^3)
    specificHeat = 300;   % Specific heat capacity (J/kg·K)
    
    % Set the time step and time range for the transient analysis
    deltaTime = 100;      % Time step (s)
    tStart = 1;           % Start time (s)
    tEnd = 10000;         % End time (s)

end

numTimeSteps = (tEnd-(tStart-1))/deltaTime; % number  of time steps

nodalTempTimeStep = zeros(totalTemperaturedof,numTimeSteps+1); %Store temperature value for every time step.

%Exact Solution for test case =1 is referred from
% FINITE ELEMENT METHOD FOR HEAT TRANSFER
% PHENOMENON ON A CLOSED RECTANGULAR PLATE
% Collins O. Akeremale,Olusegun A Olaiju, Yeak Su Hoe 2020
if testCase == 1

    %Material property for thermal field.
    thermalConductivity = [1,0;0,1];

    %Boundary condition
    TGlobal(leftBoundaryNodeId) = nodeCoordinates(leftBoundaryNodeId,2)./(1 + nodeCoordinates(leftBoundaryNodeId,2).^2);
    TGlobal(rightBoundaryNodeId) = nodeCoordinates(rightBoundaryNodeId,2)./(4 + nodeCoordinates(rightBoundaryNodeId,2).^2);
    TGlobal(bottomBoundaryNodeId) = 0;
    TGlobal(topBoundaryNodeId) = 1./((1 + nodeCoordinates(topBoundaryNodeId,1)).^2 +1);

    %Exact solution for Temperature distribution
    TExact =((nodeCoordinates(:,2))./(((1+nodeCoordinates(:,1)).^2) +(nodeCoordinates(:,2).^2)));

    nodalTempTimeStep(:,1) = TGlobal;
    TGlobalinitial = nodalTempTimeStep(:,1);

end

%Exact Solution for test case = 2 is referred from
% FINITE ELEMENT METHOD FOR HEAT TRANSFER
% PHENOMENON ON A CLOSED RECTANGULAR PLATE
% Collins O. Akeremale,Olusegun A Olaiju, Yeak Su Hoe 2020
if testCase == 2

    %Material property for thermal field.
    thermalConductivity = [1,0;0,1];

    %Boundary condition
    TGlobal(leftBoundaryNodeId) = 0;
    TGlobal(rightBoundaryNodeId) =0;
    TGlobal(bottomBoundaryNodeId) = 0;
    TGlobal(topBoundaryNodeId) =(400*sin(pi.*nodeCoordinates(topBoundaryNodeId,1)))/pi;

    %Exact solution for Temperature distribution
    TExact =(((-400*exp(pi.*nodeCoordinates(:,2)))/(pi*(exp(-pi) - exp(pi))) + (400*exp(-pi.*nodeCoordinates(:,2)))/(pi*(exp(-pi) - exp(pi))))).*sin(pi.*nodeCoordinates(:,1));

    nodalTempTimeStep(:,1) = TGlobal;
    TGlobalinitial = nodalTempTimeStep(:,1);
end


%Patch test for thermal field is referred from 
%Abaqus https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/eif/ec24dfp4.inp
if testCase == 3

    %From Abaqus
    % https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/eif/ec24dfp4.inp
    nodeCoordinates = [0,0;
        0.24,0;
        0.24,0.12;
        0,0.12;
        0.04,0.02;
        0.18,0.03;
        0.08,0.08;
        0.16,0.08;
        ];

    % Element connectivity:
    elementsConnectivity = [
        1,2,6,5;  % Element 1
        1,5,7,4;  % Element 2
        5,6,8,7;  % Element 3
        6,2,3,8;  % Element 4
        7,8,3,4   % Element 5
        ];

    %Material property
    thermalConductivity = [4.85e-4,0;0,4.85e-4];

    TGlobal = zeros(size(nodeCoordinates,1),1);
    nElement = size(elementsConnectivity,1); %Number of elements
    totalTemperaturedof = size(nodeCoordinates,1); %Total Degree of freedom for Temperature.
    
    %Know temperatue values
    TGlobal(1) = 0;
    TGlobal(2) = 48;
    TGlobal(3) = 60;
    TGlobal(4) = 12;
    
    %Exact solution constant
    a = 0;
    b = 200;
    c = 100;

    %Exact solution
    TExact = a + b * nodeCoordinates(:, 1) + c * nodeCoordinates(:, 2);

    TGlobalinitial = TGlobal;
    
    %Plot Mesh
    plotMesh(nodeCoordinates,elementsConnectivity);

end

%Exact Solution for test case = 4 is referred \
% from solving direct and inverse heat conduction book
%exercise: 6.8
%Material properties are referred from Finite Element Method in
% Steady-State and Transient Heat Conduction
% Jan Taler and Paweł Ocłon
if testCase == 4 %Convective boundary condition

    %Covective Boundary condition.
    heatCoeff = 1000;
    heatFlux = 200000;
    ambientTemp = 20;

    thermalConductivity = [50,0;0,50];

    nodalTempTimeStep(:,1) = 20;
    TGlobalinitial = nodalTempTimeStep(:,1);

    lengthDomain = 1; %Already know


    % Extract the x-coordinates at y = 0.5
    y_fixed = 0.5;
    indices_at_y05 = find(abs(nodeCoordinates(:,2) - y_fixed) < tolerance);
    x_at_y05 = sort(nodeCoordinates(indices_at_y05, 1));

    %Refrence from solving direct and inverse heat conduction book
    %exercise: 6.8
    TExact = (heatFlux*(lengthDomain - x_at_y05))/thermalConductivity(1,1) + (heatFlux/heatCoeff) + ambientTemp;

end


%Abaqus verification of code for couple conduction problem.
%Validation manual from ”SIMSolid” for material properties
if LoadCase == 5

    %Material property for thermal field.
    thermalConductivity = [43,0;0,43];

    %Material property for mechanical field
    E = 2e11;% Youngs modulus
    nu = 0.3; %Poisson ratio

    %Thermal expansion.
    alpha = 1.2e-05;

    %Temperature Boundary condition
    WGlobal(topBoundaryNodeId+totalDisplacementdof) = 100;
    WGlobal(bottomBoundaryNodeId+totalDisplacementdof) = 0;

    %Displacement Boundary condition
    disprightnodelIds = [2*rightBoundaryNodeId-1; 2*rightBoundaryNodeId];
    displeftnodelIds = [2*leftBoundaryNodeId-1; 2*leftBoundaryNodeId];
    WGlobal(disprightnodelIds(:)) = 0;
    WGlobal(displeftnodelIds(:)) = 0;

    nodalWglobalTimeStep = zeros(Totaldofs,numTimeSteps+1);

    nodalWglobalTimeStep(:,1) = WGlobal;
    WGlobalinitial = nodalWglobalTimeStep(:,1);

end

%Abaqus coupled problem verification of code
% for heat flux and convective boundary condition.
%Material properties for convection boundary condition are referred from Finite Element Method in
% Steady-State and Transient Heat Conduction
% Jan Taler and Paweł Ocłon
%Validation manual from ”SIMSolid” for Mechanical material properties
if LoadCase == 6

    %Convective Boundary condition values.
    heatCoeff = 1000;
    heatFlux = 200000;
    ambientTemp = 20;

    %Material property for thermal field.
    thermalConductivity = [50,0;0,50];


    %Material property for mechanical field
    E = 2e11;% Youngs modulus
    nu = 0.3; %Poisson ratio

    %Thermal expansion.
    alpha = 1.2e-05;

    %Displacement Boundary condition
    disprightnodelIds = [2*rightBoundaryNodeId-1; 2*rightBoundaryNodeId];

    WGlobal(disprightnodelIds(:)) = 0;

    nodalWglobalTimeStep = zeros(Totaldofs,numTimeSteps+1);

    nodalWglobalTimeStep(:,1) = WGlobal;
    WGlobalinitial = nodalWglobalTimeStep(:,1);


end

%Abaqus coupled problem verification of code from research paper
%For radiation boundary condition
% Material properties are referred from ZhihuiLiu, ZhihuiLi , Qiang Ma. ”Nonlinear finite element algorithm for
% solving fully coupled thermomechanical problems under strong
% aerothermodynamic environment” 2023
if LoadCase == 7

    %Material property for thermal field.
    thermalConductivity = [400,0;0,400];


    %Material property for mechanical field
    E = 2e7;% Youngs modulus
    nu = 0.3; %Poisson ratio

    %Thermal expansion.
    alpha = 0.001;

    emissivity = 0.5;
    stefanBoltzmann = 5.67e-08;
    ambientTemp = 20;

    %Temperature boundary condition.
    WGlobal(leftBoundaryNodeId+totalDisplacementdof) = 200;

    %Displacement Boundary condition
    disprightnodelIds = [2*rightBoundaryNodeId-1; 2*rightBoundaryNodeId];
    disptopnodelIds = [2*topBoundaryNodeId-1; 2*topBoundaryNodeId];
    dispbottomnodelIds = [2*bottomBoundaryNodeId-1; 2*bottomBoundaryNodeId];
    WGlobal(dispbottomnodelIds(2,:)) = 0;
    WGlobal(disprightnodelIds(1,:)) = 0;
    WGlobal(disptopnodelIds(2,:)) = 0;


    % Extract the x-coordinates at y = 0.5
    y_fixed = 0.5;
    indices_at_y05 = find(abs(nodeCoordinates(:,2) - y_fixed) < tolerance);
    x_at_y05 = sort(nodeCoordinates(indices_at_y05, 1));


    nodalWglobalTimeStep(:,1) = WGlobal;
    WGlobalinitial = nodalWglobalTimeStep(:,1);


end

% Check if the problem is numerical Testing
if strcmp (numericalTesting,'yes')
    nodeCoordinates = [0 0; 1 0; 1 1; 0 1];
    elementsConnectivity = [1,2,3,4];
    nElement = 1;
    plotMesh(nodeCoordinates,elementsConnectivity);
end

%Time integration loop.
for t = 1:numTimeSteps

    fprintf('Time Increment %d: Step Time = %ds \n', t, t*deltaTime);

    if testCase == 1 || testCase == 2 || testCase == 3 || testCase == 4
        TGlobal = TGlobalinitial;
    else
        WGlobal = WGlobalinitial;
    end

    %newton raphson loop
    k =1;
    while 1

        KGlobalT = zeros(totalTemperaturedof,totalTemperaturedof);
        KGlobalM = zeros(totalDisplacementdof);
        KGlobalMT = zeros(totalDisplacementdof,totalTemperaturedof);

        RGlobalT = zeros(totalTemperaturedof,1);
        RGlobalM = zeros(totalDisplacementdof,1);
        RGlobalMT = RGlobalM;

        %Element loop
        for e = 1:nElement

            nodelIds = (elementsConnectivity(e,:));

            nodeEleCor = nodeCoordinates(nodelIds,:);

            if testCase == 1 || testCase == 2 || testCase == 3 || testCase == 4

                Te = TGlobalinitial(nodelIds,:);
                tempUpdate = TGlobal(nodelIds,:);
                
                %Thermal Element routine
                [KeTemp,ReTemp] = elementThermalfem(nodelIds,leftBoundaryNodeId,rightBoundaryNodeId,topBoundaryNodeId,...
                    nodeEleCor,thermalConductivity,Te,testCase,heatFlux,heatCoeff,ambientTemp,...
                    rho,specificHeat,deltaTime,tempUpdate,temperatureDependent,emissivity,stefanBoltzmann,numericalTesting);

                %Assembly for stiffness matrix and Residual force
                KGlobalT(nodelIds,nodelIds) = KGlobalT(nodelIds,nodelIds) + KeTemp;
                RGlobalT(nodelIds,:) = RGlobalT(nodelIds,:) + ReTemp;

            elseif LoadCase == 5 || LoadCase == 6 || LoadCase == 7


                disp_nodelIds = [2*nodelIds-1; 2*nodelIds];


                Te = WGlobalinitial(nodelIds+totalDisplacementdof,:);
                tempUpdate = WGlobal(nodelIds+totalDisplacementdof,:);

                Ue = WGlobal(disp_nodelIds(:),:);
                
                %Thermal Element routine
                [KeTemp,ReTemp] = elementThermalfem(nodelIds,leftBoundaryNodeId,rightBoundaryNodeId,topBoundaryNodeId,...
                    nodeEleCor,thermalConductivity,Te,LoadCase,heatFlux,heatCoeff,ambientTemp,...
                    rho,specificHeat,deltaTime,tempUpdate,temperatureDependent,emissivity,stefanBoltzmann,numericalTesting);
                
                %Mechanical Element routine
                [KeMech,ReMech] = elementMechfem(nodeEleCor,Ue,tempUpdate,E,nu,temperatureDependent);
                
                %Coupled Thermal Mechanical Element routine
                [KeMechTemp,ReMechTemp] = elementMechThermfem(nodeEleCor,tempUpdate,E,nu,alpha,temperatureDependent);

                %Assembly for stiffness matrix and Residual force
                KGlobalT(nodelIds,nodelIds) = KGlobalT(nodelIds,nodelIds) + KeTemp;
                KGlobalM(disp_nodelIds(:),disp_nodelIds(:)) = KGlobalM(disp_nodelIds(:),disp_nodelIds(:)) + KeMech;
                KGlobalMT(disp_nodelIds(:),nodelIds) = KGlobalMT(disp_nodelIds(:),nodelIds) + KeMechTemp;

                RGlobalT(nodelIds,:) = RGlobalT(nodelIds,:) + ReTemp;
                RGlobalM(disp_nodelIds(:),:) = RGlobalM(disp_nodelIds(:),:) + ReMech;
                RGlobalMT(disp_nodelIds(:),:) = RGlobalMT(disp_nodelIds(:),:) + ReMechTemp;

            end


        end

        if testCase == 1 || testCase == 2
            %Handling constrained Boundary conditions
            mergeConstrainNodeIds = [leftBoundaryNodeId,rightBoundaryNodeId,topBoundaryNodeId,bottomBoundaryNodeId];
            constrainNodeIds = unique(mergeConstrainNodeIds);

            %Handling unconstrained Boundary conditions
            unconstrainNodeIds = setdiff(1:size(nodeCoordinates,1), constrainNodeIds);

            % %Sparse matrix
            KGlobalTSparse = sparse(KGlobalT);
            RGlobalTSparse = sparse(RGlobalT);

            %Reducing rows and columns for Stiffness matrix
            KReduce = KGlobalTSparse(unconstrainNodeIds,unconstrainNodeIds);


            %Reducing rows and columns for Residual matrix
            ResReduce = RGlobalTSparse(unconstrainNodeIds);
            
            %Solution
            deltaW = KReduce\-ResReduce;
        end

        if testCase == 3

            constrainNodeIds = [1,2,3,4];
            %Handling unconstrained Boundary conditions
            unconstrainNodeIds = setdiff(1:size(nodeCoordinates,1), constrainNodeIds);

            %Reducing rows and columns for Stiffness matrix
            KReduce = KGlobalT(unconstrainNodeIds,unconstrainNodeIds);


            %Reducing rows and columns for Residual matrix
            ResReduce = RGlobalT(unconstrainNodeIds);
            
            %Solution
            deltaW = KReduce\-ResReduce;

        end


        if testCase == 4
            
            % %Sparse matrix
            KReduce = sparse(KGlobalT);
            ResReduce = sparse(RGlobalT);

            %solution
            deltaW = KReduce \-ResReduce;

        end


        if LoadCase == 5

            %Constrained Nodal Ids
            mergeConstrainNodeIds = [topBoundaryNodeId+totalDisplacementdof, bottomBoundaryNodeId+totalDisplacementdof, disprightnodelIds(:)', displeftnodelIds(:)'];
            constrainNodeIds = unique(mergeConstrainNodeIds);

            %Unconstrained Nodal Ids
            unconstrainNodeIds = setdiff(1:size(nodeCoordinates,1)*3, constrainNodeIds);

            %Generating sparse matrix
            KGlobalMSparse = sparse(KGlobalM);
            KGlobalMTSparse = sparse(KGlobalMT);
            KGlobalTSparse = sparse(KGlobalT);
            RGlobalUSparse = sparse(RGlobalM-RGlobalMT);
            RGlobalTSparse = sparse(RGlobalT);

            % Build the Fully coupled matrix and the residual vector
            coupleK = [KGlobalMSparse, -KGlobalMTSparse; KGlobalMTSparse', KGlobalTSparse];
            coupleRes = [RGlobalUSparse;RGlobalTSparse];

            % Reduce the matrix and residual vector for boundary conditions.
            KReduce = coupleK(unconstrainNodeIds, unconstrainNodeIds);  % Remove rows and columns
            ResReduce = coupleRes(unconstrainNodeIds);  % Remove corresponding entries in the residual vector

            %Solution
            deltaW = KReduce\-ResReduce;
        end

        if LoadCase == 6

            %Constrained Nodal Ids
            mergeConstrainNodeIds = disprightnodelIds(:)';
            constrainNodeIds = unique(mergeConstrainNodeIds);

            %Unconstrained Nodal Ids
            unconstrainNodeIds = setdiff(1:size(nodeCoordinates,1)*3, constrainNodeIds);

            %Generating sparse matrix
            KGlobalMSparse = sparse(KGlobalM);
            KGlobalMTSparse = sparse(KGlobalMT);
            KGlobalTSparse = sparse(KGlobalT);
            RGlobalUSparse = sparse(RGlobalM-RGlobalMT);
            RGlobalTSparse = sparse(RGlobalT);

            % Build the  Fully coupled matrix and the residual vector
            coupleK = [KGlobalMSparse, -KGlobalMTSparse; KGlobalMTSparse', KGlobalTSparse];
            coupleRes = [RGlobalUSparse;RGlobalTSparse];


            % Reduce the matrix and residual vector for boundary conditions.
            KReduce = coupleK(unconstrainNodeIds, unconstrainNodeIds);  % Remove rows and columns
            ResReduce = coupleRes(unconstrainNodeIds);  % Remove corresponding entries in the residual vector

            %Solution
            deltaW = sparse(KReduce)\-ResReduce;

        end

        if LoadCase == 7

            mergeConstrainNodeIds = [leftBoundaryNodeId+totalDisplacementdof, dispbottomnodelIds(2,:), disprightnodelIds(1,:),disptopnodelIds(2,:)];
            constrainNodeIds = unique(mergeConstrainNodeIds);

            %Unconstrained Nodal Ids
            unconstrainNodeIds = setdiff(1:size(nodeCoordinates,1)*3, constrainNodeIds);


            %Generating sparse matrix
            KGlobalMSparse = sparse(KGlobalM);
            KGlobalMTSparse = sparse(KGlobalMT);
            KGlobalTSparse = sparse(KGlobalT);
            RGlobalUSparse = sparse(RGlobalM-RGlobalMT);
            RGlobalTSparse = sparse(RGlobalT);

            % Build the Fully coupled matrix and the residual vector
            coupleK = [KGlobalMSparse, -KGlobalMTSparse; KGlobalMTSparse', KGlobalTSparse];
            coupleRes = [RGlobalUSparse;RGlobalTSparse];


            % Reduce the matrix and residual vector for boundary conditions.
            KReduce = coupleK(unconstrainNodeIds, unconstrainNodeIds);  % Remove rows and columns
            ResReduce = coupleRes(unconstrainNodeIds);  % Remove corresponding entries in the residual vector

            %Solution
            deltaW = sparse(KReduce)\-ResReduce;

        end


        if testCase == 1 || testCase == 2 || testCase == 3
            %Print newton raphson iterations and norms values.
            fprintf('Iteration %d: Temperature Norm = %.6e  Residual Norm = %.6e\n', k, norm(deltaW),norm(ResReduce));
            TGlobal(unconstrainNodeIds,:) = TGlobal(unconstrainNodeIds,:) + deltaW;
        elseif testCase == 4
            %Print newton raphson iterations and norms values.
            fprintf('Iteration %d: Temperature Norm = %.6e  Residual Norm = %.6e\n', k, norm(deltaW),norm(ResReduce));
            TGlobal = TGlobal + deltaW;
        elseif LoadCase == 5 || LoadCase == 6 || LoadCase == 7
            %Print newton raphson iterations and norms values.
            fprintf('Iteration %d: Couple Norm = %.6e  Residual Norm = %.6e\n', k, norm(deltaW),norm(ResReduce));
            WGlobal(unconstrainNodeIds,:) = WGlobal(unconstrainNodeIds,:) + deltaW;
        end

        %Checking criteria for newton raphson loop.
        if norm(deltaW) < tolerance ||  norm(ResReduce) < tolerance

            fprintf('Convergence achieved with tolerance = %e\n',tolerance);

            break;

        else
            k = k+1;
        end
    end

    if testCase == 1 || testCase == 2 || testCase == 4

        nodalTempTimeStep(:,t+1) = TGlobal;

        TGlobalinitial = TGlobal;

    else

        nodalWglobalTimeStep(:,t+1) = WGlobal;
        WGlobalinitial = WGlobal;

    end

end

%---------------------Stress, Strain Recovery-----------------------------%
% Check if the LoadCase corresponds to 5, 6, or 7
if LoadCase == 5 || LoadCase == 6 || LoadCase == 7

    % Stress and Strain recovery for all elements
    % Initialize arrays to store total strain, thermal strain, elastic strain, and thermal stress
    totalStrain = zeros(3,4,nElement);
    thermalStrain = zeros(3,4,nElement);
    elasticStrain = zeros(3,4,nElement);
    thermalStress = zeros(3,4,nElement);
    
     % Loop through each element in the mesh
    for e = 1:nElement

        % Get the node IDs for the current element from the connectivity matrix
        nodelIds = (elementsConnectivity(e,:));

        % Get the coordinates of the nodes for the current element
        nodeEleCor = nodeCoordinates(nodelIds,:);

        % Retrieve the displacement degrees of freedom (DOFs) for the nodes in the current element
        disp_nodelIds = [2*nodelIds-1; 2*nodelIds];
            
        % Retrieve temperature values (Te) for the current element's nodes
        Te = WGlobal(nodelIds+totalDisplacementdof,:);

        % Retrieve displacement values (Ue) for the current element's nodes
        Ue = WGlobal(disp_nodelIds(:),:);
        
        % Call the recovery function to compute strains and stresses for the current element
        [totalStrainEle,thermalStrainEle,elasticStrainEle,thermalStressEle] = recoveryStressStrain(nodeEleCor,Ue,Te,E,nu,alpha,temperatureDependent);
        
        % Store the computed strains and stresses for the current element
        totalStrain(:,:,e) = totalStrainEle;
        thermalStrain(:,:,e) = thermalStrainEle;
        elasticStrain(:,:,e) = elasticStrainEle;
        thermalStress(:,:,e) = thermalStressEle;

    end
end

%--------------------------------------------------------------------------
%--------------------------Plotting----------------------------------------
%Calling plotMesh function to plot discretize domain.
%plotMesh(nodeCoordinates,elementsConnectivity);
%
%


x=nodeCoordinates(:,1); %X coordinates for nodal Ids
y= nodeCoordinates(:,2);  %Y coordinates for nodal Ids

% Check if the test case is 1, 2, or 3 and numerical testing is set to 'no'
if (testCase == 1 || testCase == 2 || testCase == 3) && strcmp(numericalTesting,'no')
    
    % Plot a line plot for both Exact and FEM Global solution of temperature values
    figure(2); % Create a new figure for the line plot
    
    % Plot the exact solution of temperature at nodes
    plot(1:size(nodeCoordinates,1), TExact, 'o', 'DisplayName', 'Exact Solution');
    hold on; % Retain current plot when adding new plots
    
    % Plot the FEM solution of temperature at nodes
    plot(1:size(nodeCoordinates,1), TGlobal, '-x', 'DisplayName', 'FEM Solution');
    
    % Add labels and title to the plot
    title('Temperature (°C) vs. Nodal Numbers'); % Title of the plot
    xlabel('Nodal Numbers'); % Label for the x-axis
    ylabel('Temperature (°C)'); % Label for the y-axis

    % Add a legend to identify the Exact and FEM solution lines
    legend;

    % Display the error norm (difference) between the exact solution and the FEM solution
    disp("Error norm between Exact solution and FEM solution for Temperature distribution:");
    disp(norm(TExact - TGlobal)); % Compute and display the error norm

    % If testCase is 3, display the FEM temperature distribution
    if testCase == 3
        disp("Temperature distribution");
        disp(TGlobal); % Display the FEM solution of temperature at all nodes
    end

    % Create a new figure for visualizing the temperature distribution in 2D
    figure(3);
    
    % Call the function to visualize the temperature distribution using trisurf
    temperatureVisualization(x, y, elementsConnectivity, TGlobal);
end


% Check if the test case is 4
if testCase == 4

    % Sort the global temperature values at y = 0.5 in descending order
    sortTempdescend = sort(TGlobal(indices_at_y05), 'descend');

    % Create a new figure for plotting the temperature distribution along y = 0.5
    figure(2);
    
    % Plot the exact temperature solution along the x coordinate at y = 0.5
    plot(x_at_y05, TExact, 'o', 'DisplayName', 'Exact Solution');
    hold on; % Retain current plot when adding new plots
    
    % Plot the FEM solution along the x coordinate at y = 0.5 (sorted in descending order)
    plot(x_at_y05, sortTempdescend, '-x', 'DisplayName', 'FEM Solution');
    
    % Add labels and title to the plot
    title('Temperature (°C) at constant y = 0.5[m]'); % Title of the plot
    xlabel('X Coordinate[m]'); % Label for the x-axis
    ylabel('Temperature (°C)'); % Label for the y-axis

    % Add a legend to identify the Exact and FEM solution lines
    legend;

    % Display the error norm (difference) between the exact solution and the sorted FEM solution
    disp("Error norm between Exact solution and FEM solution for Temperature distribution:");
    disp(norm(TExact - sortTempdescend)); % Compute and display the error norm

    % Create a new figure for visualizing the 2D temperature distribution
    figure(3);
    
    % Call the function to visualize the temperature distribution using trisurf
    temperatureVisualization(x, y, elementsConnectivity, TGlobal);

end


if LoadCase == 5 || LoadCase == 6 || LoadCase == 7

    scale_factor = 0; % Adjust the scale if needed

    % Compute deformed positions
    X_deformed = x + scale_factor*WGlobal(1:2:totalDisplacementdof);
    Y_deformed = y + scale_factor*WGlobal(2:2:totalDisplacementdof);
    
    %Plot displacement distribution in X direction U1
    figure(3);
    trisurf(elementsConnectivity, X_deformed,Y_deformed,  WGlobal(1:2:totalDisplacementdof), 'EdgeColor', 'k');
    % Create a custom colormap with a gradient from blue to cyan to yellow
    colormap turbo;
    % Add a colorbar
    colorbar;

    % Add labels and title
    title('2D U1 Distribution (m)');
    xlabel('X Coordinate(m)');
    ylabel('Y Coordinate(m)');
    view(2); % 3D view
    shading interp;

    %Plot displacement distribution in y direction U2
    figure(4);
    trisurf(elementsConnectivity, X_deformed,Y_deformed,WGlobal(2:2:totalDisplacementdof), 'EdgeColor', 'k');

    % Create a custom colormap with a gradient from blue to cyan to yellow
    colormap turbo;
    % Add a colorbar
    colorbar;

    % Add labels and title
    title('2D U2 Distribution (m)');
    xlabel('X Coordinate(m)');
    ylabel('Y Coordinate(m)');
    view(2); % 3D view
    shading interp;
    
    %Plot Temperature distribution
    figure(5);
    trisurf(elementsConnectivity, X_deformed,Y_deformed,  WGlobal(totalDisplacementdof+1:end), 'EdgeColor', 'k');
    % Create a custom colormap with a gradient from blue to cyan to yellow
    colormap turbo;
    % Add a colorbar
    colorbar;

    % Add labels and title
    title('2D Temperature Distribution (°C)');
    xlabel('X Coordinate(m)');
    ylabel('Y Coordinate(m)');
    view(2); % 3D view
    shading interp;
    
    %Plot displacement Magnitude distribution
    figure(6);
    % Compute total displacement magnitude
    U = sqrt(WGlobal(1:2:totalDisplacementdof).^2 + WGlobal(2:2:totalDisplacementdof).^2);
    
    trisurf(elementsConnectivity, X_deformed, Y_deformed, U, 'EdgeColor', 'none');
    title('2D Displacement Magnitude');
    xlabel('X Coordinate(m)');
    ylabel('Y Coordinate(m)');
    colormap turbo;
    colorbar;
    view(2); % Set the view to 2D for a contour plot
    shading interp;

    hold off;

end

% Check if the load case is 7
if LoadCase == 7

    % If it's a transient case (not steady-state)
    if strcmp(steadyState, 'no')

        % Sort the temperature values at y = 0.5 for specific time steps (10th, 25th, 50th time steps)
        sortTempdescend10 = sort(nodalWglobalTimeStep(indices_at_y05 + totalDisplacementdof, 10), 'descend');
        sortTempdescend25 = sort(nodalWglobalTimeStep(indices_at_y05 + totalDisplacementdof, 25), 'descend');
        sortTempdescend50 = sort(nodalWglobalTimeStep(indices_at_y05 + totalDisplacementdof, 50), 'descend');

        % Sort the temperature values at y = 0.5 for the final time step (10000s)
        sortTempdescend = sort(WGlobal(indices_at_y05 + totalDisplacementdof), 'descend');

        % Create a new figure to plot the temperature distribution at different time steps
        figure(2);
        
        % Plot the FEM solution at 1000 seconds
        plot(x_at_y05, sortTempdescend10, '-x', 'DisplayName', 'FEM Solution 1000s');
        hold on; % Retain current plot when adding new plots
        
        % Plot the FEM solution at 2500 seconds
        plot(x_at_y05, sortTempdescend25, '-x', 'DisplayName', 'FEM Solution 2500s');
        
        % Plot the FEM solution at 5000 seconds
        plot(x_at_y05, sortTempdescend50, '-x', 'DisplayName', 'FEM Solution 5000s');
        
        % Plot the FEM solution at 10000 seconds
        plot(x_at_y05, sortTempdescend, '-x', 'DisplayName', 'FEM Solution 10000s');

        % Add labels and title to the plot
        title('Temperature (°C) at constant y = 0.5[m]');
        xlabel('X Coordinate[m]'); % Label for the x-axis
        ylabel('Temperature (°C)'); % Label for the y-axis

        % Add a legend to identify the lines for different time steps
        legend('Location', 'southwest'); % Place the legend at the bottom-left corner

        hold off; % Release the plot

        % Define the displacement node IDs for the left boundary
        displeftnodelIds = [2*leftBoundaryNodeId-1; 2*leftBoundaryNodeId];

        % Extract the displacement values in the x and y directions for the left boundary nodes
        UxLeft = nodalWglobalTimeStep(displeftnodelIds(1,:), :); % Displacements in the x-direction
        UyLeft = nodalWglobalTimeStep(displeftnodelIds(2,:), :); % Displacements in the y-direction

        % Calculate the magnitude of displacement for each node at every time step
        Umagnitude = sqrt(UxLeft.^2 + UyLeft.^2);

        % Create a new figure for plotting the displacement variation over time
        figure(7);
        
        % Plot the displacement magnitude at the second left boundary node over time
        plot((1:size(Umagnitude, 2)) * deltaTime, Umagnitude(2,:), '-x', 'DisplayName', 'FEM Solution');
        
        % Add labels and title to the plot
        title('Displacement Variation on the left surface');
        xlabel('Time [s]'); % Label for the x-axis (time)
        ylabel('Displacement [m]'); % Label for the y-axis (displacement magnitude)

        % Add a legend to identify the FEM solution line
        legend('Location', 'southeast'); % Place the legend at the bottom-right corner

    end
end


% Check if the eigenvalue test is enabled
if strcmp(eigenValueTest, 'yes')

    % Compute the eigenvalues of the global thermal stiffness matrix
    eigenValue = eig(KGlobalT);
    
    % Display the computed eigenvalues
    disp("Eigen values for thermal element routine:");
    disp(eigenValue);

    % Calculate the number of eigenvalues that are approximately zero (less than 1e-5)
    ZeroEigSum = sum(eig(KGlobalT) < 1e-5);
    
    % Display the total number of zero eigenvalues
    disp("Total Zero eigen value for thermal element routine:");
    disp(ZeroEigSum);

end