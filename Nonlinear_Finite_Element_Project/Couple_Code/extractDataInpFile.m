% extractDataInpFile reads data from an Gmsh .inp file and extracts
% node coordinates, element connectivity, and boundary node IDs for
% specified sets.
%
% This function processes the .inp file to gather the relevant data 
% for finite element analysis, including node coordinates and 
% connectivity of elements. It also identifies node sets corresponding 
% to different boundary conditions.
%
% Syntax:
%   [nodeCoordinates, elementsConnectivity, leftBoundaryNodeId, rightBoundaryNodeId, topBoundaryNodeId, bottomBoundaryNodeId] = extractDataInpFile(readFile)
%
% Inputs:
%   readFile           - A file identifier for the .inp file to be read.
%
% Outputs:
%   nodeCoordinates     - An Nx2 matrix containing the coordinates of the nodes.
%
%   elementsConnectivity - An MxN matrix specifying the connectivity of the elements.
%
%   leftBoundaryNodeId  - Node IDs on the left boundary.
%
%   rightBoundaryNodeId - Node IDs on the right boundary.
%
%   topBoundaryNodeId   - Node IDs on the top boundary.
%
%   bottomBoundaryNodeId - Node IDs on the bottom boundary.
%
% Limitations:
%   - When creating the domain in GMSH, the physical group names for the edges 
%     must be specified as follows: 'inelt' for the left boundary, 'outlet' 
%     for the right boundary, 'top' for the top boundary, and 'bottom' for 
%     the bottom boundary. This naming convention is necessary for correct 
%     extraction of boundary node IDs.

function [nodeCoordinates, elementsConnectivity, leftBoundaryNodeId, rightBoundaryNodeId, topBoundaryNodeId, bottomBoundaryNodeId] = extractDataInpFile(readFile)

% Initialize empty arrays to extract data from .inp file
nodeCoordinates = [];
elementsConnectivity = [];
leftBoundaryNodeId = [];
rightBoundaryNodeId = [];
topBoundaryNodeId = [];
bottomBoundaryNodeId = [];

% Initialize variables for storing node set names and IDs
setName = '';
ids = [];

% Read the .inp file line by line until the end of the file
while ~feof(readFile)
    line = strtrim(fgets(readFile)); % Read a line and trim whitespace

    % Check for the nodes section
    if startsWith(line, '*NODE', 'IgnoreCase', true)
        line = strtrim(fgets(readFile)); % Read the next line for node data
        while ~contains(line, "*")
            % Extract node data (Node ID, X, Y) and append to nodeCoordinates
            nodeData = sscanf(line, '%f, %f, %f, %f')';
            nodeCoordinates = [nodeCoordinates; nodeData(2:3)];
            line = strtrim(fgets(readFile)); % Read next line
        end
        continue;
    end

    % Check for the elements section of type CPS4 and ELSET=Surface1
    if contains(line, 'type=CPS4', 'IgnoreCase', true) && contains(line, 'ELSET=Surface1', 'IgnoreCase', true)
        line = strtrim(fgets(readFile)); % Read the next line for element data
        while ~contains(line, "*ELSET")
            % Extract element connectivity data and append to elementsConnectivity
            elementData = sscanf(line, '%f, %f, %f, %f, %f')';
            elementsConnectivity = [elementsConnectivity; elementData(2:end)];
            line = strtrim(fgets(readFile)); % Read next line
        end
        continue;
    end

    % Check if the line starts a new node set
    if startsWith(line, '*NSET,NSET=')
        % If there was a previous node set, store its IDs in the correct array
        if strcmp(setName, 'inlet')
            leftBoundaryNodeId = ids;
        elseif strcmp(setName, 'outlet')
            rightBoundaryNodeId = ids;
        elseif strcmp(setName, 'top')
            topBoundaryNodeId = ids;
        elseif strcmp(setName, 'bottom')
            bottomBoundaryNodeId = ids;
        end

        % Extract the new set name and reset the IDs for the next set
        setName = extractAfter(line, '*NSET,NSET=');
        ids = [];

        % Process lines that contain node IDs
    elseif ~isempty(line) && ~startsWith(line, '*')
        % Append node IDs from the current line to the ids array
        ids = [ids, sscanf(line, '%f,')'];
    end
end

% Store the last node set if it was the bottom boundary
if strcmp(setName, 'bottom')
    bottomBoundaryNodeId = ids;
end

end
