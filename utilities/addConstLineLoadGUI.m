function addConstLineLoad(startNodeId, endNodeId, nodeArray, fx, fy)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

load = [fx fy];

startCoord = nodeArray(startNodeId).getCoords();
endCoord = nodeArray(endNodeId).getCoords();
tol = 1e-10;

ii = 1;
    if  abs(startCoord(1)-endCoord(1)) < tol
        lengthY = abs(startCoord(2)-endCoord(2));
        for i = 1:length(nodeArray)
            if abs(nodeArray(i).getX-startCoord(1)) < tol
                loadNodes(ii) = nodeArray(i);
                ii = ii + 1;
            end
        end
        number_nodes = length(loadNodes);
        nodalLoad = lengthY * load / (number_nodes - 1);
        loadNodes.setDofLoad('DISPLACEMENT_X', nodalLoad(1));
        loadNodes.setDofLoad('DISPLACEMENT_Y', nodalLoad(2));
        nodeArray(startNodeId).setDofLoad('DISPLACEMENT_X', (nodalLoad(1) / 2));
        nodeArray(startNodeId).setDofLoad('DISPLACEMENT_Y', (nodalLoad(2) / 2));
        nodeArray(endNodeId).setDofLoad('DISPLACEMENT_X', (nodalLoad(1) / 2));
        nodeArray(endNodeId).setDofLoad('DISPLACEMENT_Y', (nodalLoad(2) / 2));
        
    elseif abs(startCoord(2)-endCoord(2)) < tol
        lengthX = abs(startCoord(1)-endCoord(1));
        for i = 1:length(nodeArray)
            if abs(nodeArray(i).getY-startCoord(2)) < tol
                loadNodes(ii) = nodeArray(i);
                ii = ii + 1;
            end
        end
        number_nodes = length(loadNodes);
        nodalLoad = lengthX * load / (number_nodes - 1);
        loadNodes.setDofLoad('DISPLACEMENT_X', nodalLoad(1));
        loadNodes.setDofLoad('DISPLACEMENT_Y', nodalLoad(2));
        nodeArray(startNodeId).setDofLoad('DISPLACEMENT_X', (nodalLoad(1) / 2));
        nodeArray(startNodeId).setDofLoad('DISPLACEMENT_Y', (nodalLoad(2) / 2));
        nodeArray(endNodeId).setDofLoad('DISPLACEMENT_X', (nodalLoad(1) / 2));
        nodeArray(endNodeId).setDofLoad('DISPLACEMENT_Y', (nodalLoad(2) / 2));

    else
        warndlg('Line Load must be parallel to global X or Y Axis.','Warning');
    end
end
