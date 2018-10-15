function addConstLineLoad(startNodeId, endNodeId, nodeArray, modulus, direction)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

direction = direction ./ norm(direction);
load = direction .* modulus;

startCoord = nodeArray(startNodeId).getCoords();
endCoord = nodeArray(endNodeId).getCoords();

ii = 1;
    if startCoord(1) == endCoord(1)
        lengthY = abs(startCoord(2)-endCoord(2));
        for i = 1:length(nodeArray)
            if nodeArray(i).getX == startCoord(1)
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
        
    elseif startCoord(2) == endCoord(2)
        lengthX = abs(startCoord(1)-endCoord(1));
        for i = 1:length(nodeArray)
            if nodeArray(i).getY == startCoord(2)
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
        error('Line Load must be parallel to global X or Y Axsis.')
    end
end
