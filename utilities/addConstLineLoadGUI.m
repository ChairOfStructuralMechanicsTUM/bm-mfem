function addConstLineLoad(startNodeId, endNodeId, nodeArray, fx, fy)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

load = [fx fy];

startCoord = nodeArray(startNodeId).getCoords();
endCoord = nodeArray(endNodeId).getCoords();
tol = 1e-10;

ii = 1;
    if  abs(startCoord(1)-endCoord(1)) < tol %vertical
        lengthY = abs(startCoord(2)-endCoord(2));
        for i = 1:length(nodeArray)
            if abs(nodeArray(i).getX-startCoord(1)) < tol
                selectedNodes(ii) = nodeArray(i);
                ii = ii + 1;
            end
        end
        [y,idx] = sort(selectedNodes.getY);
        lengthA = 0;
        for i = 1:length(y)
            if i < length(y)
                lengthB = (y(i+1)-y(i))/2;
            else
                lengthB = 0;
            end
            factor = (lengthA + lengthB)/lengthY;
            selectedNodes(idx(i)).setDofLoad('DISPLACEMENT_X', factor * load(1));
            selectedNodes(idx(i)).setDofLoad('DISPLACEMENT_Y', factor * load(2));
            lengthA = lengthB;
        end
        
    elseif abs(startCoord(2)-endCoord(2)) < tol %horizontal
        lengthX = abs(startCoord(1)-endCoord(1));
        for i = 1:length(nodeArray)
            if abs(nodeArray(i).getY-startCoord(2)) < tol
                selectedNodes(ii) = nodeArray(i);
                ii = ii + 1;
            end
        end
        [x,idx] = sort(selectedNodes.getX);
        lengthA = 0;
        for i = 1:length(x)
            if i < length(x)
                lengthB = (x(i+1)-x(i))/2;
            else
                lengthB = 0;
            end
            factor = (lengthA + lengthB)/lengthX;
            selectedNodes(idx(i)).setDofLoad('DISPLACEMENT_X', factor * load(1));
            selectedNodes(idx(i)).setDofLoad('DISPLACEMENT_Y', factor * load(2));
            lengthA = lengthB;
        end
%         number_nodes = length(selectedNodes);
%         nodalLoad = lengthX * load / (number_nodes - 1);
%         selectedNodes.setDofLoad('DISPLACEMENT_X', nodalLoad(1));
%         selectedNodes.setDofLoad('DISPLACEMENT_Y', nodalLoad(2));
%         nodeArray(startNodeId).setDofLoad('DISPLACEMENT_X', (nodalLoad(1) / 2));
%         nodeArray(startNodeId).setDofLoad('DISPLACEMENT_Y', (nodalLoad(2) / 2));
%         nodeArray(endNodeId).setDofLoad('DISPLACEMENT_X', (nodalLoad(1) / 2));
%         nodeArray(endNodeId).setDofLoad('DISPLACEMENT_Y', (nodalLoad(2) / 2));

    else
        warndlg('Line Load must be parallel to global X or Y Axis.','Warning');
    end
end
