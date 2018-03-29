function addLineBC(startNodeId, endNodeId, nodeArray)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

startCoord = nodeArray(startNodeId).getCoords();
endCoord = nodeArray(endNodeId).getCoords();

ii = 1;
    if startCoord(1) == endCoord(1)
        lengthY = abs(startCoord(2)-endCoord(2));
        for i = 1:length(nodeArray)
            if nodeArray(i).getX == startCoord(1)
                boundary(ii) = nodeArray(i);
                ii = ii + 1;
            end
        end
        boundary.fixDof('DISPLACEMENT_X');
        boundary.fixDof('DISPLACEMENT_Y');
        
    elseif startCoord(2) == endCoord(2)
        lengthX = abs(startCoord(1)-endCoord(1));
        for i = 1:length(nodeArray)
            if nodeArray(i).getY == startCoord(2)
                boundary(ii) = nodeArray(i);
                ii = ii + 1;
            end
        end
        boundary.fixDof('DISPLACEMENT_X');
        boundary.fixDof('DISPLACEMENT_Y');

    else
        error('LineBC must be parallel to global X or Y Axsis.')
    end
end

