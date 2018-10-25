function addLineBC(startNodeId, endNodeId, nodeArray)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

startCoord = nodeArray(startNodeId).getCoords();
endCoord = nodeArray(endNodeId).getCoords();
tol = 1e-10;

ii = 1;
    if abs(startCoord(1)-endCoord(1)) < tol
        for i = 1:length(nodeArray)
            if abs(nodeArray(i).getX-startCoord(1)) < tol
                boundary(ii) = nodeArray(i);
                ii = ii + 1;
            end
        end
        boundary.fixDof('DISPLACEMENT_X');
        boundary.fixDof('DISPLACEMENT_Y');
        
    elseif abs(startCoord(2)-endCoord(2)) < tol
        for i = 1:length(nodeArray)
            if abs(nodeArray(i).getY-startCoord(2)) < tol
                boundary(ii) = nodeArray(i);
                ii = ii + 1;
            end
        end
        boundary.fixDof('DISPLACEMENT_X');
        boundary.fixDof('DISPLACEMENT_Y');

    else
        error('LineBC must be parallel to global X or Y Axis.')
    end
end

