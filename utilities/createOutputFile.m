function [ outputFile ] = createOutputFile( Model ) 
%CREATEOUTPUTFILE Creates an output file 
%   including information about geometry, deformation and stresses


% Create TextFile
fileID = fopen('femCalculations.txt','w');


% Print Node Coordinates
fprintf(fileID,'Node Coordinates: \r\n \r\n');

nodeArray = getAllNodes(Model);

ii = 1:1:size(Model.getAllNodes,2);
A = [getId(nodeArray(ii)); getX(nodeArray(ii)); getY(nodeArray(ii)); getZ(nodeArray(ii))];

fprintf(fileID,'%6s %12s %12s %12s\r\n','node','x-coord','y-coord','z-coord');
fprintf(fileID,'%6.0f %12.2f %12.2f %12.2f\r\n',A);
fprintf(fileID,'\r\n \r\n');

% Print Element Data
fprintf(fileID,'Element Data: \r\n \r\n');
fprintf(fileID,'%6s %15s %12s %12s\r\n','elem','nodes','E-modulus','area');

elementArray = getAllElements(Model);

for ii = 1:1:size(Model.getAllElements,2);
elementNodes = getId(getNodes(elementArray(ii)));
beginNode=elementNodes(1,1);
endNode=elementNodes(1,2);
id = getId(elementArray(ii));
A = [getId(elementArray(ii)); beginNode; endNode; ...
    getParameterValue(getMaterial(elementArray(ii)), 'YOUNGS_MODULUS'); ...
    getCrossSectionArea(elementArray(ii))];
fprintf(fileID,'%6.0f \t [%5.0f,%5.0f] %12.2f %12.2f\r\n',A);
end

fprintf(fileID,'\r\n \r\n');

% Print DOF Activity
fprintf(fileID,'DOF Activity (DOFs fixed, DOF load): \r\n \r\n');
fprintf(fileID,'%6s %6s %6s %6s %10s %10s %10s  %10s \r\n',...
    'node','x-tag','y-tag','z-tag','x-load','y-load','z-load');
fprintf(fileID,'\r\n');

for ii = 1:1:size(Model.getAllNodes,2)
    
    
    dofArray = getDofArray(nodeArray(ii));
    
    for jj = 1:length(dofArray)
        fixedDof(jj) = isFixed(dofArray(jj));
        if fixedDof(jj) == 1
            loadDof(jj) = 0;
        else
            loadDof(jj) = getValue(dofArray(jj));
        end
    end
    
    A = [getId(nodeArray(ii)); fixedDof(1); fixedDof(2); fixedDof(3);...
        loadDof(1); loadDof(2); loadDof(3)];
    fprintf(fileID,'%6.0f %6.0f %6.0f %6.0f %10.2f %10.2f %10.2f\r\n',A);
end

fclose(fileID);

end

