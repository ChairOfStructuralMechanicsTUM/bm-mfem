function [ outputFile ] = createOutputFile( Model ) 
%CREATEOUTPUTFILE Creates an output file 
%   including information about geometry, deformation and stresses


% Create TextFile
fileID = fopen('femCalculations.txt','w');

% Specify output step
step = 1;


%% Print Node Coordinates

fprintf(fileID,'Node Coordinates: \r\n \r\n');

nodeArray = getAllNodes(Model);

ii = 1:1:size(Model.getAllNodes,2);
A = [getId(nodeArray(ii)); getX(nodeArray(ii)); getY(nodeArray(ii)); getZ(nodeArray(ii))];

fprintf(fileID,'%6s %12s %12s %12s\r\n','node','x-coord','y-coord','z-coord');
fprintf(fileID,'%6.0f %12.2f %12.2f %12.2f\r\n',A);

fprintf(fileID,'\r\n \r\n');


%% Print Element Data

fprintf(fileID,'Element Data: \r\n \r\n');
fprintf(fileID,'%6s %15s %12s %12s\r\n','elem','nodes','E-modulus','area');

elementArray = getAllElements(Model);

for ii = 1:1:size(Model.getAllElements,2)
    elementNodes = getId(getNodes(elementArray(ii)));
    beginNode=elementNodes(1,1);
    endNode=elementNodes(1,2);
    % id = getId(elementArray(ii));
    A = [getId(elementArray(ii)); beginNode; endNode; ...
        getPropertyValue(elementArray(ii), 'YOUNGS_MODULUS'); ...
        getPropertyValue(elementArray(ii), 'CROSS_SECTION')];
    fprintf(fileID,'%6.0f \t [%5.0f,%5.0f] %12.2f %12.2f\r\n',A);
end

fprintf(fileID,'\r\n \r\n');


%% Print DOF Activity

fprintf(fileID,'DOF Activity (DOFs fixed, DOF load): \r\n \r\n');
fprintf(fileID,'%6s %6s %6s %6s %10s %10s %10s  %10s \r\n',...
    'node','x-tag','y-tag','z-tag','x-load','y-load','z-load');
fprintf(fileID,'\r\n');

for ii = 1:1:size(Model.getAllNodes,2)

    dofArray = getDofArray(nodeArray(ii));

    for jj = 1:length(dofArray)
        fixedDof(jj) = isFixed(dofArray(jj));
    end
    
    for jj = 1:length(dofArray)
        loadDof(jj)=zeros;
        if isempty(getDofLoad(dofArray(jj)))
            loadDof(jj)=0;
        else
            
            loadDof(jj) = getDofLoad(dofArray(jj));
        end
    end
    
    A = [getId(nodeArray(ii)); fixedDof(1); fixedDof(2); fixedDof(3);...
        loadDof(1); loadDof(2); loadDof(3)];
    fprintf(fileID,'%6.0f %6.0f %6.0f %6.0f %10.2f %10.2f %10.2f\r\n',A);
end

fprintf(fileID,'\r\n \r\n');


%% Print Node Displacements

fprintf(fileID,'Node Displacements in Step %d: \r\n \r\n', step);
fprintf(fileID,'%6s %14s %14s %14s \r\n','node','x-displ','y-displ','z-displ');
fprintf(fileID,'\r\n');

for ii = 1:1:size(Model.getAllNodes,2)

    dofArray = getDofArray(nodeArray(ii));

    for jj = 1:length(dofArray)
        fixedDof(jj) = isFixed(dofArray(jj));
        if fixedDof(jj) == 1
            valueDof(jj) = 0;
        else
            valueDof(jj) = getValue(dofArray(jj), step);
        end
    end
    
    A = [getId(nodeArray(ii)); valueDof(1); valueDof(2); valueDof(3)];
    fprintf(fileID,'%6.0f %14.6f %14.6f %14.6f\r\n',A);
end

fprintf(fileID,'\r\n \r\n');


%% Print Node Forces inclusive Reactions

fprintf(fileID,'Node Forces (inclusive Reactions) in Step %d: \r\n \r\n', step);
fprintf(fileID,'%6s %14s %14s %14s \r\n','node','x-force','y-force','z-force');
fprintf(fileID,'\r\n');
%TODO: ugly workaround
solver = SimpleSolvingStrategy(Model);

for ii = 1:1:size(Model.getAllNodes,2)

    dofArray = getDofArray(nodeArray(ii));
    nodeForces = solver.getNodalForces(step);
     
  
    for jj = 1:length(dofArray)
        dofForce(jj) = zeros;
        dofForce(jj) = nodeForces((ii-1)*length(dofArray)+jj);
    end
       

    A = [getId(nodeArray(ii)); dofForce(1); dofForce(2); dofForce(3)];
    fprintf(fileID,'%6.0f %14.6f %14.6f %14.6f\r\n',A);
end

fprintf(fileID,'\r\n \r\n');


%% Print Internal Element Forces and Stresses

fprintf(fileID,'Internal Element Stresses and Forces in Step %d: \r\n \r\n', step);
fprintf(fileID,'%6s %18s %18s  \r\n','elem','axial stress','axial force');
fprintf(fileID,'\r\n');

elementArray = getAllElements(Model);

for ii = 1:1:size(Model.getAllElements,2)
    
    elementStress(ii) = computeElementStress(elementArray(ii), step);
    elementForce(ii) = elementStress(ii) * getPropertyValue(elementArray(ii), 'CROSS_SECTION');
    
    
    A = [getId(elementArray(ii)); elementStress(ii); elementForce(ii)];
    fprintf(fileID,'%6.0f %18.6f %18.6f\r\n',A);
end


fprintf(fileID,'\r\n \r\n');


%% Close TextFile
fclose(fileID);

end

