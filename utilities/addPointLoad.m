function addPointLoad( nodeArray, modulus, direction )
%ADDPOINTLOAD Adds a point load on an array of nodes
%   nodeArray : nodes where the load is imposed
%   modulus : load modulus
%   direction : load direction with [x, y(, z)]

for itNode = 1:length(nodeArray)
    currentNode = nodeArray(itNode);
   
    currentDofArray = getDofArray(currentNode);
    
    if length(currentDofArray) ~= length(direction)
        error('the dimensions of the dofs and the load are not the same')
    end
    
    direction = direction./norm(direction);
    
    for itDof = 1:length(currentDofArray)
        currentDofArray(itDof).setLoad(modulus * direction(itDof));
 %       setLoad(currentDof(itDof), modulus * direction(1)); 
    end
   
    
   
%    currentNode.setDofValue('DISPLACEMENT_X', modulus * direction(1));
%    currentNode.setDofValue('DISPLACEMENT_Y', modulus * direction(2));
%
%    if length(direction) == 3
%        currentNode.setDofValue('DISPLACEMENT_Z', modulus * direction(3));
%
%    end

end
end

