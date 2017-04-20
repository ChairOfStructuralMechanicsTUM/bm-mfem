function addPointLoad( nodeArray, modulus, direction )
%ADDPOINTLOAD Adds a point load on an array of nodes
%   nodeArray : nodes where the load is imposed
%   modulus : load modulus
%   direction : load direction with [x, y(, z)]
for itNode = 1:length(nodeArray)
   currentNode = nodeArray(itNode);
   currentNode.setDofValue('DISPLACEMENT_X', modulus * direction(1));
   currentNode.setDofValue('DISPLACEMENT_Y', modulus * direction(2));
   if length(direction) == 3
       currentNode.setDofValue('DISPLACEMENT_Z', modulus * direction(3));
   end
end

end

