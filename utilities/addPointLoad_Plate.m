function addPointLoad( nodeArray, modulus, direction )
%ADDPOINTLOAD Adds a point load on an array of nodes
%   nodeArray : nodes where the load is imposed
%   modulus : load modulus
%   direction : load direction with [x, y(, z)]

direction = direction ./ norm(direction);
load = direction .* modulus;

nodeArray.setDofLoad('DISPLACEMENT_Z', load(1));
nodeArray.setDofLoad('ROTATION_X', load(2));
if length(load) == 3
    nodeArray.setDofLoad('ROTATION_Y', load(3));
end

end



