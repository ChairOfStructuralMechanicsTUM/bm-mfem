function [displacements, displacementsIntf] = orderDisplacements(substructures)
    %function which orders the degrees of freedom by coordinates and
    %gives their displacements
    
    nodes = [];
    %get all nodes of all substructures
    for numSubs = 1:length(substructures)
        nodes = [nodes substructures(numSubs).getAllNodes];
    end
    
    %find dimensions
    dim = Substructure.findDim(substructures(1,1));
    
    switch dim
        case 3
            %order nodes by z-direction
            [~,idx] = sort([nodes.getX]);
            nodes = nodes(idx);
            
            ii = 2;
            n = 0;
            %order nodes by x-direction
            while ii <= length(nodes)
                if nodes(ii).getZ < nodes(ii-1).getZ 
                    temp = nodes(ii);
                    nodes(ii) = nodes(ii-1);
                    nodes(ii-1) = temp;
                    ii = ii-n;
                    n = 0;
                else
                    ii = ii+1;
                    n = n+1;
                end
            end
     
            
            ii = 2;
            n = 0;
            %order nodes by y-direction
            while ii <= length(nodes)
                if nodes(ii).getY < nodes(ii-1).getY && nodes(ii).getZ == nodes(ii-1).getZ 
                    temp = nodes(ii);
                    nodes(ii) = nodes(ii-1);
                    nodes(ii-1) = temp;
                    ii = ii-n;
                    n = 0;
                elseif nodes(ii).getZ ~= nodes(ii-1).getZ
                    ii = ii+1;
                    n = 0;
                else
                    ii = ii+1;
                    n = n+1;
                end
            end
            
        case 2
            %order nodes by x-direction
            [~,idx] = sort([nodes.getX]);
            nodes = nodes(idx);
            
            ii = 2;
            n = 0;
            %order nodes by y-direction
            while ii <= length(nodes)
                if nodes(ii).getY < nodes(ii-1).getY && nodes(ii).getX == nodes(ii-1).getX
                    temp = nodes(ii);
                    nodes(ii) = nodes(ii-1);
                    nodes(ii-1) = temp;
                    ii = ii-n;
                    n = 0;
                elseif nodes(ii).getX ~= nodes(ii-1).getX
                    ii = ii+1;
                    n = 0;
                else
                    ii = ii+1;
                    n = n+1;
                end
            end
            
        case 1
            [~,idx] = sort([nodes.getX]);
            nodes = nodes(idx);  
    end
    
    %filter out the interface dofs
    numNodes = 2;
    intf = [];
    while numNodes <= length(nodes)
        if nodes(numNodes).getCoords == nodes(numNodes-1).getCoords
            intf = [intf nodes(numNodes)];
            nodes(numNodes) = [];
        else
            numNodes = numNodes+1;
        end
    end
    
    %find displacements
    displacements = zeros(dim*length(nodes),1);
    for numNodes = 1:length(nodes)
        dofs = nodes(numNodes).getDofArray;
        for dof = 1:dim
            displacements(numNodes*dim+dof-dim,1) = dofs(dof).getValue;
        end
    end
    %find interface displacements
    displacementsIntf = zeros(2*length(intf),1);
    for numNodes = 1:length(intf)
        dofs = intf(numNodes).getDofArray;
        for dof = 1:dim
            displacementsIntf(numNodes*dim+dof-dim,1) = dofs(dof).getValue;
        end
    end
end

