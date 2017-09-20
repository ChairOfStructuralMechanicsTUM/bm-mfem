classdef Substructure < FemModel
    %Substructure The class of the individual substructure models
    %   This class is a substructure of the mesh in a FemModel. It can be 
    %   used to create substructes which are later used for FETI methods
    
    properties (Access = private)
        %nodes at the interface
        nodesIntf
    end 
    
    methods
        %constructor
        function substructure = Substructure(nodeArray, elementArray, nodesIntf, idt, femModelParts)
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                super_args = {nodeArray, elementArray};
            elseif nargin == 4
                super_args = {nodeArray, elementArray};
            elseif nargin == 5
                super_args = {nodeArray, elementArray, femModelParts};
            end
            substructure@FemModel(super_args{:});
            substructure.nodesIntf = [nodesIntf];
            substructure.nodesIntf.getId;

            Substructure.checkNodeIds(nodeArray, idt);
        end    
        
        %getter
        function nodesIntf = getNodesIntf(substructur)
            nodesIntf = substructur.nodesIntf;
        end
    end
    
    methods (Static)
        function dim = findDim(substructures)
            %find dimension of the problem: 1D, 2D, 3D
            nodes = substructures(1,1).getAllNodes;
            count = 0;
            for numNodes = 1:length(nodes)
                dofs = nodes(numNodes).getDofArray;
                for dof = 1:length(dofs)
                    if ~dofs(dof).isFixed
                        count = count+1;
                    end
                end
            end
            
            n = count/numel(nodes);
            %assign dimension
            switch (1)
                case n > 2
                    dim = 3;
                case n > 1
                    dim = 2;
                case n > 0
                    dim = 1;
            end
        end
        
        function checkNodeIds(nodeArray, idt)
            %functions which checks whether the node Ids in a substructure
            %are all directly increasing (just one higher than previous
            %node)
            for ii = 1:(length(nodeArray)-1)
                if nodeArray(ii+1).getId - nodeArray(ii).getId > 1
                    Substructure.changeNodeIds(nodeArray,idt);
                end
            end
        end
        
        function changeNodeIds(nodeArray, idt)
            %functions changes the Ids of the nodes in a substructure to be
            %directly increasing. otherwise computation of the stiffness
            %matrix is not possible in the simple assembler.
            
            id = idt.getNodeId;
            for ii = 1:length(nodeArray)
                nodeArray(ii).setId(id+ii);
            end
            idt.setNodeId(id+ii);
        end
    end 
end
    


