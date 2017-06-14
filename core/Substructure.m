classdef Substructure < FemModel
    %Substructure The class of the individual substructure models
    %   This class is a substructure of the mesh in a FemModel. It can be 
    %   used to create substructes which are later used for FETI methods
    
    properties (Access = private)
        nodesIntf
    end 
    
    methods
        %constructor
        function substructur = Substructure(nodeArray, elementArray, nodesIntf, femModelParts)
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                super_args = {nodeArray, elementArray};
            elseif nargin == 3
                super_args = {nodeArray, elementArray};
            elseif nargin == 4
                super_args = {nodeArray, elementArray, femModelParts};
            end
            substructur@FemModel(super_args{:});
            substructur.nodesIntf = nodesIntf;
        end    
        
        %getter
        function nodesIntf = getNodesIntf(substructur)
            nodesIntf = substructur.nodesIntf;
        end
    end
        
        
end
    


