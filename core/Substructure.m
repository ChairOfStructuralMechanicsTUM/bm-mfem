classdef Substructure < FemModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
   interfaceNodes
    end
    
    methods
        function substructure = Substructure(nodeArray,elementArray,femModelParts)
            if nargin == 0
               super_args = {};
           elseif nargin == 2
              super_args = {nodeArray, elementArray};
            elseif nargin == 3
               super_args = {nodeArray, elementArray, femModelParts};
            end
        substructure@FemModel(super_args{:});
        
        end
    
        
    end
    
end

