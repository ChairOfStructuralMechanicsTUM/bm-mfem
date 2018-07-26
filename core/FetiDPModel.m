classdef FetiDPModel < FemModel
    properties (Access= private)
        %FetiDPModel inherits all properties of FemModel
    end
    
    methods
        % call constructor of superclass FemModel
        function fetidpModel=FetiDPModel(~)
            fetidpModel=fetidpModel@FemModel(nodeArray, elementArray, femModelParts); 
        end
        % FetiDP substructuring
        function substructure(~)
        end
    end
end