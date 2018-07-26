classdef FetiDPModel < FemModel
    properties (Access= private)
    end
    
    methods
        % call constructor of superclass FemModel
        function fetidpModel=FetiDPModel(~)
            fetidpModel=fetidpModel@FemModel(nodeArray, elementArray, femModelParts); 
        end
    end
end