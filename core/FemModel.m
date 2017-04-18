classdef FemModel < handle
    %FEMMODEL The class of the complete model
    %   This class keeps track over all entities in the model 
    
    properties (Access = private)
        nodeArray
        elementArray
        
        
    end
    
    methods
        function loadFemModel(femModel, nodeArray, elementArray)
           
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            if (hasDuplicates(nodeIds))
                error('problem with node ids');
            else
                femModel.nodeArray = nodeArray;
            end
            
            elementIds = arrayfun(@(node) node.getId, nodeArray);
            if (hasDuplicates(elementIds))
                error('problem with element ids');
            else
                femModel.elementArray = elementArray;
            end
            
        end
    end
    
end

