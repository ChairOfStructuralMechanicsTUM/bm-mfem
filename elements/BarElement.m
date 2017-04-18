classdef BarElement < Element
    %BAR_ELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        crossSectionArea
    end
    
    methods
        % constructor
        function barElement = BarElement(id, nodeArray, material, crossSectionArea)
           barElement@Element(id, material);
           
           if (length(nodeArray) == 2 && isa(nodeArray,'Node'))
               barElement.nodeArray = nodeArray;
           else
               error('problem with the nodes in element %d', id);
           end
           
           barElement.crossSectionArea = crossSectionArea;
           
        end
        
        % getter functions
        
    end
    
end

