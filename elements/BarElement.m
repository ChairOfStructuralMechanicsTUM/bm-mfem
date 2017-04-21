classdef (Abstract) BarElement < Element
    %BAR_ELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        crossSectionArea
        length
    end
    
    methods
        % constructor
        function barElement = BarElement(id, material, crossSectionArea)
           barElement@Element(id, material);
           
           barElement.crossSectionArea = crossSectionArea;
           
        end
        
        % getter functions
        function len = getLength(barElement)
            len = barElement.length;
        end
        
        function area = getCrossSectionArea(barElement)
            area = barElement.crossSectionArea;
        end
        
        % member functions
        
    end
    
end

