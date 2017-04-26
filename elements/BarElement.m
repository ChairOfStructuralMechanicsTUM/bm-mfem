classdef (Abstract) BarElement < Element
    %BARELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        crossSectionArea
        length
    end
    
    methods
        % constructor
        function barElement = BarElement(id, material, crossSectionArea)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id; material};
            end
            
            barElement@Element(super_args{:});
            
            if nargin > 0
                barElement.crossSectionArea = crossSectionArea;
            end
            
        end
        
        % getter functions
        function len = getLength(barElement)
            len = barElement.length;
        end
        
        function area = getCrossSectionArea(barElement)
            area = barElement.crossSectionArea;
        end
        
        % member functions
        function c = barycenter(barElement)
            nodes = barElement.getNodes;
            c = (nodes(1).getCoords + nodes(2).getCoords) ./ 2;
        end
        
    end
    
end

