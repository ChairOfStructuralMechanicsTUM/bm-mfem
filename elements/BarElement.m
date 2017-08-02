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
            if barElement.length == 0
                barElement.length = computeLength(barElement.nodeArray(1).getCoords, ...
                    barElement.nodeArray(2).getCoords);
            end
            len = barElement.length;
        end
        
        function area = getCrossSectionArea(barElement)
            area = barElement.crossSectionArea;
        end
        
        % setter functions
        function setCrossSectionArea(barElement, area)
            barElement.crossSectionArea = area;
        end
        
        % member functions
        function update(barElement)
            barElement.length = computeLength(barElement.nodeArray(1).getCoords, ...
                    barElement.nodeArray(2).getCoords);
        end
        
        function c = barycenter(barElement)
            nodes = barElement.getNodes;
            c = (nodes(1).getCoords + nodes(2).getCoords) ./ 2;
        end
        
        function pl = draw(barElement)
            pl = line(barElement.nodeArray.getX, barElement.nodeArray.getY);
        end
        
        function pl = drawDeformed(barElement, step, scaling)
            pl = line(barElement.nodeArray.getX' + scaling * barElement.nodeArray.getDofValue('DISPLACEMENT_X', step), ...
                barElement.nodeArray.getY' + scaling * barElement.nodeArray.getDofValue('DISPLACEMENT_Y', step));
        end
        
    end
    
    methods (Access = protected)
       function cp = copyElement(obj)
           cp = copyElement@Element(obj);
           cp.crossSectionArea = obj.crossSectionArea;
           cp.length = obj.length;
        end 
    end
    
end

