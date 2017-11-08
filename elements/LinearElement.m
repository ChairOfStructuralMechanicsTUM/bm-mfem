classdef (Abstract) LinearElement < Element
    %BARELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        crossSectionArea
        length
    end
    
    methods
        % constructor
        function e = LinearElement(id, requiredProperties)
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                super_args = {id, requiredProperties};
            end
            
            e@Element(super_args{:});            
        end
        
        % getter functions
        function len = getLength(element)
            if element.length == 0
                element.length = computeLength(element.nodeArray(1).getCoords, ...
                    element.nodeArray(2).getCoords);
            end
            len = element.length;
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

