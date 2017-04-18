classdef Element < handle
    %ELEMENT The element class
    %   Detailed explanation goes here
    
    properties (Access = private)
        id
        material
    end
    properties (Access = protected)
        nodeArray
        dofArray
    end
    
    methods
        % constructor
        function element = Element(id, material)
           element.id = id; 
           if (isa(material,'Material'))
               element.material = material;
           else
               error('problem with the material in element %d', id);
           end
        end
        
        % getter functions
        function material = getMaterial(element)
            material = element.material;
        end
        
        function nodeArray = getNodeArray(element)
            nodeArray = element.nodeArray;
        end
        
    end
    
end

