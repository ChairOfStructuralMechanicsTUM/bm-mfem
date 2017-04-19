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
           element.nodeArray = {};
           element.dofArray = {};
        end
        
        % getter functions
        function material = getMaterial(element)
            material = element.material;
        end
        
        function nodeArray = getNodeArray(element)
            nodeArray = element.nodeArray;
        end
        
        function dofArray = getDofArray(element)
            dofArray = element.dofArray;
        end
        
    end
    
    methods (Access = protected)
        
        function addDofs(element, dofNames)
            for itNode = 1:length(element.nodeArray)
                for itDof = 1:length(dofNames)
                    newDof = Dof(element.nodeArray(itNode),0.0,dofNames(itDof));
                    element.dofArray = [element.dofArray, newDof];
                end
            end
        end
        
    end
    
end

