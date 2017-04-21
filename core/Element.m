classdef (Abstract) Element < handle
    %ELEMENT The element class
    %   Abstract base class for all element implementations
    
    properties (Access = private)
        id
        material
    end
    properties (Access = protected)
        nodeArray
    end
    
    methods
        % constructor
        function element = Element(id, material)
            if (nargin > 0)
                element.id = id;
                if (isa(material,'Material'))
                    element.material = material;
                else
                    error('problem with the material in element %d', id);
                end
                element.nodeArray = {};
            end
        end
        
        % getter functions
        function id = getId(element)
            id = element.id;
        end
        
        function material = getMaterial(element)
            material = element.material;
        end
        
        function nodeArray = getNodeArray(element)
            nodeArray = element.nodeArray;
        end
        
    end
    
    methods (Access = protected)
        
        function addDofs(element, dofNames)
            for itNode = 1:length(element.nodeArray)
                nodalDofs(1, length(dofNames)) = Dof;
                for itDof = 1:length(dofNames)
                    newDof = Dof(element.nodeArray(itNode),0.0,dofNames(itDof));
                    nodalDofs(itDof) = newDof;
                end
                element.nodeArray(itNode).setDofArray(nodalDofs);
            end
        end
        
    end
    
end

