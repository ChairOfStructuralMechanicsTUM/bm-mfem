classdef FemModelPart < handle
    %FEMMODELPART Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        nodeArray
        elementArray
        dofArray
        fixedDofs
        freeDofs
        parentFemModel
    end
    
    methods
        
        function obj = FemModelPart(name, nodes, elements, parentFemModel)
            obj.name = name;
            obj.nodeArray = nodes;
            obj.elementArray = elements;
            if nargin > 3
                obj.parentFemModel = parentFemModel;
            end
        end
        
        function initialize(obj)
            if ~isempty(obj.nodeArray)
                tmp = arrayfun(@(node) node.getDofArray, obj.nodeArray, 'UniformOutput', false);
                obj.dofArray = [tmp{:}];

                fixed = obj.dofArray.isFixed();
                obj.fixedDofs = obj.dofArray(fixed);
                obj.freeDofs = obj.dofArray(~fixed);
            end
        end
        
        function n = getName(obj)
            n = obj.name;
        end
        
        function n = getNodes(obj)
            n = obj.nodeArray;
        end
        
        function e = getElements(obj)
            e = obj.elementArray;
        end
        
        function d = getDofArray(obj)
            d = obj.dofArray;
        end
        
        function [fix, free] = getDofConstraints(obj)
            fix = obj.fixedDofs;
            free = obj.freeDofs;
        end
        
        function p = getParentFemModel(obj)
            p = obj.parentFemModel;
        end
        
        function addNode(obj, node)
            if isa(node,'Node')
                obj.nodeArray = [obj.nodeArray node];
            else
                obj.nodeArray = [obj.nodeArray obj.getNode(node)];
            end
        end
        
        function addElement(obj, element)
            if isa(element,'Element')
                obj.elementArray = [obj.elementArray element];
            else
                obj.elementArray = [obj.elementArray obj.getElement(element)];
            end
        end
    end
    
end

