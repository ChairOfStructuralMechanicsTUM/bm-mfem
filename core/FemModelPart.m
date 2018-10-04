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
        
        function n = getNodeById(obj, id)
            n = Node.empty;
            if ~ isempty(obj.nodeArray)
                nodeIds = obj.nodeArray.getId();
                n = obj.nodeArray(id==nodeIds);
            end
            if isempty(n)
                msg = ['FemModelPart: Node with id ', num2str(id), ...
                    ' not found in modelpart ', obj.name];
                e = MException('MATLAB:bm_mfem:nodeNotFound',msg);
                throw(e);
            end
        end
        
        function e = getElements(obj)
            e = obj.elementArray;
        end
        
        function e = getElementById(obj, id)
            e = Element.empty;
            if ~ isempty(obj.elementArray)
                elementIds = obj.elementArray.getId();
                e = obj.elementArray(id==elementIds);
            end
            if isempty(e)
                msg = ['FemModelPart: Element with id ', num2str(id), ...
                    ' not found in modelpart ', obj.name];
                err = MException('MATLAB:bm_mfem:elementNotFound',msg);
                throw(err);
            end
        end
        
        function d = getDofArray(obj)
            d = obj.dofArray;
        end
        
        function [free, fix] = getDofConstraints(obj)
            free = obj.freeDofs;
            fix = obj.fixedDofs;
        end
        
        function p = getParentFemModel(obj)
            p = obj.parentFemModel;
        end
        
        function addNode(obj, node)
            if isa(node,'Node')
                obj.nodeArray = [obj.nodeArray node];
            else
                obj.nodeArray = [obj.nodeArray obj.parentFemModel.getNode(node)];
            end
        end
        
        function addElement(obj, element)
            if isa(element,'Element')
                obj.elementArray = [obj.elementArray element];
            else
                obj.elementArray = [obj.elementArray obj.parentFemModel.getElement(element)];
            end
        end
        
        function n = addNewNode(obj, id, x, y, z)
            if ~isempty(obj.parentFemModel)
                n = obj.parentFemModel.addNewNode(id, x, y, z);
                obj.nodeArray = [obj.nodeArray n];
            else
                msg = 'FemModelPart: No parent model defined';
                e = MException('MATLAB:bm_mfem:noParentModel',msg);
                throw(e);
            end
        end
        
        function e = addNewElement(obj, elementName, id, nodes, props)
            if ~isempty(obj.parentFemModel)
                if nargin == 4
                    e = obj.parentFemModel.addNewElement(elementName, id, nodes);
                elseif nargin == 5
                    e = obj.parentFemModel.addNewElement(elementName, id, nodes, props);
                else
                    msg = 'FemModelPart: Wrong number of input arguments';
                    err = MException('MATLAB:bm_mfem:invalidArguments',msg);
                    throw(err);
                end
                obj.elementArray = [obj.elementArray e];
            else
                msg = 'FemModelPart: No parent model defined';
                err = MException('MATLAB:bm_mfem:noParentModel',msg);
                throw(err);
            end
        end
    end
    
end

