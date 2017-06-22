classdef FemModel < handle
    %FEMMODEL The class of the complete model
    %   femModel = FemModel(nodeArray, elementArray, femModelParts)
    %   This class keeps track over all entities in the model
    
    properties (Access = private)
        nodeArray = Node.empty
        elementArray = Element.empty
        femModelParts = containers.Map
        dofArray = {}
        initialized = false
    end
    
    methods
        % constructor
        function femModel = FemModel(nodeArray, elementArray, femModelParts)
            
            if nargin > 0
                nodeIds = arrayfun(@(node) node.getId, nodeArray);
                duplicated_nodes = findDuplicates(nodeIds);
                if ~ isempty(duplicated_nodes)
                    error('multiple nodes with id %d exist',duplicated_nodes);
                else
                    femModel.nodeArray = nodeArray;
                end
                
                elementIds = arrayfun(@(element) element.getId, elementArray);
                duplicated_elements = findDuplicates(elementIds);
                if ~ isempty(duplicated_elements)
                    error('multiple elements with id %d exist',duplicated_elements);
                else
                    femModel.elementArray = elementArray;
                end
            end
            
            if nargin == 3
                femModel.femModelParts = femModelParts;
            end
            
        end
        
        % getter functions
        function nodeArray = getAllNodes(femModel)
            nodeArray = femModel.nodeArray;
        end
        
        function elementArray = getAllElements(femModel)
            elementArray = femModel.elementArray;
        end
        
        function dofArray = getDofArray(femModel)
            if ~ femModel.initialized
                femModel.initialize;
            end
            
            dofArray = femModel.dofArray;
        end
        
        function node = getNode(femModel, id)
            node = femModel.nodeArray(id);
        end
        
        function nodes = getNodes(femModel, ids)
            nodes = Node.empty;
            for ii = 1:length(ids)
                nodes(ii) = femModel.nodeArray(ids(ii));
            end
        end
        
        function element = getElement(femModel, id)
            element = femModel.elementArray(id);
        end
        
        function elements = getElements(femModel, ids)
            elements = Element.empty;
            for ii = 1:length(ids)
                elements(ii) = femModel.elementArray(ids(ii));
            end
        end
        
        function femModelParts = getAllModelParts(femModel)
            femModelParts = femModel.femModelParts;
        end
        
        function modelPart = getModelPart(femModel, name)
            modelPart = femModel.femModelParts(name);
        end
        
        % setter functions
        function addModelPart(femModel, name, entityArray)
            femModel.femModelParts(name) = entityArray;
        end
        
        % member functions
        function initialize(femModel)
            femModel.dofArray = arrayfun(@(node) node.getDofArray, femModel.nodeArray, 'UniformOutput', false)';
            femModel.dofArray = [femModel.dofArray{:}];
            femModel.dofArray = reshape(femModel.dofArray,1,size(femModel.dofArray,1)*size(femModel.dofArray,2));
            for ii = 1:length(femModel.dofArray)
               femModel.dofArray(ii).setId(ii); 
            end
            
            femModel.initialized = true;
        end
        
        function node = addNewNode(femModel, id, x, y, z)
            nodeIds = arrayfun(@(node) node.getId, femModel.nodeArray);
            if any(id == nodeIds)
                error('a node with id %d already exists in the model', id)
            end
            
            if nargin == 4
                node = Node(id, x, y);
                femModel.nodeArray(id) = node;
            elseif nargin == 5
                node = Node(id, x, y, z);
                femModel.nodeArray(id) = node;
            else
                error('wrong input parameters')
            end
        end
        
        function element = addNewElement(femModel, elementName, id, nodes, properties)
            elementIds = arrayfun(@(element) element.getId, femModel.elementArray);
            if any(id == elementIds)
                error('an element with id %d already exists in the model', id)
            end
            
            if ~ isa(nodes,'Node')
                nodes = femModel.getNodes(nodes);
            end
            
            switch elementName
                case 'BarElement2d2n'
                    crossArea = properties.getValue('CROSS_SECTION');
                    element = BarElement2d2n(id, nodes, properties, crossArea);
                case 'BarElement3d2n'
                    crossArea = properties.getValue('CROSS_SECTION');
                    element = BarElement3d2n(id, nodes, properties, crossArea);
                case 'ConcentratedMassElement3d1n'
                    element = ConcentratedMassElement3d1n(id, nodes, properties);
                case 'SpringDamperElement3d2n'
                    element = SpringDamperElement3d2n(id, nodes, properties);
                    
                otherwise
                    error('unknown element %s',elementName)
            end %switch
            
            femModel.elementArray(id) = element;
            
        end
        
    end
    
end

