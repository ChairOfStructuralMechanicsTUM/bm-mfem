classdef FemModel < handle
    %FEMMODEL The class of the complete model
    %   femModel = FemModel(nodeArray, elementArray, femModelParts)
    %   This class keeps track over all entities in the model
    
    properties (Access = public) % currently changed from private to public
        nodeArray
        elementArray
        %         dofArray = {}
        femModelParts = containers.Map
        
    end
    
    properties (Access = public)
        dofArray = {}
    end
    
    methods
        % constructor
        function femModel = FemModel(nodeArray, elementArray, femModelParts)
            % dim=2; % dim=1 for X; =2 for Y; =3 for Z 
            
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            if (hasDuplicates(nodeIds))
                error('problem with node ids');
            else
                femModel.nodeArray = nodeArray;
            end
            
            femModel.dofArray = arrayfun(@(node) node.getDofArray, nodeArray, 'UniformOutput', false);
            femModel.dofArray = [femModel.dofArray{:}];
            
            elementIds = arrayfun(@(element) element.getId, elementArray);
            if (hasDuplicates(elementIds))
                error('problem with element ids');
            else
                femModel.elementArray = elementArray;
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
            if nargin == 4
                properties = PropertyContainer();
            end
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

