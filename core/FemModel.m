classdef FemModel < handle
    %FEMMODEL The class of the complete model
    %   femModel = FemModel(nodeArray, elementArray, femModelParts)
    %   This class keeps track over all entities in the model
    
    properties (Access = private)
        nodeArray = Node.empty
        elementArray = Element.empty
        femModelParts
        dofArray
        initialized = false
        fixedDofs
        freeDofs
        mProperties
    end
    
    methods
        % constructor
        function obj = FemModel(nodeArray, elementArray, femModelParts)
            
            obj.mProperties = PropertyContainer();
            obj.femModelParts = containers.Map;
            
            if nargin > 0
                nodeIds = arrayfun(@(node) node.getId, nodeArray);
                duplicated_nodes = findDuplicates(nodeIds);
                if ~ isempty(duplicated_nodes)
                    error('multiple nodes with id %d exist',duplicated_nodes);
                else
                    obj.nodeArray = nodeArray;
                end
                
                elementIds = arrayfun(@(element) element.getId, elementArray);
                duplicated_elements = findDuplicates(elementIds);
                if ~ isempty(duplicated_elements)
                    error('multiple elements with id %d exist',duplicated_elements);
                else
                    obj.elementArray = elementArray;
                end
            end
            
            if nargin == 3
                obj.femModelParts = femModelParts;
            end
            
        end
        
        % getter functions
        function init = isInitialized(obj)
            init = obj.initialized;
        end
        
        function nodeArray = getAllNodes(obj)
            nodeArray = obj.nodeArray;
        end
        
        function elementArray = getAllElements(obj)
            elementArray = obj.elementArray;
        end
        
        function dofArray = getDofArray(obj)
            if ~ obj.initialized; obj.initialize; end
            dofArray = obj.dofArray;
        end
        
        function [freeDofs, fixedDofs] = getDofConstraints(obj)
            if ~ obj.initialized; obj.initialize; end
            freeDofs = obj.freeDofs;
            fixedDofs = obj.fixedDofs;
        end
        
        function node = getNode(obj, id)
            node = obj.nodeArray(id);
        end
        
        function nodes = getNodes(obj, ids)
            nodes = Node.empty;
            for ii = 1:length(ids)
                nodes(ii) = obj.nodeArray(ids(ii));
            end
        end
        
        function element = getElement(obj, id)
            element = obj.elementArray(id);
        end
        
        function elements = getElements(obj, ids)
            elements = Element.empty;
            for ii = 1:length(ids)
                elements(ii) = obj.elementArray(ids(ii));
            end
        end
        
        function mp = getAllModelParts(obj)
            %GETALLMODELPARTS returns the names of all model parts
            mp = obj.femModelParts.keys;
        end
        
        function mp = getModelPart(obj, name)
        %GETMODELPART returns a model part
        %   mp = GETMODELPART(name)
        %   see also FEMMODELPART
            mp = obj.femModelParts(name);
        end
        
        function mp = getMainModelPart(obj)
        %GETMAINMODELPART returns the main model part
            if ~ obj.initialized; obj.initialize; end
            mp = obj.femModelParts('MainModelPart');
        end
        
        function mProperties = getProperties(femModel)
            mProperties = femModel.mProperties;
        end
        
        % setter functions
        function addNewModelPart(obj, name, nodeIds, elementIds)
        %ADDNEWMODELPART add a new model part
        %   ADDNEWMODELPART(name, nodeIds, elementIds) adds a modelpart to the 
        %   FemModel containing nodes with nodeIds and/or elements with
        %   elementIds.
        %    
        %   see also FEMMODELPART
        
            nodes = obj.getNodes(nodeIds);
            elements = obj.getElements(elementIds);
            
            obj.femModelParts(name) = FemModelPart(name, ...
                nodes, elements, obj);
            
            obj.initialized = false;
        end
        
        % member functions
        function initialize(obj)
        %INITIALIZE makes the model ready for computation by performing the
        %   following tasks:
        %       (1) check if all elements are initialized correctly
        %       (2) assign unique ids to the dofs
        %       (3) determine free and fixed dofs
        %       (4) sets the step to 1
        %       (5) add all nodes and elements to the main model part
        %       (6) initialize the modelparts
            if obj.initialized
                return;
            end
%             
%             % (1) check for double ids
%             
%             nodeIds = arrayfun(@(node) node.getId, obj.nodeArray);
%             [~, uids] = unique(nodeIds);
%             if length(nodeIds) ~= length(uids)
%                 duplicate_ids = setdiff(1:length(nodeIds), uids);
%                 msg = ['FemModel: Multiple nodes with id(s) \"', ...
%                     nodeIds(duplicate_ids), '\" exist'];
%                 e = MException('MATLAB:bm_mfem:duplicateId',msg);
%                 throw(e);
%             end
            
            % (1) check all elements
            arrayfun(@(e) e.check(), obj.elementArray);
            
            % (2) assign dof ids
            obj.dofArray = arrayfun(@(node) node.getDofArray, obj.nodeArray, 'UniformOutput', false);
            obj.dofArray = [obj.dofArray{:}];
            
            for ii = 1:length(obj.dofArray)
                obj.dofArray(ii).setId(ii);
            end
            
            % (3) determine free and fixed dofs
            fixed = obj.dofArray.isFixed();
            obj.fixedDofs = obj.dofArray(fixed);
            obj.freeDofs = obj.dofArray(~fixed);
            
            % (4) set step to 1
            obj.getProperties.addValue('STEP',1);
            
            % (5) build main model part
            obj.femModelParts('MainModelPart') = FemModelPart('MainModelPart', ...
                obj.nodeArray, obj.elementArray, obj);
            
            % (6) initialize the modelparts
            cellfun(@(mp) mp.initialize(), obj.femModelParts.values);
            
            obj.initialized = true;
        end
        
        function node = addNewNode(obj, id, x, y, z)
        %ADDNEWNODE inserts a new node in the fem model with id and x,y,z
        %   coordinates
            if id <= length(obj.nodeArray)
                if obj.nodeArray(id).getId ~= -1
                    msg = ['FemModel: A node with id \"', num2str(id), ...
                        '\" already exists in the model'];
                    e = MException('MATLAB:bm_mfem:duplicateId',msg);
                    throw(e);
                end
            end

            if nargin == 4
                node = Node(id, x, y);
                obj.nodeArray(id) = node;
            elseif nargin == 5
                node = Node(id, x, y, z);
                obj.nodeArray(id) = node;
            else
                error('wrong input parameters')
            end
            
            obj.initialized = false;
        end
        
        function element = addNewElement(obj, elementName, id, nodes, props)
        %ADDNEWELEMENT inserts a new element in the fem model with
        %   elementName, id, and an array of nodes

            if id <= length(obj.elementArray)
                if obj.elementArray(id).getId ~= -1
                    msg = ['FemModel: An element with id \"', num2str(id), ...
                        '\" already exists in the model'];
                    e = MException('MATLAB:bm_mfem:duplicateId',msg);
                    throw(e);
                end
            end
            
            if ~ isa(nodes,'Node')
                nodes = obj.getNodes(nodes);
            end
            
            switch elementName
                case 'BarElement2d2n'
                    element = BarElement2d2n(id, nodes);
                case 'BarElement3d2n'
                    element = BarElement3d2n(id, nodes);
                case 'BeamElement3d2n'
                    element = BeamElement3d2n(id, nodes);
                case 'ConcentratedMassElement3d1n'
                    element = ConcentratedMassElement3d1n(id, nodes);
                case 'SpringDamperElement3d2n'
                    element = SpringDamperElement3d2n(id, nodes);
                case 'ReissnerMindlinElement3d4n'
                    element = ReissnerMindlinElement3d4n(id, nodes);
                case 'ShellElement3d4n'
                    element = ShellElement3d4n(id, nodes);
                case 'DiscreteKirchhoffElement3d4n'
                    element = DiscreteKirchhoffElement3d4n(id, nodes);
                case 'QuadrilateralElement2d4n'
                    element = QuadrilateralElement2d4n(id, nodes);
                case 'QuadrilateralElement2d4nPlaneStrain'
                    element = QuadrilateralElement2d4nPlaneStrain(id, nodes);
                case 'HexahedronElement3d8n'
                    element = HexahedronElement3d8n(id, nodes);
                    
                otherwise
                    error('unknown element %s',elementName)
            end %switch
            
            if nargin == 5
                element.setProperties(props);
            end
            
            obj.elementArray(id) = element;
            obj.initialized = false;
        end
        
    end
    
end

