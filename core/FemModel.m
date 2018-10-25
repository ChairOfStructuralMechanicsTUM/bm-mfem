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
        function mp = addNewModelPart(obj, name, nodes, elements)
        %ADDNEWMODELPART add a new model part
        %   ADDNEWMODELPART(name, nodes, elements) adds a modelpart to the
        %   FemModel containing nodes and/or elements.
        %
        %   ADDNEWMODELPART(name, nodeIds, elementIds) adds a modelpart to
        %   the FemModel containing nodes with nodeIds and/or elements with
        %   elementIds.
        %
        %   ANNEWMODELPART(name) adds an empty model part
        %
        %   see also FEMMODELPART
            if nargin == 2
                obj.femModelParts(name) = FemModelPart(name, [], [], obj);
            elseif nargin == 4
                if ~isa(nodes,'Node')
                    nodes = obj.getNodes(nodes);
                end
                if ~isa(elements,'Element')
                    elements = obj.getElements(elements);
                end
                
                obj.femModelParts(name) = FemModelPart(name, ...
                    nodes, elements, obj);
            else
                msg = 'FemModel: Wrong number of input arguments';
                e = MException('MATLAB:bm_mfem:invalidArguments',msg);
                throw(e);
            end
            
            mp = obj.femModelParts(name);
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
        
        function addNodes(obj, nodes)
            %ADDNODES adds an array of nodes to the model
            nodeIds = nodes.getId();
            if ~isempty(obj.nodeArray)
                nodeIds = [nodeIds obj.nodeArray.getId()];
            end
            duplicated_nodes = findDuplicates(nodeIds);
            if ~ isempty(duplicated_nodes)
                error('multiple nodes with id %d exist',duplicated_nodes);
            else
                obj.nodeArray(nodeIds) = nodes;
            end
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
            
            % this creates elements based on their string name
            try element = feval(elementName, id, nodes);
                if ~isa(element,'Element')
                    msg = [class(obj), ': Only elements can be created.'];
                    e = MException('MATLAB:bm_mfem:invalidInput',msg);
                    throw(e);
                end
            catch
                msg = [class(obj), ': Unknown element \"', elementName, '\".'];
                e = MException('MATLAB:bm_mfem:unknownElement',msg);
                throw(e);
            end
            
            if nargin == 5
                element.setProperties(props);
            end
            
            obj.elementArray(id) = element;
            obj.initialized = false;
        end
        
    end
    
end
