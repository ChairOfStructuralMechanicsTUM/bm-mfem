classdef GmshInput < ModelIO
    %GMSHINPUT Read input from meshes generated with gmsh
    %   Detailed explanation goes here
    
    properties (Access = private)
        modelPartNames = strings
        props = PropertyContainer
        
        % default element types
        line2n = 'BarElement3d2n'
        triangle3n = 'not yet implemented'
        %etc
    end
    
    methods
        
        %Constructor
        function gmshInput = GmshInput(file)
            if nargin == 0
                super_args = {};
            elseif nargin == 1
                if ~ exist(file, 'file')
                    msg = ['GmshInput: File ', file, ' not found.'];
                    e = MException('MATLAB:bm_mfem:fileNotFound',msg);
                    throw(e);
                end
                super_args = {file};
            end
            gmshInput@ModelIO(super_args{:});
        end
        
        function setElementTypeLine2n(obj, typeName)
            obj.line2n = typeName;
        end
        
        function printElementTypes(obj)
            fprintf('Currently defined element types:\n')
            fprintf('Line2n: %s\n',obj.line2n)
            fprintf('Triangle3n: %s\n',obj.triangle3n)
        end
        
        % member functions
        function model = readModel(obj)
            obj.readProperties;
            modelParts = obj.readModelParts();
            nodes = obj.readNodes();
            elements = obj.readElements(modelParts, nodes);
            model = FemModel(nodes, elements, modelParts);
        end
        
        function nodeArray = readNodes(obj)
            nodeArray = Node.empty;
            fid = fopen(obj.file);
            tline = fgetl(fid);
            
            while ischar(tline)
                if strcmp(tline,'$Nodes')
                    % skip to the first line with nodes
                    tline = fgetl(fid);
                    if (length(str2num(tline)) == 1)
                        tline = fgetl(fid);
                    end
                    
                    % read the data and create the nodes
                    while ~strcmp(tline,'$EndNodes')
                        nodeData = cell2mat(textscan(tline,'%f'))';
                        cNode = Node(nodeData(1),nodeData(2),nodeData(3),nodeData(4));
                        try
                            nodeArray(nodeData(1))
                            idExists = true;
                        catch
                            idExists = false;
                            nodeArray(nodeData(1)) = cNode;
                        end
                        
                        if idExists
                            error('multiple nodes with id %d exist',nodeData(1))
                        end
                        
                        tline = fgetl(fid);
                    end
                    
                end
                tline = fgetl(fid);
            end
            
        end
        
        function elementArray = readElements(obj, modelParts, nodes)
            elementArray = Element.empty;
            fid = fopen(obj.file);
            tline = fgetl(fid);
            
            while ischar(tline)
                if strcmp(tline,'$Elements')
                    % skip to the first line with elements
                    tline = fgetl(fid);
                    if (length(str2num(tline)) == 1)
                        tline = fgetl(fid);
                    end
                    
                    while ~strcmp(tline,'$EndElements')
                        %elm-number elm-type number-of-tags < tag > … node-number-list
                        clear('cElement');
                        elementData = cell2mat(textscan(tline,'%f','TreatAsEmpty',{'x'}));
                        elementType = elementData(2);
                        
                        switch elementType
                            
                            case 1 % 2-node line
                                switch obj.line2n
                                    case 'BarElement2d2n' %id, nodeArray, material, crossSectionArea)
                                        cProperties = obj.props(elementData(4));
                                        
                                        cElement = BarElement2d2n(elementData(1), ...
                                            [nodes(elementData(end-1)) ...
                                            nodes(elementData(end))]);
                                        cElement.setProperties(cProperties);
                                        
                                    case 'BarElement3d2n'
                                        cProperties = obj.props(elementData(4));
                                        
                                        cElement = BarElement3d2n(elementData(1), ...
                                            [nodes(elementData(end-1)) ...
                                            nodes(elementData(end))]);
                                        cElement.setProperties(cProperties);
                                        
                                    case 'BeamElement3d2n'
                                        cProperties = obj.props(elementData(4));
                                        
                                        cElement = BeamElement3d2n(elementData(1), ...
                                            [nodes(elementData(end-1)) ...
                                            nodes(elementData(end))]);
                                        cElement.setProperties(cProperties);
                                        
                                    otherwise
                                        error('2-node element type %s not available',obj.line2n)
                                end
                                
                                
                                
                            case 2 % 3-node triangle
                                cElement = PlaneStressElement2d3n(elementData(1),...
                                    [nodes(elementData(end-2)) nodes(elementData(end-1)) ...
                                    nodes(elementData(end))]);
                                
                            case 3 %4-node quadrangle
                                cElement = PlaneStressElement2d4n(elementData(1),...
                                    [nodes(elementData(end-3)) nodes(elementData(end-2)) ...
                                    nodes(elementData(end-1)) nodes(elementData(end))]);
                                
                            case 9 %6-node second order triangle
                                cElement = PlaneStressElement2d6n(elementData(1),...
                                    [nodes(elementData(end-5)) nodes(elementData(end-4)) ...
                                    nodes(elementData(end-3)) nodes(elementData(end-2)) ...
                                    nodes(elementData(end-1)) nodes(elementData(end))]);
                                
                            case 10 %9-node second order quadrangle
                                cElement = PlaneStressElement2d9n(elementData(1),...
                                    [nodes(elementData(end-8)) nodes(elementData(end-7))...
                                    nodes(elementData(end-6)) nodes(elementData(end-5)) ...
                                    nodes(elementData(end-4)) nodes(elementData(end-3)) ...
                                    nodes(elementData(end-2)) nodes(elementData(end-1)) ...
                                    nodes(elementData(end))]);
                                
                            case 15 %1-node point
                                % here, cElement points to a node; gmsh
                                % distinguishes between nodes and 1-d point
                                % elements, we don't
                                cElement = nodes(elementData(end));
                                
                            case 16 %8-node second order quadrangle
                                cElement = PlaneStressElement2d8n(elementData(1),...
                                    [nodes(elementData(end-7)) nodes(elementData(end-6))...
                                    nodes(elementData(end-5)) nodes(elementData(end-4)) ...
                                    nodes(elementData(end-3)) nodes(elementData(end-2)) ...
                                    nodes(elementData(end-1)) nodes(elementData(end))]);
                                
                        end %switch
                        
                        
                        
                        %add the element to the array and the modelParts map
                        if exist('cElement','var') == 1
                            nModelPartId = elementData(4);
                            modelPartName = char(obj.modelPartNames(nModelPartId));
                            nModelPart = modelParts(modelPartName);
                            
                            if elementType == 15 % don't add nodes to the array
                                nModelPart.addNode(cElement);
                                modelParts(modelPartName) = nModelPart;
                            else
                                elementArray = [elementArray cElement];
                                nModelPart.addElement(cElement);
                                modelParts(modelPartName) = nModelPart;
                            end
                            
                        end
                        
                        tline = fgetl(fid);
                    end %while element
                end
                
                tline = fgetl(fid);
            end %while
        end
        
        function modelParts = readModelParts(obj)
            % the "physical tag" will be used to assign material properties
            % and boundary conditions. 2d: points -> bc, lines -> material
            % 3d: lines -> bc, tetras.. -> material
            modelParts = containers.Map;
            fid = fopen(obj.file);
            tline = fgetl(fid);
            
            while ischar(tline)
                if strcmp(tline,'$PhysicalNames')
                    % generates an array with all model part names, which
                    % is stored as member variable; the place in the
                    % modelPartNames array is determined by its id; the
                    % returned modelParts is an empty map with all names,
                    % which is to be populated in the readElements function
                    
                    tline = fgetl(fid);
                    if (length(str2num(tline)) == 1)
                        tline = fgetl(fid);
                    end
                    
                    while ~strcmp(tline,'$EndPhysicalNames')
                        modelPartNameString = strsplit(tline);
                        nModelPartName = str2double(modelPartNameString(end-1));
                        modelPartName = char(strrep(modelPartNameString(end),'"',''));
                        
                        obj.modelPartNames(nModelPartName) = modelPartName;
                        modelParts(modelPartName) = FemModelPart(modelPartName,[],[]);
                        %                         modelParts(modelPartName) = struct('nodes',[],...
                        %                             'elements',[]);
                        tline = fgetl(fid);
                    end
                    
                end
                
                tline = fgetl(fid);
            end
        end
        
        function readProperties(obj)
            fid = fopen(obj.file);
            tline = fgetl(fid);
            
            while ischar(tline)
                if regexp(tline,'Properties') == 2
                    % generates the properties as a map with key: variable name, value:
                    % variable value; the place in the propertyData array is determined
                    % by the property id
                    nProperty = str2double(strsplit(tline));
                    nProperty = nProperty(2);
                    property = PropertyContainer;
                    
                    tline = fgetl(fid);
                    
                    while ~strcmp(tline,'$EndProperties')
                        prop = strsplit(tline);
                        property.addValue(cell2mat(prop(1)), str2double(prop(2)));
                        tline = fgetl(fid);
                    end
                    
                    obj.props(nProperty) = property;
                    
                end
                tline = fgetl(fid);
            end
        end
        
        
    end
    
end

