classdef ModelIO < handle
    %MODELIO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        inputFile
        modelPartNames = strings
        props = PropertyContainer
        
        % default element types
        line2n = 'BarElement3d2n'
        triangle3n = 'not yet implemented'
        %etc
    end
    
    methods
        
        %Constructor
        function modelIO = ModelIO(file)
            modelIO.inputFile = file;
        end
        
        function setElementTypeLine2n(modelIO, typeName)
           modelIO.line2n = typeName; 
        end
        
        function printElementTypes(modelIO)
           fprintf('Currently defined element types:\n')
           fprintf('Line2n: %s\n',modelIO.line2n)
           fprintf('Triangle3n: %s\n',modelIO.triangle3n)
        end
        
        % member functions
        function model = readModel(modelIO)
            modelIO.readProperties;
            modelParts = modelIO.readModelParts();
            nodes = modelIO.readNodes();
            elements = modelIO.readElements(modelParts, nodes);
            model = FemModel(nodes, elements, modelParts);
        end
        
        function nodeArray = readNodes(modelIO)
            nodeArray = Node.empty;
            fid = fopen(modelIO.inputFile);
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
        
        function elementArray = readElements(modelIO, modelParts, nodes)
            elementArray = Element.empty;
            fid = fopen(modelIO.inputFile);
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
                                switch modelIO.line2n
                                    case 'BarElement2d2n' %id, nodeArray, material, crossSectionArea)
                                        cProperties = modelIO.props(elementData(4));
                                        
                                        cElement = BarElement2d2n(elementData(1), ...
                                            [nodes(elementData(end-1)) ...
                                            nodes(elementData(end))]);
                                        cElement.setProperties(cProperties);
                                        
                                    case 'BarElement3d2n'
                                        cProperties = modelIO.props(elementData(4));
                                        
                                        cElement = BarElement3d2n(elementData(1), ...
                                            [nodes(elementData(end-1)) ...
                                            nodes(elementData(end))]);
                                        cElement.setProperties(cProperties);
                                        
                                    otherwise
                                        error('2-node element type %s not available',modelIO.line2n)
                                end
                                
                                
                                
                            case 2 % 3-node triangle
                                
                            case 3 %4-node quadrangle
                                
                            case 15 %1-node point
                                % here, cElement points to a node; gmsh
                                % distinguishes between nodes and 1-d point
                                % elements, we don't
                                cElement = nodes(elementData(end));
                                
                        end %switch
                        
                        
                        
                        %add the element to the array and the modelParts map
                        if exist('cElement','var') == 1
                            % don't add nodes to the array
                            if elementType ~= 15
                                elementArray = [elementArray cElement];
                            end
                            nModelPart = elementData(4);
                            
                            if modelParts(char(modelIO.modelPartNames(nModelPart))) == 0
                                modelParts(char(modelIO.modelPartNames(nModelPart))) = cElement;
                            else
                                modelParts(char(modelIO.modelPartNames(nModelPart))) = ...
                                    [modelParts(char(modelIO.modelPartNames(nModelPart))) cElement];
                            end
                            
                        end
                        
                        tline = fgetl(fid);
                    end %while element
                end
                
                tline = fgetl(fid);
            end %while
        end
        
        function modelParts = readModelParts(modelIO)
            % the "physical tag" will be used to assign material properties
            % and boundary conditions. 2d: points -> bc, lines -> material
            % 3d: lines -> bc, tetras.. -> material
            modelParts = containers.Map;
            fid = fopen(modelIO.inputFile);
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
                        
                        modelIO.modelPartNames(nModelPartName) = modelPartName;
                        modelParts(modelPartName) = 0;
                        tline = fgetl(fid);
                    end
                    
                end
                
                tline = fgetl(fid);
            end
        end
        
        function readProperties(modelIO)
            fid = fopen(modelIO.inputFile);
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
                        property.setValue(cell2mat(prop(1)), str2double(prop(2)));
                        tline = fgetl(fid);
                    end
                    
                    modelIO.props(nProperty) = property;
                    
                end
                tline = fgetl(fid);
            end
        end
        
        
    end
    
end

