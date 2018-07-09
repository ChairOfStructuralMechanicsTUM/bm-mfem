classdef MdpaInput < ModelIO
    %MDPAINPUT Read model input from mdpa
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        
        function mdpaInput = MdpaInput(file)
            if nargin == 0
                super_args = {};
            elseif nargin == 1
                if ~ exist(file, 'file')
                    msg = ['MdpaInput: File ', file, ' not found.'];
                    e = MException('MATLAB:bm_mfem:fileNotFound',msg);
                    throw(e);
                end
                super_args = {file};
            end
            mdpaInput@ModelIO(super_args{:});
        end
        
        function model = readModel(obj)
            fid = fopen(obj.file);
            props = PropertyContainer;
            model = FemModel;
            
            tline = fgetl(fid);
            
            while ischar(tline)
                if startsWith(strtrim(tline),'//'); tline = fgetl(fid); continue; end
            
                if contains(tline,'Begin Properties')
                    nProp = strsplit(tline);
                    nProp = str2double(nProp{end});
                    [fid,props] = obj.readProperties(fid,props,nProp);
                end
                
                if contains(tline,'Begin Nodes')
                    [fid,model] = obj.readNodes(fid,model);
                end
                
                if contains(tline,'Begin Elements')
                    etype = strsplit(tline);
                    etype = cell2mat(etype(3));
                    etype = erase(etype,'//');
                    [fid,model] = obj.readElements(fid,model,etype,props);
                end
                
                if contains(tline,'Begin SubModelPart')
                    name = strsplit(tline);
                    name = cell2mat(name(3));
                    [fid,model] = obj.readSubModelParts(fid,model,name);
                end
                
                tline = fgetl(fid);
            end
            
        end
        
    end
    
    methods (Static)
        
        function [fid,props] = readProperties(fid,props,nProp)
            tline = fgetl(fid);
            property = PropertyContainer;
            while ~ strcmp(tline,'End Properties')
                if startsWith(strtrim(tline),'//'); tline = fgetl(fid); continue; end
                
                propData = strsplit(tline);
                propData(cellfun('isempty',propData)) = [];
                property.addValue(cell2mat(propData(1)),str2double(propData(2)));
                tline = fgetl(fid);
            end
            props(nProp) = property;
        end
        
        function [fid,model] = readNodes(fid,model)
            tline = fgetl(fid);
            while ~ strcmp(tline,'End Nodes')
                if startsWith(strtrim(tline),'//'); tline = fgetl(fid); continue; end
                
                n = cell2mat(textscan(tline,'%f'))';
                model.addNewNode(n(1),n(2),n(3),n(4));
                tline = fgetl(fid);
            end
        end
        
        function [fid,model] = readElements(fid,model,etype,props)
            tline = fgetl(fid);
            while ~ strcmp(tline,'End Elements')
                if startsWith(strtrim(tline),'//'); tline = fgetl(fid); continue; end
                
                e = cell2mat(textscan(tline,'%f'));
                model.addNewElement(etype,e(1),e(3:end),props(e(2)));
                tline = fgetl(fid);
            end
        end
        
        function [fid,model] = readSubModelParts(fid,model,name)
            tline = fgetl(fid);
            nodes = [];
            elements = [];
            while ~ strcmp(tline,'End SubModelPart')
                if startsWith(strtrim(tline),'//'); tline = fgetl(fid); continue; end
                
                if contains(tline,'Begin SubModelPartNodes')
                    tline = fgetl(fid);
                    while ~ contains(tline,'End SubModelPartNodes')
                        nodes = [nodes str2double(tline)]; %#ok<AGROW>
                        tline = fgetl(fid);
                    end
                end
                
                if contains(tline,'Begin SubModelPartElements')
                    tline = fgetl(fid);
                    while ~ contains(tline,'End SubModelPartElements')
                        elements = [elements str2double(tline)]; %#ok<AGROW>
                        tline = fgetl(fid);
                    end
                end
                
                tline = fgetl(fid);
            end
            
            model.addNewModelPart(name,nodes,elements);
        end
        
    end
    
end

