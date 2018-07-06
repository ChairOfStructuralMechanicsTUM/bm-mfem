classdef MdpaInput < ModelIO
    %MDPAINPUT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Access = private)
        props = PropertyContainer
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
            tline = fgetl(fid);
            
            while ischar(tline)
                if contains(tline,'Begin Properties')
                    nProp = strsplit(tline);
                    nProp = str2double(nProp{end});
                    fid = obj.readProperties(fid, nProp);
                end
                
                tline = fgetl(fid);
            end
%             obj.readProperties;
%             modelParts = obj.readModelParts();
%             nodes = obj.readNodes();
%             elements = obj.readElements(modelParts, nodes);
%             model = FemModel(nodes, elements, modelParts);
        end
        
        function fid = readProperties(obj, fid, nProp)
            tline = fgetl(fid);
            property = PropertyContainer;
            while ~ strcmp(tline, 'End Properties')
                prop = strsplit(tline);
                property.addValue(cell2mat(prop(1)), str2double(prop(2)));
                tline = fgetl(fid);
            end
            obj.props(nProp) = property;
        end
        
    end
    
end

