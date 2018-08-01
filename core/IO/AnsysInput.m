classdef AnsysInput < ModelIO
    %AnsysInput Read model input from mdpa
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        
        function ansysInput = AnsysInput(file)
            % file - full file path with filetype as extensions
            % examample: file='F:\ALMA_simulations\AnsysModels\plate.inp';
            % ansysExecutable - full file path of Ansys executable
            %
            if nargin == 0
                super_args = {};
            elseif nargin == 1
                if ~ exist(file, 'file')
                    msg = ['AnsysInput: File ', file, ' not found.'];
                    e = MException('MATLAB:bm_mfem:fileNotFound',msg);
                    throw(e);
                end
                super_args = {file};
            end
            ansysInput@ModelIO(super_args{:});
            
        end
        
        function model = readModel(obj,ansysExecutable)
                        
            A=strsplit(obj.file,'\');
            
            B=strsplit(A{end},'.');
            
            C=strsplit(obj.file,B{1});
            
            if strcmp(A{1},obj.file)
                folder = [pwd '\'];
            else
                folder=C{1};
            end
            
            if strcmp(B{1},A{end})
                msg = 'AnsysInput: file type not found.';
                e = MException('MATLAB:bm_mfem:filetypeNotFound',msg);
                throw(e);
            else
                fileName=B{1};
                extension=B{end};
            end
            
            disp('Ansys input file selected')
            disp(['folder: ' folder])
            disp(['filename: ' fileName])
            disp(['filetype: ' extension])
            
            data=obj.runAnsys(ansysExecutable,folder,fileName,extension);
            
            props = PropertyContainer;
            model = FemModel;

            
%% Continue here
%             tline = fgetl(fid);
%             
%             while ischar(tline)
%                 if startsWith(strtrim(tline),'//'); tline = fgetl(fid); continue; end
%                 
%                 if contains(tline,'Begin Properties')
%                     nProp = strsplit(tline);
%                     nProp = str2double(nProp{end});
%                     [fid,props] = obj.readProperties(fid,props,nProp);
%                 end
%                 
%                 if contains(tline,'Begin Nodes')
%                     [fid,model] = obj.readNodes(fid,model);
%                 end
%                 
%                 if contains(tline,'Begin Elements')
%                     etype = strsplit(tline);
%                     etype = cell2mat(etype(3));
%                     etype = erase(etype,'//');
%                     [fid,model] = obj.readElements(fid,model,etype,props);
%                 end
%                 
%                 if contains(tline,'Begin SubModelPart')
%                     name = strsplit(tline);
%                     name = cell2mat(name(3));
%                     [fid,model] = obj.readSubModelParts(fid,model,name);
%                 end
%                 
%                 tline = fgetl(fid);
            end
            
        end
        
    end
    
    methods (Static)
        
        function data=runAnsys(ansysExecutable,folder,file,extension)
            if(nargin > 0)
                % make directory for files
                mkdir('DataAnsys')
                % Header
                fidl=fopen('DataAnsys/modelFile.txt','w');
                % Read input file
                fprintf(fidl,'/INPUT,%s,%s,%s,,0 \r\n',file...
                    ,extension...
                    ,folder);
                % Enters the preprocessor
                fprintf(fidl,'/PREP7 \r\n');
                % Set EMAWRITE to "yes" to obtain binary files
                fprintf(fidl,'EMATWRITE,YES \r\n');
                % Create external file with the nodal information
                fprintf(fidl,'/output,DataAnsys/nodeCoor,txt \r\n');
                fprintf(fidl,'NLIST,ALL,,,XYZ,NODE,NODE,NODE \r\n');
                fprintf(fidl,'/output \r\n');
                % Create external file with the restricted DoF
                fprintf(fidl,'/output,DataAnsys/nodeRest,txt \r\n');
                fprintf(fidl,'DLIST,ALL,\r\n');
                fprintf(fidl,'/output \r\n');
                % Create external file with the element type
                fprintf(fidl,'/output,DataAnsys/elemTyp,txt \r\n');
                fprintf(fidl,'ETLIST,ALL,\r\n');
                fprintf(fidl,'/output \r\n');
                % Create external file with the element type
                fprintf(fidl,'/output,DataAnsys/elemNodes,txt \r\n');
                fprintf(fidl,'ELIST,ALL,\r\n');
                fprintf(fidl,'/output \r\n');
                % Create external file with the material information
                fprintf(fidl,'/output,DataAnsys/materialPro,txt \r\n');
                fprintf(fidl,'MPLIST,ALL,,,EVLT\r\n');
                fprintf(fidl,'/output \r\n');
                % Set solver
                fprintf(fidl,'/SOLU\r\n');
                fprintf(fidl,'ANTYPE,MODAL\r\n');
                fprintf(fidl,'MODOPT,DAMP,1\r\n');
                fprintf(fidl,'MXPAND,1,,,YES\r\n');
                fprintf(fidl,'WRFULL,1\r\n');
                fprintf(fidl,'SOLVE\r\n');
                fprintf(fidl,'FINISH\r\n');
                % Extract the mass, damping and stifness matrix
                fprintf(fidl,'/AUX2\r\n');
                fprintf(fidl,'file,,full\r\n');
                
                fprintf(fidl,...
                    'HBMAT,DataAnsys/HBMstiff,txt,,ascii,stiff,yes,yes\r\n');
                
                fprintf(fidl,...
                    'HBMAT,DataAnsys/HBMmass,txt,,ascii,mass,yes,yes\r\n');
                
                fprintf(fidl,...
                    'HBMAT,DataAnsys/HBMdamp,txt,,ascii,damp,yes,yes\r\n');
                
                fprintf(fidl,'FINISH\r\n');
                % Get the modified nodal information
                fprintf(fidl,'/AUX2\r\n');
                fprintf(fidl,'/output,DataAnsys/record_5,txt \r\n');
                fprintf(fidl,'form,long\r\n');
                fprintf(fidl,'fileaux2,,emat,%s\r\n',folder);
                fprintf(fidl,'dump,5,5\r\n');
                fprintf(fidl,'/output \r\n');
                fprintf(fidl,'FINISH\r\n');
                fclose(fidl);
                % Run ANSYS
                % Warning: modify path according to the user
                %!"C:\Program Files\ANSYS Inc\v171\ansys\bin\winx64\ANSYS171.exe" -b  -i DataAnsys/modelFile.txt -o DataAnsys/result.out
                eval(['!"' ansysExecutable '" -b  -i DataAnsys/modelFile.txt -o DataAnsys/result.out'])
                
                
                
                % See http://people.sc.fsu.edu/~jburkardt/m_src/hb_to_msm/hb_to_msm.html
                % The function hb_to_msm is freely distributed according to the website
                data.Mansys = hb_to_msm('DataAnsys/HBMmass.txt');
                
                %                data.M = full(Mansys)+transpose(full(Mansys))...
                %                    -diag(diag(full(Mansys)));
                
                data.Kansys = hb_to_msm('DataAnsys/HBMstiff.txt');
                %                data.K = full(Kansys)+transpose(full(Kansys))...
                %                    -diag(diag(full(Kansys)));
                
                data.Cansys = hb_to_msm('DataAnsys/HBMdamp.txt');
                %                data.C = full(Cansys)+transpose(full(Cansys))...
                %                    -diag(diag(full(Cansys)));
                %                if  sum(sum(abs(data.C))) ~= 0
                %                    data.damped = true;
                %                end
                % Delete files
                try
                    delete('*.emat');
                    delete('*.esav');
                    delete('*.full');
                    delete('*.mlv');
                    delete('*.err');
                    delete('*.db');
                    delete('*.log');
                catch
                end
                try
                    delete('*.BCS');
                    delete('*.ce');
                    delete('*.mode');
                    delete('*.stat');
                    delete('*.xml');
                catch
                end
                
                % Read restriction on the nodes
                data.nodeRest = readRestrictions();
                % Read record 5 of ANSYS to get the position of the entries
                % of each node with the matrices
                nodeEquiv = readRecord_5();
                % Read the coordinates of each node
                data.nodeList = readCoord();
                % Order the coordinates acoording to the record 5
                nodesC = data.nodeList(nodeEquiv,:);
                data.nodesOrderByDofs=nodesC(:,1)';
                
                % Create objects "Node" and assign them to a model
                arrayAux(1,size(nodesC,1)) = Node;
                for i = 1 : size(nodesC,1)
                    referenceNode =  nodesC(i,1);
                    x  =  nodesC(i,2);
                    y  =  nodesC(i,3);
                    z  =  nodesC(i,4);
                    arrayAux(1,i) = Node(referenceNode,x,y,z);
                end
                data.nodes = arrayAux;
                data.numNodes = length(arrayAux);
                
                % Create available elements
                %elementTypes = createElements(data.dimension);
                [data.elementsOfModel,data.nodeElementList,data.nodeConnectivity] = readElements();
                
                rmdir('DataAnsys', 's')
            end
        end
        
        
        
    end
    
end

