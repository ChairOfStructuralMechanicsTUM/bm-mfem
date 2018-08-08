classdef AnsysInput < ModelIO
    %AnsysInput Read model input from mdpa
    %   Detailed explanation goes here
    
    properties (Access = private)
        ansysExecutable
    end
    
    methods
        
        function obj = AnsysInput(file, ansysExecutable)
            % file - full file path with filetype as extensions
            % examample: file='F:\ALMA_simulations\AnsysModels\plate.inp';
            % ansysExecutable - full file path of Ansys executable
            %
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~ exist(file, 'file')
                    msg = ['AnsysInput: File ', file, ' not found.'];
                    e = MException('MATLAB:bm_mfem:fileNotFound',msg);
                    throw(e);
                end
                super_args = {file};
            else
                msg = 'AnsysInput: Wrong number of input parameters.';
                e = MException('MATLAB:bm_mfem:invalidInput',msg);
                throw(e);
            end
            
            obj@ModelIO(super_args{:});
            
            obj.ansysExecutable = ansysExecutable;
        end
        
        function model = readModel(obj)
            
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
            
            model = FemModel;
            mp = model.addNewModelPart('ANSYS_import');
            data = obj.runAnsys(obj.ansysExecutable,folder,fileName,extension);
            
            % Read restriction on the nodes
            data.nodeRest = obj.readRestrictions();
            % Read record 5 of ANSYS to get the position of the entries
            % of each node with the matrices
            nodeEquiv = obj.readRecord_5();
            % Read the coordinates of each node
            data.nodeList = obj.readCoord();
            % Order the coordinates acoording to the record 5
            nodesC = data.nodeList(nodeEquiv,:);
            data.nodesOrderByDofs=nodesC(:,1)';
            
            % Create objects "Node" and assign them to a model
            for ii = 1:size(nodesC,1)
                id = nodesC(ii,1);
                x = nodesC(ii,2);
                y = nodesC(ii,3);
                z = nodesC(ii,4);
                mp.addNewNode(id,x,y,z);
            end
            data.numNodes = length(mp.getNodes());
            
            % Read available element data
            [data.elementsOfModel,data.nodeElementList,data.nodeConnectivity] = AnsysInput.readElements();
            
            % Assign dofs to nodes
            for ii = 1 :size(data.elementsOfModel,1)
                nodeData=cell2mat(data.nodeElementList(ii));
                dofs = obj.getDofsOfAnsysElements(data.elementsOfModel{ii,1}, data.elementsOfModel{ii,2});
                nodeData(isnan(nodeData)) = [];
                for j = 1:length(nodeData)
                    n = mp.getNodeById(nodeData(j));
                    n.addDof(dofs);
                end
            end
            % Set dof ids according to node numbering (bm-mfem style,
            % reordering is done in the DummyElement
            dofArray = arrayfun(@(node) node.getDofArray, mp.getNodes(), 'UniformOutput', false);
            dofArray = [dofArray{:}];
            for ii = 1:length(dofArray)
                dofArray(ii).setId(ii);
            end
            
            % Create the dummy element
            e = mp.addNewElement('DummyElement',1,mp.getNodes);
            
            % Set the matrices
            systemSize = size(data.Mansys,1);
            Mdiag = spdiags(data.Mansys,0);
            M = data.Mansys + data.Mansys.' - spdiags(Mdiag(:),0,systemSize,systemSize);
            Cdiag = spdiags(data.Cansys,0);
            C = data.Cansys + data.Cansys.' - spdiags(Cdiag(:),0,systemSize,systemSize);
            Kdiag = spdiags(data.Kansys,0);
            K = data.Kansys + data.Kansys.' - spdiags(Kdiag(:),0,systemSize,systemSize);
            e.setMatrices(M, C, K);
            
            % Set the dof ordering in the dummy element
            e.setDofOrder(data.nodesOrderByDofs);
            
            % Set restrictions
            obj.setRestrictions(data.nodeRest, mp);
            
            % Set node connectivity
            e.setNodeConnectivity(data.nodeConnectivity);
            
            % Delete auxiliary files
            
%             try
%                 rmdir('DataAnsys', 's')
%             catch
%             end
            
        end
        
    end
    
    methods (Static)
        
        function data = runAnsys(ansysExecutable,folder,file,extension)
            % make directory for files
            if exist('DataAnsys','dir') ~= 7; mkdir('DataAnsys'); end
            
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
            eval(['!"' ansysExecutable '" -b  -i DataAnsys/modelFile.txt -o DataAnsys/result.out'])            
            
            % Import sparse matrices form harwell boeing format to matlab
            % See http://people.sc.fsu.edu/~jburkardt/m_src/hb_to_msm/hb_to_msm.html
            % The function hb_to_msm is freely distributed according to the website
            data.Mansys = hb_to_msm('DataAnsys/HBMmass.txt');
            data.Kansys = hb_to_msm('DataAnsys/HBMstiff.txt');
            data.Cansys = hb_to_msm('DataAnsys/HBMdamp.txt');
            
            % Delete files
            try
                delete('*.emat');
                delete('*.esav');
                delete('*.full');
                delete('*.mlv');
                delete('*.err');
                delete('*.db');
                delete('*.log');
                delete('*.BCS');
                delete('*.ce');
                delete('*.mode');
                delete('*.stat');
                delete('*.xml');
            catch
            end
            
        end
        
        function A = readCoord()
            % function A = readCoord()
            %
            % Function   : readCoord
            %
            % Description: This function gets the nodes coordinates from txt files
            %              created in ANSYS
            %
            % Parameters :
            %
            % Return     : A                   - matrix with nodal information
            %
            fid=fopen('DataAnsys/nodeCoor.txt') ;
            fidd=fopen('DataAnsys/nodeCoor_modified.dat','w') ;
            if fid < 0, error('Cannot open file'); end
            % Discard some line to read the data from the txt files
            for j = 1 : 13
                fgetl(fid) ;
            end
            
            while ~feof(fid)
                tline=fgets(fid);
                if isspace(tline)
                    for j = 1 : 9
                        fgetl(fid) ;
                    end
                else
                    fwrite(fidd,tline) ;
                end
            end
            
            fclose all ;
            filename = 'DataAnsys/nodeCoor_modified.dat';
            delimiterIn = ' ';
            % Get data in matlab
            A = importdata(filename,delimiterIn);
        end
        
        function [elementList,nodesArrays,nodeConnectivity] = readElements()
            % function [elementList,nodesArrays] = readElements()
            %
            % Function   : readElements
            %
            % Description: This function gets the element information from the txt
            %              files generated by ANSYS
            %
            % Parameters :
            %
            % Return     : elementList  - cell array with element names in
            %                             elementList{:,1} and keyopts in
            %                             elementList{:,2}
            %              nodesArrays  - array with the nodes related
            %                             to certain element
            fid=fopen('DataAnsys/elemNodes.txt') ;
            fidd=fopen('DataAnsys/elemNodes_modified.dat','w') ;
            if fid < 0, error('Cannot open file'); end
            % Discard some line to read the data from the txt files
            for j = 1 : 13
                fgetl(fid) ;
            end
            while ~feof(fid)
                tline=fgets(fid);
                if isspace(tline)
                    for j = 1 : 10
                        fgetl(fid) ;
                    end
                else
                    fwrite(fidd,tline) ;
                end
            end
            fclose all ;
            filename = 'DataAnsys/elemNodes_modified.dat';
            delimiterIn = ' ';
            % Get data in matlab
            A = importdata(filename,delimiterIn);
            
            nodeConnectivity=A(:,7:end);
            
            A(:,1) = [];
            A(:,1) = [];
            A(:,2) = [];
            A(:,2) = [];
            A(:,2) = [];
            % Check how many element types there are
            elements = unique(A(:,1));
            %nodesArrays = cell(length(elements),1);
            aux1 = [];
            aux2 = [];
            nodesArrays = cell(length(elements),1);
            for i = 1 : length(elements)
                for j = 1 : size(A,1)
                    if A(j,1) == elements(i)
                        for k = 2:size(A,2)
                            aux1 = [aux1 A(j,k)];
                        end
                        aux2 = [aux2 aux1];
                    end
                end
                nodesArrays(i) = {unique(aux2)};
                aux1 = [];
                aux2 = [];
            end
            
            % read element types with their keyopts
            fid=fopen('DataAnsys/elemTyp.txt');
            fgetl(fid);
            tline = fgetl(fid);
            tmp = strsplit(strtrim(tline),' ');
            elementList = cell(str2double(tmp{7}),2);    %array for element type and keyopts
            tline = fgetl(fid);
            
            while ~ feof(fid)
                if contains(tline,'ELEMENT TYPE')
                    tmp = strsplit(strtrim(tline),' ');
                    n_etype = str2double(tmp{3});
                    elementList{n_etype,1} = tmp{5};
                    keyopts = zeros(1,18);
                    for ii = 0:2
                        tline = fgetl(fid);
                        tmp = str2double(strsplit(strtrim(tline),' '));
                        keyopts(ii*6+1:ii*6+6) = tmp(end-5:end);
                    end
                    elementList{n_etype,2} = keyopts;
                end
                tline = fgetl(fid);
            end
            
            
        end
        
        function nodeNum = readRecord_5()
            % function nodeNum = readRecord_5()
            %
            %
            % Function   : readRecord_5
            %
            % Description: This function get the nodal equivalence between the
            %              original number of node given in ANSYS and the distribution
            %              within the stiffness matrix through reading record5
            %
            % Parameters :
            %
            % Return     : nodeNum                   - array
            %
            fid = fopen('DataAnsys/record_5.txt', 'r') ;
            if fid < 0, error('Cannot open file'); end
            for i = 1 : 10
                fgetl(fid) ;
            end
            buffer = fread(fid, Inf) ;
            fclose(fid);
            fid = fopen('DataAnsys/record_5_modified.txt', 'w')  ;
            fwrite(fid, buffer) ;
            fclose(fid) ;
            filename = 'DataAnsys/record_5_modified.txt';
            delimiterIn = ' ';
            A = importdata(filename,delimiterIn);
            %wenn nur ein Element in FEmodel dann ist A kein struct
            if ~isstruct(A)
                B.textdata=strtrim(cellstr(num2str(A'))');
                A=B;
            else
                A.textdata(:,6)=[];
            end
            nodeNum = [];
            for i = 1 : size(A.textdata,1)
                %for j = 1 : 5
                for j = 1 : size(A.textdata,2)
                    if isnan(str2double(A.textdata(i,j)))
                        break;
                    end
                    nodeNum = [nodeNum str2double(A.textdata(i,j))];
                end
            end
        end
        
        function A = readRestrictions()
            % function A = readRestrictions()
            %
            % Function   : readRestrictions
            %
            % Description: This function get the restrictions on the coordinates from
            %              ANSYS
            %
            % Parameters :
            %
            % Return     : A cell with restrictions. 
            %               A{1}: node numbers
            %               A{2}: resticted dof name
            %               A{3}: real values the dof is restricted to
            %               A{4}: imaginary values the dof is restricted to
            
            fid=fopen('DataAnsys/nodeRest.txt');
            for j = 1:13
                fgetl(fid); % Read/discard line.
            end
            
            A = textscan(fid,'%u%s%f%f');
            
        end
        
        function setRestrictions(restrictions, mp)
           for ii = 1:length(restrictions{1})
                n = mp.getNodeById(restrictions{1}(ii));
                if restrictions{3}(ii) == 0
                    switch restrictions{2}{ii}
                        case 'UX'
                            mp.getElementById(1).setDofRestriction(n, 'DISPLACEMENT_X');
                        case 'UY'
                            mp.getElementById(1).setDofRestriction(n, 'DISPLACEMENT_Y');
                        case 'UZ'
                            mp.getElementById(1).setDofRestriction(n, 'DISPLACEMENT_Z');
                        case 'ROTX'
                            mp.getElementById(1).setDofRestriction(n, 'ROTATION_X');
                        case 'ROTY'
                            mp.getElementById(1).setDofRestriction(n, 'ROTATION_Y');
                        case 'ROTZ'
                            mp.getElementById(1).setDofRestriction(n, 'ROTATION_Z');
                        otherwise
                            error('required restriction not yet implemented')
                    end
                else
                    switch restrictions{2}{ii}
                        case 'UX'
                            n.setDofValue('DISPLACEMENT_X', restrictions{3}(ii));
                        case 'UY'
                            n.setDofValue('DISPLACEMENT_Y', restrictions{3}(ii));
                        case 'UZ'
                            n.setDofValue('DISPLACEMENT_Z', restrictions{3}(ii));
                        case 'ROTX'
                            n.setDofValue('ROTATION_X', restrictions{3}(ii));
                        case 'ROTY'
                            n.setDofValue('ROTATION_Y', restrictions{3}(ii));
                        case 'ROTZ'
                            n.setDofValue('ROTATION_Z', restrictions{3}(ii));
                        otherwise
                            error('required restriction not yet implemented')
                    end
                end
            end 
        end
        
        function dofs = getDofsOfAnsysElements(elementName, keyopts)
            ux = "DISPLACEMENT_X";
            uy = "DISPLACEMENT_Y";
            uz = "DISPLACEMENT_Z";
            rx = "ROTATION_X";
            ry = "ROTATION_Y";
            rz = "ROTATION_Z";
            
            switch 1
                case strcmp(elementName,'BEAM3')
                    dofs = [ux uy rz];
                    
                case strcmp(elementName,'COMBIN14')
                    if keyopts(2) == 0
                        if keyopts(3) == 0
                            dofs = [ux uy uz];
                        elseif keyopts(3) == 1
                            dofs = [rx ry rz];
                        elseif keyopts(3) == 2
                            dofs = [ux uy];
                        elseif keyopts(3) == 4
                            dofs = [ux uy];
                        else
                            msg = ['AnsysInput: Invalid keyopts for element type ', ...
                                elementName];
                            e = MException('MATLAB:bm_mfem:invalidKeyopts',msg);
                            throw(e);
                        end
                    elseif keyopts(2) == 1
                        dofs = ux;
                    elseif keyopts(2) == 2
                        dofs = uy;
                    elseif keyopts(2) == 3
                        dofs = uz;
                    elseif keyopts(2) == 4
                        dofs = rx;
                    elseif keyopts(2) == 5
                        dofs = ry;
                    elseif keyopts(2) == 6
                        dofs = rz;
                    else
                        msg = ['AnsysInput: Invalid keyopts for element type ', ...
                            elementName];
                        e = MException('MATLAB:bm_mfem:invalidKeyopts',msg);
                        throw(e);
                    end
                    
                case strcmp(elementName,'MASS21')
                    if keyopts(3) == 0
                        dofs = [ux uy uz rx ry rz];
                    elseif keyopts(3) == 2
                        dofs = [ux uy uz];
                    elseif keyopts(3) == 3
                        dofs = [ux uy rz];
                    elseif keyopts(3) == 4
                        dofs = [ux uy];
                    else
                        msg = ['AnsysInput: Invalid keyopts for element type ', ...
                            elementName];
                        e = MException('MATLAB:bm_mfem:invalidKeyopts',msg);
                        throw(e);
                    end
                    
                case strcmp(elementName,'SHELL63')
                    dofs = [ux uy uz rx ry rz];
                    
                case strcmp(elementName,'SHELL181')
                    if keyopts(1) == 0
                        dofs = [ux uy uz rx ry rz];
                    elseif keyopts(1) == 1
                        dofs = [ux uy uz];
                    else
                        msg = ['AnsysInput: Invalid keyopts for element type ', ...
                            elementName];
                        e = MException('MATLAB:bm_mfem:invalidKeyopts',msg);
                        throw(e);
                    end
                    
                case strcmp(elementName,'PLANE182')
                    dofs = [ux uy];
                    
                case strcmp(elementName,'SOLID185')
                    dofs = [ux uy uz];
                    
                case strcmp(elementName,'SOLID186')
                    dofs = [ux uy uz];
                    
                case strcmp(elementName,'SOLID187')
                    dofs = [ux uy uz];
                    
                otherwise
                    msg = ['AnsysInput: Available dofs for element ', ...
                        elementName, ' not defined.'];
                    e = MException('MATLAB:bm_mfem:undefinedElement',msg);
                    throw(e);
            end
        end
        
    end
    
end

