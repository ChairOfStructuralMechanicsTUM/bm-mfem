classdef AnsysInput < ModelIO
    %AnsysInput Read model input from ANSYS
    %   Detailed explanation goes here
    
    properties (Access = private)
        ansysExecutable
        printOutput = false
    end
    
    methods
        
        function obj = AnsysInput(file, ansysExecutable)
            % file - full file path with filetype as extensions
            % exaample: file='F:\ALMA_simulations\AnsysModels\plate.inp';
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
                msg = 'AnsysInput: Wrong number of input arguments';
                err = MException('MATLAB:bm_mfem:invalidArguments',msg);
                throw(err);
            end
            
            obj@ModelIO(super_args{:});
            
            if exist(ansysExecutable, 'file') == 2
                obj.ansysExecutable = ansysExecutable;
            else
                msg = ['AnsysInput: Ansys executable could not be found at' ...
                    ' specified location ', ansysExecutable];
                err = MException('MATLAB:bm_mfem:ansysNotFound',msg);
                throw(err);
            end
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
            
            if obj.printOutput
                disp('Ansys input file selected')
                disp(['folder: ' folder])
                disp(['filename: ' fileName])
                disp(['filetype: ' extension])
            end
            
            model = FemModel;
            data=obj.runAnsys(obj.ansysExecutable,folder,fileName,extension);
            
            % Read restriction on the nodes
            data.nodeRest = obj.readRestrictions();
            % Read loads on the nodes
            data.nodeLoads = obj.readLoads();
            % Read record 5 of ANSYS to get the position of the entries
            % of each node with the matrices
            nodeEquiv = obj.readRecord_5();
            % Read the coordinates of each node
            data.nodeList = obj.readCoord();
            % Order the coordinates acoording to the record 5
            [~, order] = ismember(nodeEquiv.', data.nodeList(:,1));
            nodesC = data.nodeList(order,:);
            data.nodesOrderByDofs=nodesC(:,1)';
            
            % Create objects "Node" and assign them to a model
            for i = 1 : size(nodesC,1)
                id =  nodesC(i,1);
                x  =  nodesC(i,2);
                y  =  nodesC(i,3);
                z  =  nodesC(i,4);
                model.addNewNode(id,x,y,z);
            end
            data.numNodes = length(model.getAllNodes());
            
            % Read available element data
            [data.elementsOfModel,data.nodeElementList,data.nodeConnectivity] = AnsysInput.readElements();
            
            % Assign dofs to node
            for i = 1 :size(data.elementsOfModel,1)
                nodeData=cell2mat(data.nodeElementList(i));
                dofs = AnsysInput.getDofsOfAnsysElements(data.elementsOfModel{i,1}, data.elementsOfModel{i,2});
                nodeData(isnan(nodeData)) = [];
                for j = 1:length(nodeData)
                    n = model.getNode(nodeData(j));
                    n.addDof(dofs);
                end
            end
            % Order dofs by nodes (bm-mfem style; the dummy element takes
            % care of reordering the matrices)
            dofArray = arrayfun(@(node) node.getDofArray, model.getAllNodes(), 'UniformOutput', false);
            dofArray = [dofArray{:}];
            for ii = 1:length(dofArray)
                dofArray(ii).setId(ii);
            end
            
            % Create the dummy element
            e = model.addNewElement('DummyElement',1,model.getAllNodes);
            
            % Set the matrices
            systemSize = size(data.Mansys,1);
            Mdiag = spdiags(data.Mansys,0);
            M = data.Mansys + data.Mansys.' - spdiags(Mdiag(:),0,systemSize,systemSize);
            Cdiag = spdiags(data.Cansys,0);
            C = data.Cansys + data.Cansys.' - spdiags(Cdiag(:),0,systemSize,systemSize);
            Kdiag = spdiags(data.Kansys,0);
            K = data.Kansys + data.Kansys.' - spdiags(Kdiag(:),0,systemSize,systemSize);
            e.setMatrices(M, C, K);
            
            % Set the dof order
            e.setDofOrder(data.nodesOrderByDofs);
            
            % Set restrictions
            e.setDofRestrictions(data.nodeRest);
            
            % Set loads
            for ii=1:length(data.nodeLoads{1})
                n = model.getNode(data.nodeLoads{1}(ii));
                n.setDofLoad(data.nodeLoads{2}{ii}, data.nodeLoads{3}(ii));
            end
            
            % Set node connectivity
            e.setNodeConnectivity(data.nodeConnectivity);
            
            % Add everything to a model part
            model.addNewModelPart('ANSYS_model', ...
                model.getAllNodes().getId(), model.getAllElements().getId());
            
            try
                rmdir('DataAnsys', 's')
            catch
            end
            
        end
        
        function setPrintOutput(obj, print)
            if isa(print, 'logical')
                obj.printOutput = print;
            else
                msg = 'AnsysInput: Input argument must be a logical';
                err = MException('MATLAB:bm_mfem:invalidArguments',msg);
                throw(err);
            end
        end
        
    end
    
    methods (Static)
        
        function data = runAnsys(ansysExecutable,folder,file,extension)
            % make directory for files
            if ~exist('DataAnsys','dir'); mkdir('DataAnsys'); end
            
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
            % Create external file with loads
            fprintf(fidl,'/output,DataAnsys/nodeLoads,txt \r\n');
            fprintf(fidl,'FLIST,ALL,\r\n');
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
            
            % Read matrices
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
                            aux1 = [aux1 A(j,k)]; %#ok<AGROW>
                        end
                        aux2 = [aux2 aux1]; %#ok<AGROW>
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
            fclose(fid);
            
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
                    nodeNum = [nodeNum str2double(A.textdata(i,j))]; %#ok<AGROW>
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
            
            for ii = 1:length(A{2})
                A{2}{ii} = AnsysInput.replaceANSYSDofName(A{2}{ii});
            end
            
        end
        
        function A = readLoads()
            % function A = readLoads()
            %
            % Function   : readLoads
            %
            % Description: This function gets the loads on coordinates from
            %              ANSYS
            %
            % Parameters :
            %
            % Return     : A cell with loads.
            %               A{1}: node numbers
            %               A{2}: resticted dof name
            %               A{3}: real values the dof is restricted to
            %               A{4}: imaginary values the dof is restricted to
            
            fid=fopen('DataAnsys/nodeLoads.txt');
            for j = 1:13
                fgetl(fid); % Read/discard line.
            end
            
            A = textscan(fid,'%u%s%f%f');
            
            for ii = 1:length(A{2})
                A{2}{ii} = AnsysInput.replaceANSYSDofName(A{2}{ii});
            end
            
        end
        
        function name = replaceANSYSDofName(ansysName)
            switch ansysName
                case {'UX', 'FX'}
                    name = 'DISPLACEMENT_X';
                case {'UY', 'FY'}
                    name = 'DISPLACEMENT_Y';
                case {'UZ', 'FZ'}
                    name = 'DISPLACEMENT_Z';
                case 'ROTX'
                    name = 'ROTATION_X';
                case 'ROTY'
                    name = 'ROTATION_Y';
                case 'ROTZ'
                    name = 'ROTATION_Z';
                otherwise
                    msg = ['AnsysInput: Unknown ANSYS dof name ', ansysName];
                    err = MException('MATLAB:bm_mfem:invalidArguments',msg);
                    throw(err);
            end
        end
        
    end
    
end

