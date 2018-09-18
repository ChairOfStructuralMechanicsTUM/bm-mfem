classdef VisualizationParaviewXML < handle
    properties (Access = private)
        model
        filename
        dofname
    end
    
    methods
        
        % constructor
        function obj = VisualizationParaviewXML(femModel, filename,varargin)
            obj.model = femModel;
            obj.filename = filename;
            obj.dofname = varargin;
        end
        
        function vtkWrite(obj)
            % Header
            
            % VTK Header using ASCII Format
            dofs = obj.model.getAllNodes.getDofMap.keys;
            nTsteps = length(obj.model.getNode(1).getDofValue(dofs{1},'all'));
            
            if exist(obj.filename, 'dir')
                rmdir(obj.filename, 's');
            end
            mkdir(obj.filename);
            
            for tstep = 1 : nTsteps
                fid = fopen(fullfile(obj.filename,[obj.filename, '_', num2str(tstep), '.vtu']), 'w');
                
                nPoints = length(obj.model.getAllNodes());
                nElements = length(obj.model.getAllElements());
                
                fprintf(fid, ['<VTKFile type="UnstructuredGrid" version="1.0" '...
                    'byte_order="LittleEndian" header_type="UInt64">\n',...
                    '\t<UnstructuredGrid>\n', ...
                    '\t\t<Piece NumberOfPoints="%d" NumberOfCells="%d">\n '], nPoints, nElements);
                
                % Points
                
                fprintf(fid, ['\t\t\t<Points>\n',...
                    '\t\t\t\t<DataArray type="Float32" Name="Points" NumberOfComponents="3" '...
                    'format="ascii" RangeMin="0" RangeMax="0">\n']);
                
                % Coordinates of the Vertices
                for i = 1 : nPoints
                    fprintf(fid, ['\t\t\t\t\t' num2str(obj.model.getNode(i).getCoords()) '\n']);
                end
                
                fprintf(fid, ['\t\t\t\t</DataArray>\n'...
                    '\t\t\t</Points>\n']);
                
                % Cells
                
                fprintf(fid, ['\t\t\t<Cells>\n',...
                    '\t\t\t\t<DataArray type="Int64" Name="connectivity" format="ascii" ',...
                    'RangeMin="0" RangeMax="0">\n']);
                
                % Connectivity List
                for i = 1 : nElements
                    fprintf(fid , ['\t\t\t\t\t' num2str(obj.model.getElement(i).getNodes().getId() - 1) '\n']);
                end
                
                fprintf(fid,['\t\t\t\t</DataArray>\n'...
                    '\t\t\t\t<DataArray type="UInt8" Name="types" format="ascii" RangeMin="0" RangeMax="0">\n']);
                
                % Cell Type
                for i = 1 : nElements
                    switch class(obj.model.getElement(i))
                        case 'BarElement2d2n'
                            fprintf(fid, ['\t\t\t\t\t' num2str(4) '\n']);
                        case 'BarElement3d2n'
                            fprintf(fid, ['\t\t\t\t\t' num2str(4) '\n']);
                        case 'BeamElement3d2n'
                            fprintf(fid, ['\t\t\t\t\t' num2str(4) '\n']);
                        case 'ConcentratedMassElement3d1n'
                            fprintf(fid, ['\t\t\t\t\t' num2str(1) '\n']);
                        case 'SpringDamperElement3d2n'
                            fprintf(fid, ['\t\t\t\t\t' num2str(4) '\n']);
                        case 'ReissnerMindlinElement3d4n'
                            fprintf(fid, ['\t\t\t\t\t' num2str(9) '\n']);
                        case 'ShellElement3d4n'
                            fprintf(fid, ['\t\t\t\t\t' num2str(9) '\n']);
                        case 'DiscreteKirchhoffElement3d4n'
                            fprintf(fid, ['\t\t\t\t\t' num2str(9) '\n']);
                        case 'QuadrilateralElement2d4n'
                            fprintf(fid, ['\t\t\t\t\t' num2str(4) '\n']);
                        case 'HexahedronElement3d8n'
                            fprintf(fid, ['\t\t\t\t\t' num2str(11) '\n']);
                    end
                end
                
                % Offsets
                fprintf(fid,['\t\t\t\t</DataArray>\n'...
                    '\t\t\t\t<DataArray type="Int64" Name="offsets" format="ascii" RangeMin="0" RangeMax="0">\n']);
                counter = 0;
                for i = 1 : nElements
                    switch class(obj.model.getElement(i))
                        case 'BarElement2d2n'
                            counter = counter + 2;
                            fprintf(fid, ['\t\t\t\t\t' num2str(counter) '\n']);
                        case 'BarElement3d2n'
                            counter = counter + 2;
                            fprintf(fid, ['\t\t\t\t\t' num2str(counter) '\n']);
                        case 'BeamElement3d2n'
                            counter = counter + 2;
                            fprintf(fid, ['\t\t\t\t\t' num2str(counter) '\n']);
                        case 'ConcentratedMassElement3d1n'
                            counter = counter + 1;
                            fprintf(fid, ['\t\t\t\t\t' num2str(counter) '\n']);
                        case 'SpringDamperElement3d2n'
                            counter = counter + 2;
                            fprintf(fid, ['\t\t\t\t\t' num2str(counter) '\n']);
                        case 'ReissnerMindlinElement3d4n'
                            counter = counter + 4;
                            fprintf(fid, ['\t\t\t\t\t' num2str(counter) '\n']);
                        case 'ShellElement3d4n'
                            counter = counter + 4;
                            fprintf(fid, ['\t\t\t\t\t' num2str(counter) '\n']);
                        case 'DiscreteKirchhoffElement3d4n'
                            counter = counter + 4;
                            fprintf(fid, ['\t\t\t\t\t' num2str(counter) '\n']);
                        case 'QuadrilateralElement2d4n'
                            counter = counter + 4;
                            fprintf(fid, ['\t\t\t\t\t' num2str(counter) '\n']);
                        case 'HexahedronElement3d8n'
                            counter = counter + 8;
                            fprintf(fid, ['\t\t\t\t\t' num2str(counter) '\n']);
                    end
                end
                
                fprintf(fid,['\t\t\t\t</DataArray>\n'...
                    '\t\t\t</Cells>\n']);
                
                % Point Data
                
                fprintf(fid, '\t\t\t<PointData>\n');
                
                for i = 1 : length(obj.dofname)
                    
                    fprintf(fid, ['\t\t\t\t<DataArray type="Float32" Name="',obj.dofname{i} , '" '...
                        'NumberOfComponents="3" format="ascii" RangeMin="0" RangeMax="0">\n']);
                    
                    sol = zeros(nPoints,3);
                    
                    if  ismember([obj.dofname{i} , '_X'], dofs)
                        sol(:,1) = obj.model.getAllNodes.getDofValue([obj.dofname{i} , '_X'], tstep);
                    end
                    if  ismember([obj.dofname{i} , '_Y'], dofs)
                        sol(:,2) = obj.model.getAllNodes.getDofValue([obj.dofname{i} , '_Y'], tstep);
                    end
                    if  ismember([obj.dofname{i} , '_Z'], dofs)
                        sol(:,3) = obj.model.getAllNodes.getDofValue([obj.dofname{i} , '_Z'], tstep);
                    end
                    
                    for j = 1 : nPoints
                        fprintf(fid, '\t\t\t\t\t%f \t %f \t %f \n', sol(j,1), sol(j,2), sol(j,3));
                    end    
                    fprintf(fid, '\t\t\t\t</DataArray>\n');
                end
                
                fprintf(fid, '\t\t\t</PointData>\n');
                
                % Appendix
                
                fprintf(fid,['\t\t</Piece>\n'...
                    '\t</UnstructuredGrid>\n'...
                    '</VTKFile>']);
                fclose(fid);
            end
        end
        
        function pvdWrite(obj)
            
            fid = fopen([obj.filename, '.pvd'], 'w');
            
            fprintf(fid, ['<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" '...
                'header_type="UInt64">\n'...
                '\t<Collection>\n']);
            
            dofs = obj.model.getAllNodes.getDofMap.keys;
            nTsteps = length(obj.model.getNode(1).getDofValue(dofs{1},'all'));
            
            for i = 1 : nTsteps
                
                fprintf(fid,['\t\t<DataSet 	timestep="', num2str(i),'" group="" part="0" '...
                    'file="', obj.filename,'/',obj.filename, '_', num2str(i), '.vtu', '"/>\n']);
            end
            
            fprintf(fid, ['\t</Collection>\n'...
                '</VTKFile>']);
            
            obj.vtkWrite();
            
        end
        
    end
end


