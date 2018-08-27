classdef VisualizationParaviewXML < handle
    properties (Access = private)
        model
        filename
    end
    
    methods
        
        % constructor
        function obj = VisualizationParaviewXML(femModel, filename)
            obj.model = femModel;
            obj.filename = filename;
        end
        
        function vtkWrite(obj)
            % Header
            
            % VTK Header using ASCII Format
            dofs = obj.model.getAllNodes.getDofMap.keys;
            nTsteps = length(obj.model.getNode(1).getDofValue(dofs{1},'all'));
            
            if ~exist(obj.filename, 'dir')
                mkdir(obj.filename);
            else
                warning('Folder already exists');
            end
            
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
                
                % Displacements
                fprintf(fid, ['\t\t\t<PointData>\n'...
                    '\t\t\t\t<DataArray type="Float32" Name="Displacement" '...
                    'NumberOfComponents="3" format="ascii" RangeMin="0" RangeMax="0">\n']);
                
                disp = zeros(nPoints,3);
                
                if  ismember( 'DISPLACEMENT_X', dofs)
                    disp(:,1) = obj.model.getAllNodes.getDofValue('DISPLACEMENT_X',tstep);
                end
                if  ismember( 'DISPLACEMENT_Y', dofs)
                    disp(:,2) = obj.model.getAllNodes.getDofValue('DISPLACEMENT_Y',tstep);
                end
                if  ismember( 'DISPLACEMENT_Z', dofs)
                    disp(:,3) = obj.model.getAllNodes.getDofValue('DISPLACEMENT_Z',tstep);
                end
                
                for j = 1 : nPoints
                    fprintf(fid, '\t\t\t\t\t%f \t %f \t %f \n', disp(j,1), disp(j,2), disp(j,3));
                end
                fprintf(fid, ['\t\t\t\t</DataArray>\n'...
                    '\t\t\t\t<DataArray type="Float32" Name="Rotation" '...
                    'NumberOfComponents="3" format="ascii" RangeMin="0" RangeMax="0">\n']);
                
                %Rotations
                rot = zeros(nPoints,3);
                if  ismember( 'ROTATION_X', dofs)
                    rot(:,1) = obj.model.getAllNodes.getDofValue('ROTATION_X',tstep);
                end
                if  ismember( 'ROTATION_Y', dofs)
                    rot(:,2) = obj.model.getAllNodes.getDofValue('ROTATION_Y',tstep);
                end
                if  ismember( 'ROTATION_Z', dofs)
                    rot(:,3) = obj.model.getAllNodes.getDofValue('ROTATION_Z',tstep);
                end
                
                for j = 1 : nPoints
                    fprintf(fid, '\t\t\t\t\t%f \t %f \t %f \n', rot(j,1), rot(j,2), rot(j,3));
                end
                
                fprintf(fid, ['\t\t\t\t</DataArray>\n'...
                    '\t\t\t</PointData>\n']);
                
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


