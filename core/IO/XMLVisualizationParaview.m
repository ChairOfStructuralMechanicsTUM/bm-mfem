classdef XMLVisualizationParaview < handle
    properties (Access = private)
        model
        filename
        dofname
    end
    
    methods
        
        % constructor
        function obj = VisualizationParaview(femModel, filename, varargin)
            obj.model = femModel;
            obj.filename = filename;
            obj.dofname = varargin;
        end
        
        function vtkWrite(obj)
            
            dofs = obj.model.getAllNodes.getDofMap.keys;
            nTsteps = length(obj.model.getNode(1).getDofValue(dofs{1},'all'));
            
            for tstep = 1:nTsteps
                
            % VTK Header using ASCII Format
            fid = fopen([num2str(tstep) obj.filename], 'w');
            
            nPoints = length(obj.model.getAllNodes());
             nElements = length(obj.model.getAllElements());
             
             
            fprintf(fid, ['<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n',...
                             '\t<UnstructuredGrid>\n', ...
                             '\t\t<Piece NumberOfPoints="%d" NumberOfCells="%d">\n ', ...
                             '\t\t\t<PointData>\n'], nPoints, nElements);
            
            
            
            % Coordinates of the Vertices
            for i = 1:nPoints
                fprintf(fid, [num2str(obj.model.getNode(i).getCoords()) '\n']);
            end
            
            % Connectivity List
            nElements = length(obj.model.getAllElements());
            fprintf(fid, 'CELLS %d %d \n', nElements, nElements*(length(obj.model.getElement(1).getNodes)+1));
            
            for i = 1:nElements
                fprintf(fid , ['%d \t' num2str(obj.model.getElement(i).getNodes().getId() - 1) '\n'], length(obj.model.getElement(1).getNodes));
            end
            
            fprintf(fid, ['CELL_TYPES ' num2str(nElements) '\n']);
            for i = 1:nElements
                switch class(obj.model.getElement(i))
                    case 'BarElement2d2n'
                        fprintf(fid, [num2str(4) '\n']);
                    case 'BarElement3d2n'
                        fprintf(fid, [num2str(4) '\n']);
                    case 'BeamElement3d2n'
                        fprintf(fid, [num2str(4) '\n']);
                    case 'ConcentratedMassElement3d1n'
                        fprintf(fid, [num2str(1) '\n']);
                    case 'SpringDamperElement3d2n'
                        fprintf(fid, [num2str(4) '\n']);
                    case 'ReissnerMindlinElement3d4n'
                        fprintf(fid, [num2str(9) '\n']);
                    case 'ShellElement3d4n'
                        fprintf(fid, [num2str(9) '\n']);
                    case 'DiscreteKirchhoffElement3d4n'
                        fprintf(fid, [num2str(9) '\n']);
                    case 'QuadrilateralElement2d4n'
                        fprintf(fid, [num2str(4) '\n']);
                    case 'HexahedronElement3d8n'
                        fprintf(fid, [num2str(11) '\n']);
                end
            end
            
            
            % Solution Vector           
            fprintf(fid, 'POINT_DATA %d \n', nPoints);
            fprintf(fid, 'FIELD Output %d \n', 2);
            

                disp = zeros(nPoints,3);
                rot = zeros(nPoints,3);

            
%             for i = 2:nTsteps
    
                
                if find(ismember(obj.model.getAllNodes.getDofMap.keys, 'DISPLACEMENT_X'))
                    disp(:,1) = obj.model.getAllNodes.getDofValue('DISPLACEMENT_X',tstep);
                end
                if find(ismember(obj.model.getAllNodes.getDofMap.keys, 'DISPLACEMENT_Y'))
                    disp(:,2) = obj.model.getAllNodes.getDofValue('DISPLACEMENT_Y',tstep);
                end
                if find(ismember(obj.model.getAllNodes.getDofMap.keys, 'DISPLACEMENT_Z'))
                    disp(:,3) = obj.model.getAllNodes.getDofValue('DISPLACEMENT_Z',tstep);
                end
                if find(ismember(obj.model.getAllNodes.getDofMap.keys, 'ROTATION_X'))
                    rot(:,1) = obj.model.getAllNodes.getDofValue('ROTATION_X',tstep);
                end
                if find(ismember(obj.model.getAllNodes.getDofMap.keys, 'ROTATION_Y'))
                    rot(:,2) = obj.model.getAllNodes.getDofValue('ROTATION_Y',tstep);
                end
                if find(ismember(obj.model.getAllNodes.getDofMap.keys, 'ROTATION_Z'))
                    rot(:,3) = obj.model.getAllNodes.getDofValue('ROTATION_Z',tstep);
                end

                fprintf(fid, 'Displacement %d %d float \n' , 3, nPoints);
                for j = 1 : nPoints                    
                    fprintf(fid, '%f \t %f \t %f \n', disp(j,1), disp(j,2), disp(j,3));                    
                end

                fprintf(fid, 'Rotation %d %d float \n', 3, nPoints);
                for j = 1 : nPoints                    
                    fprintf(fid, '%f \t %f \t %f \n', rot(j,1), rot(j,2), rot(j,3));                    
                end

%             end
 
            fclose(fid);
            end
        end
    end
end

