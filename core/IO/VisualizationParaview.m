classdef VisualizationParaview < handle
    properties (Access = private)
        model
        filename
    end
    
    methods
        
        % constructor
        function v = VisualizationParaview(femModel, filename)
            v.model = femModel;
            v.filename = filename;
        end
        
        function vtkWrite(v)
            % VTK Header using ASCII Format
            fid = fopen(v.filename, 'w');
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, 'VTK from Matlab\n');
            fprintf(fid, 'ASCII\nDATASET UNSTRUCTURED_GRID\n');
            
            % Coordinates of the Vertices
            nPoints = length(v.model.getAllNodes());
            fprintf(fid, ['POINTS ' num2str(nPoints) ' float\n']);
            for i = 1:nPoints
                fprintf(fid, [num2str(v.model.getNode(i).getCoords()) '\n']);
            end
            
            % Connectivity List
            nElements = length(v.model.getAllElements());
            fprintf(fid, ['CELLS ' num2str(nElements) ' ' num2str(nElements *5) '\n']);            
            for i = 1:nElements
                fprintf(fid , ['4\t' num2str(v.model.getElement(i).getNodes().getId() - 1) '\n']);
            end
            fprintf(fid, ['CELL_TYPES ' num2str(nElements) '\n']);
            for i = 1:nElements
                fprintf(fid, [num2str(9) '\n']);
            end
            
            % Solution Vector
            fprintf(fid, ['POINT_DATA ', num2str(nPoints), '\n']);
            fprintf(fid, 'VECTORS Displacement float\n');
            for i = 1 : nPoints
                fprintf(fid, [num2str(v.model.getNode(i).getDofValue('DISPLACEMENT_X'))...
                             '\t' num2str(v.model.getNode(i).getDofValue('DISPLACEMENT_Y'))...
                             '\t' num2str(v.model.getNode(i).getDofValue('DISPLACEMENT_Z')) '\n']);
            end
             
            fclose(fid);
            
        end
        
    end
end
