classdef VisualizationParaview < handle
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
            % VTK Header using ASCII Format
            fid = fopen(obj.filename, 'w');
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, 'VTK from Matlab\n');
            fprintf(fid, 'ASCII\nDATASET UNSTRUCTURED_GRID\n');
            
            % Coordinates of the Vertices
            nPoints = length(obj.model.getAllNodes());
            fprintf(fid, ['POINTS ' num2str(nPoints) ' float\n']);
            for i = 1:nPoints
                fprintf(fid, [num2str(obj.model.getNode(i).getCoords()) '\n']);
            end
            
            % Connectivity List
            nElements = length(obj.model.getAllElements());
            fprintf(fid, ['CELLS ' num2str(nElements) ' ' num2str(nElements *5) '\n']);
            for i = 1:nElements
                fprintf(fid , ['4\t' num2str(obj.model.getElement(i).getNodes().getId() - 1) '\n']);
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
                        fprintf(fid, [num2str(4) '\n']);
                    case 'SpringDamperElement3d2n'
                        fprintf(fid, [num2str(9) '\n']);
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
            
            
            nTsteps = length(obj.model.getNode(i).getDofValue('DISPLACEMENT_Z','all'));
            % Solution Vector
            
            fprintf(fid, ['POINT_DATA ', num2str(nPoints), '\n']);
            fprintf(fid, ['FIELD Output ', num2str((nTsteps-1) * length(obj.dofname)), '\n']);
            
            for i = 1 : length(obj.dofname)
                
                for j = 2:nTsteps
                    
                    fprintf(fid, [[obj.dofname{i} num2str(j-1)] ' ' num2str(3) ' ' num2str(nPoints) ' float\n']);
                    
                    for k = 1 : nPoints
                        
                        if obj.dofname{i}(end) == 'X'
                            
                            fprintf(fid, [num2str(obj.model.getNode(k).getDofValue(obj.dofname{i},j))  '\t' num2str(0) '\t' num2str(0) '\n']);
                        
                        elseif obj.dofname{i}(end) == 'Y'
                            
                            fprintf(fid, [num2str(0) '\t' num2str(obj.model.getNode(k).getDofValue(obj.dofname{i},j))  '\t' num2str(0) '\n']);
                        
                        else
                            
                            fprintf(fid, [num2str(0)  '\t' num2str(0) '\t' num2str(obj.model.getNode(k).getDofValue(obj.dofname{i},j)) '\n']);
                        
                        end
                    end
                end
            end
            
            
            fclose(fid);
            
        end
    end
end

