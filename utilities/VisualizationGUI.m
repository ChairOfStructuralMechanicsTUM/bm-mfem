classdef VisualizationGUI < handle
    %VISUALIZATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        model
        
        deformedPlotLines
        fieldPlotLines
        
        scalingFactor
    end
    
    methods
        
        % constructor
        function v = VisualizationGUI(femModel)
            v.model = femModel;
            v.scalingFactor = 1;
        end 

%         function plotUndeformed(visualization)
%             hold on
%             tic;
%             elements = visualization.model.getAllElements;
%             nodes = visualization.model.getAllNodes;
%             
%             for ii = 1:length(elements)
%             connect(ii,1:4) = elements(ii).getNodes.getId();
%             end
%             
%             coords(:,1) = transpose(nodes.getX);
%             coords(:,2) = transpose(nodes.getY);
%             
%             % plot the undeformed system
%             for ii = 1:length(elements)
%                 
%                 ord = elements(ii).drawOrder();
%                 
%                 x = coords(connect(ii,ord),1);
%                 y = coords(connect(ii,ord),2);
% 
%                 plot(x,y,'Tag','undeformed','Color','blue');
% 
%             end
%             toc
%             % show element numbers
%             for ii = 1:length(elements)
%                 c = elements(ii).barycenter;
%                 elemStr = strcat('\sffamily\fbox{',num2str(elements(ii).getId),'}');
%                 text(c(1), c(2), elemStr,'Interpreter','latex', ...
%                     'HorizontalAlignment','center','FontSize',10,'Tag','ElemNum','Visible','off')
%             end
%             % show node numbers
%             for ii = 1:length(nodes)
% %                 nodeStr = strcat('\sffamily\textcircled{',num2str(nodes(ii).getId),'}');
%                 nodeStr = num2str(nodes(ii).getId());
%                 text(nodes(ii).getX, nodes(ii).getY, nodeStr,'Interpreter','latex', ...
%                     'HorizontalAlignment','left','FontSize',10,'Color','blue','Tag','NodeNum','Visible','off')
%             end
%             axis equal
%             hold off
%             
%         end
        function plotUndeformed(visualization)
            hold on
            elements = visualization.model.getAllElements;
            nodes = visualization.model.getAllNodes;
            
            % plot the undeformed system
            for ii = 1:length(elements)
                elPlot = elements(ii).draw;
                elPlot.Tag = 'undeformed';
            end
            % show element numbers
            for ii = 1:length(elements)
                c = elements(ii).barycenter;
                elemStr = strcat('\sffamily\fbox{',num2str(elements(ii).getId),'}');
                text(c(1), c(2), elemStr,'Interpreter','latex', ...
                    'HorizontalAlignment','center','FontSize',10,'Tag','ElemNum','Visible','off')
            end
            % show node numbers
            for ii = 1:length(nodes)
%                 nodeStr = strcat('\sffamily\textcircled{',num2str(nodes(ii).getId),'}');
                nodeStr = num2str(nodes(ii).getId());
                text(nodes(ii).getX, nodes(ii).getY, nodeStr,'Interpreter','latex', ...
                    'HorizontalAlignment','left','FontSize',10,'Color','blue','Tag','NodeNum','Visible','off')
            end
            axis equal
            hold off
            
        end
        
        function plotDeformed(visualization, step)
            hold on
            if nargin == 1
                step = 'end';
            end
            
            visualization.clearDeformed;
            
            elements = visualization.model.getAllElements;
            nElements = length(elements);
            visualization.deformedPlotLines = matlab.graphics.primitive.Line.empty;
            
            % plot the deformed system
            for ii = 1:nElements
                elPlot = elements(ii).drawDeformed(step, visualization.scalingFactor);
                elPlot.Color = 'red';
                elPlot.Visible = 'off';
                elPlot.Tag = 'deformed';
                visualization.deformedPlotLines(ii) = elPlot;
            end
            
            data = guidata(findobj('Tag','figure1'));
            data.deformedScaling = visualization.scalingFactor;
            guidata(findobj('Tag','figure1'),data);
            
            hold off
        end
        
        function plotField(visualization,fieldType, step)
            
            hold on
            if nargin == 2
                step = 'end';
            end
            
            visualization.clearField;
            
            elements = visualization.model.getAllElements;
            nodes = visualization.model.getAllNodes;
            nElements = length(elements);
            [stressValue, element_connect] = computeElementStress(elements,nodes);
            
            disp_x = nodes.getDofValue('DISPLACEMENT_X');
            disp_y = nodes.getDofValue('DISPLACEMENT_Y');
            disp_absolute = sqrt(disp_x.^2+disp_y.^2);
            
            coords = zeros(length(nodes),3);
            scaling = visualization.scalingFactor;
            
            coords(:,1) = transpose(nodes.getX) + scaling .* disp_x;
            coords(:,2) = transpose(nodes.getY) + scaling .* disp_y;
            coords(:,3) = nodes.getZ;
            
            if     strcmp(fieldType,'displacement_x')
                    field = disp_x;
                    
            elseif     strcmp(fieldType,'displacement_y')
                    field = disp_y;
                    
            elseif     strcmp(fieldType,'displacement_absolute')
                    field = disp_absolute;
                
            elseif     strcmp(fieldType,'sigma_xx')
                    field = stressValue(1,:);
                    
            elseif strcmp(fieldType,'sigma_yy')
                    field = stressValue(2,:);
                    
            elseif strcmp(fieldType,'sigma_xy')
                    field = stressValue(3,:);
                    
            elseif strcmp(fieldType,'prin_I')
                    field = stressValue(4,:);
                    
            elseif strcmp(fieldType,'prin_II')
                    field = stressValue(5,:);
            
            elseif strcmp(fieldType,'vm_stress')
                    field = stressValue(6,:);
            end
            
            for ii = 1:nElements
                
                ord = elements(ii).drawOrder();
                
                xpt = coords(element_connect(ii,ord),1);
                ypt = coords(element_connect(ii,ord),2);
                zpt = coords(element_connect(ii,ord),3);
                fpt = field(element_connect(ii,ord));

                fieldPlot = fill3(xpt,ypt,zpt,fpt);
                fieldPlot.Tag = fieldType;
                visualization.fieldPlotLines(ii) = fieldPlot;
            end
            
            data = guidata(findobj('Tag','figure1'));
            data.fieldScaling = visualization.scalingFactor;
            guidata(findobj('Tag','figure1'),data);
                
            colorbar
            hold off
        end
        
        function close(visualization)
           close(ancestor(visualization.panel,'Figure'))
        end
        
        function setScaling(v, scalingFactor)
            v.scalingFactor = scalingFactor;
        end
        
        function scalingFactor = getScaling(v)
            scalingFactor = v.scalingFactor;
        end
        
        function clearDeformed(v)
           delete(v.deformedPlotLines);
           drawnow;
        end
        
        function clearField(v)
           delete(v.fieldPlotLines);
           drawnow;
        end
        
    end
    
end
