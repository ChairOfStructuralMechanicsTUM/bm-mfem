classdef VisualizationGUI < handle
    %VISUALIZATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        model
        
        deformedPlotLines
        
        scalingFactor
    end
    
    methods
        
        % constructor
        function v = VisualizationGUI(femModel)
            v.model = femModel;
            v.scalingFactor = 1;
        end
        
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
                nodeStr = strcat('\sffamily\textcircled{',num2str(nodes(ii).getId),'}');
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
            hold off
        end
        
        function plotField(visualization,fieldType, step)
            
            hold on
            if nargin == 2
                step = 'end';
            end            
            
            elements = visualization.model.getAllElements;
            nodes = visualization.model.getAllNodes;
            nElements = length(elements);
            [stressValue, element_connect] = computeElementStress(elements,nodes);
            coords = zeros(length(nodes),3);

            for i = 1:length(nodes)
                coords(i,1) = nodes(i).getX();
                coords(i,2) = nodes(i).getY();
                coords(i,3) = nodes(i).getZ();
            end
            
            if     fieldType == 'sigma_xx'
                    field = stressValue(1,:);
                    
            elseif fieldType == 'sigma_yy'
                    field = stressValue(2,:);
                    
            elseif fieldType == 'sigma_xy'
                    field = stressValue(3,:);
            end
            
            for ii = 1:nElements
                
                ord = elements(ii).drawOrder();
                
                xpt = coords(element_connect(ii,ord),1);
                ypt = coords(element_connect(ii,ord),2);
                zpt = coords(element_connect(ii,ord),3);
                fpt = field(element_connect(ii,ord));

                fill3(xpt,ypt,zpt,fpt,'Tag',fieldType);

            end
                
            colorbar
            hold off
        end
        
        function close(visualization)
           close(ancestor(visualization.panel,'Figure'))
        end
        
        function setScaling(v, scalingFactor)
            v.scalingFactor = scalingFactor;
        end
        
        function clearDeformed(v)
           delete(v.deformedPlotLines);
           drawnow;
        end
        
    end
    
end

