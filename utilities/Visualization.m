classdef Visualization < handle
    %VISUALIZATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        model
        panel
        
        deformedPlotLines
        
        scalingFactor
    end
    
    methods
        
        % constructor
        function v = Visualization(femModel)
            v.model = femModel;
            v.panel = axes;
            close(ancestor(v.panel,'Figure'))
            v.scalingFactor = 1;
        end
        
        function plotUndeformed(visualization, textOutput)
            
            if nargin ~= 2
                textOutput = false;
            end
            
            if ~isvalid(visualization.panel)
                visualization.panel = axes;
                visualization.panel.DataAspectRatio = [1 1 1]; % fix aspect ratio
            end
            
            elements = visualization.model.getAllElements;
            nodes = visualization.model.getAllNodes;
            
            % plot the undeformed system
            for ii = 1:length(elements)
                elPlot = elements(ii).draw;
                elPlot.Parent = visualization.panel;
            end
            
            if textOutput
                % show element numbers
                if(all(nodes.getDimension == 3))
                    for ii = 1:length(elements)
                        c = elements(ii).barycenter;
                        elemStr = strcat('\sffamily\fbox{',num2str(elements(ii).getId),'}');
                        if length(c) == 2
                            text(c(1), c(2), elemStr,'Interpreter','latex', ...
                                'HorizontalAlignment','center','FontSize',10)
                        elseif length(c) == 3
                            text(c(1), c(2), c(3), elemStr,'Interpreter','latex', ...
                                'HorizontalAlignment','center','FontSize',10)
                        end
                    end
                else
                    for ii = 1:length(elements)
                        c = elements(ii).barycenter;
                        elemStr = strcat('\sffamily\fbox{',num2str(elements(ii).getId),'}');
                        text(c(1), c(2), elemStr,'Interpreter','latex', ...
                            'HorizontalAlignment','center','FontSize',10)
                    end
                end

                % show node numbers
                if(all(nodes.getDimension == 3))
                    for ii = 1:length(nodes)
                        nodeStr = strcat('\sffamily\textcircled{',num2str(nodes(ii).getId),'}');
                        text(nodes(ii).getX, nodes(ii).getY, nodes(ii).getZ, nodeStr,'Interpreter','latex', ...
                            'HorizontalAlignment','left','FontSize',10,'Color','blue')
                    end
                else
                    for ii = 1:length(nodes)
                        nodeStr = strcat('\sffamily\textcircled{',num2str(nodes(ii).getId),'}');
                        text(nodes(ii).getX, nodes(ii).getY, nodeStr,'Interpreter','latex', ...
                            'HorizontalAlignment','left','FontSize',10,'Color','blue')
                    end
                end
            end
            
        end
        
        function plotDeformed(visualization, step)
            if nargin == 1
                step = 'end';
            end
            
            if ~isvalid(visualization.panel)
                visualization.panel = axes;
                visualization.panel.DataAspectRatio = [1 1 1]; % fix aspect ratio
            end
            
            visualization.clearDeformed;
            
            elements = visualization.model.getAllElements;
            nElements = length(elements);
            visualization.deformedPlotLines = matlab.graphics.primitive.Line.empty;
            
            % plot the deformed system
            for ii = 1:nElements
                elPlot = elements(ii).drawDeformed(step, visualization.scalingFactor);
                elPlot.Parent = visualization.panel;
                elPlot.Color = 'red';
                visualization.deformedPlotLines(ii) = elPlot;
            end
            
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

