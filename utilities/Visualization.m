classdef Visualization < handle
    %VISUALIZATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        model
        panel
%         existPanel
    end
    
    methods
        
        % constructor
        function v = Visualization(femModel)
            v.model = femModel;
            v.panel = axes;
            close(ancestor(v.panel,'Figure'))
        end
        
        function plotUndeformed(visualization)
            
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
            
            % show element numbers
            for ii = 1:length(elements)
                c = elements(ii).barycenter;
                elemStr = strcat('\sffamily\fbox{',num2str(elements(ii).getId),'}');
                text(c(1), c(2), elemStr,'Interpreter','latex', ...
                    'HorizontalAlignment','center','FontSize',10)
            end
            
            % show node numbers
            for ii = 1:length(nodes)
                nodeStr = strcat('\sffamily\textcircled{',num2str(nodes(ii).getId),'}');
                text(nodes(ii).getX, nodes(ii).getY, nodeStr,'Interpreter','latex', ...
                    'HorizontalAlignment','left','FontSize',10,'Color','blue')
            end
            
            
        end
        
        function plotDeformed(visualization)
            
            if ~isvalid(visualization.panel)
                visualization.panel = axes;
                visualization.panel.DataAspectRatio = [1 1 1]; % fix aspect ratio
            end
            
            elements = visualization.model.getAllElements;
            
            % plot the deformed system
            for ii = 1:length(elements)
                elPlot = elements(ii).drawDeformed;
                elPlot.Parent = visualization.panel;
                elPlot.Color = 'red';
            end
            
        end
        
        function close(visualization)
           close(ancestor(visualization.panel,'Figure'))
        end
        
    end
    
end

