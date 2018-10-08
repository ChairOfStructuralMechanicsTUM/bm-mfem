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
        
        function plotField(visualization,fieldType, step)
            if nargin == 2
                step = 'end';
            end

            if ~isvalid(visualization.panel)
                visualization.panel = axes;
                visualization.panel.DataAspectRatio = [1 1 1]; % fix aspect ratio
            end
            hold on
            
            
            elements = visualization.model.getAllElements;
            nodes = visualization.model.getAllNodes;
            nElements = length(elements);
            [stressValue, ~, element_connect] = computeElementStress(elements,nodes,step);
            coords = zeros(length(nodes),3);

            for i = 1:length(nodes)
                coords(i,1) = nodes(i).getX();
                coords(i,2) = nodes(i).getY();
                coords(i,3) = nodes(i).getZ();
            end
            
            if     strcmp(fieldType,'sigma_xx')
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

                fill3(xpt,ypt,zpt,fpt);

            end
                
            colorbar
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
