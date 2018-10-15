classdef Visualization < handle
    %VISUALIZATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        model
        panel
        
        deformedPlotLines
        contourPlot
        elemNumPlot
        nodeNumPlot
        loadPlotLines
        constrainPlotLines
        
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
                text(c(1), c(2), elemStr,'Interpreter','latex','Visible','off', ...
                    'HorizontalAlignment','center','FontSize',10)
            end
            
            % show node numbers
            for ii = 1:length(nodes)
                nodeStr = strcat('\sffamily\textcircled{',num2str(nodes(ii).getId),'}');
                text(nodes(ii).getX, nodes(ii).getY, nodeStr,'Interpreter','latex','Visible','off', ...
                    'HorizontalAlignment','left','FontSize',10,'Color','blue')
            end
            
            
        end
        
        function plotNumbering(visualization,numType,step)
            
            if nargin == 1
                step = 'end';
            end
            if strcmp(numType,'elements')
                visualization.clearElemNumbering;
            elseif strcmp(numType,'nodes')
                visualization.clearNodeNumbering;
            else
                warning('Available input for variable ''numType'' is ''nodes'' or ''elements''.')
            end
            
            elements = visualization.model.getAllElements;
            nodes = visualization.model.getAllNodes;
            disp_x = nodes.getDofValue('DISPLACEMENT_X',step);
            disp_y = nodes.getDofValue('DISPLACEMENT_Y',step);
            scaling = visualization.scalingFactor;
            
            
            x = nodes.getX' + scaling .* disp_x;
            y = nodes.getY' + scaling .* disp_y;
            
            if strcmp(numType,'elements')
                for ii = 1:length(elements)
                    c = elements(ii).barycenter;
                    elemStr = strcat('\sffamily\fbox{',num2str(elements(ii).getId),'}');
                    elemPlot = text(c(1), c(2), elemStr);
                    elemPlot.Interpreter = 'latex';
                    elemPlot.HorizontalAlignment = 'center';
                    elemPlot.FontSize = 10;
                    elemPlot.Tag = 'ElemNum';
                    visualization.elemNumPlot(ii) = elemPlot;
                end         
            elseif strcmp(numType,'nodes')
                for ii = 1:length(nodes)
                    nodeStr = num2str(nodes(ii).getId());
                    nodePlot = text(x(ii), y(ii), nodeStr);
                    nodePlot.Interpreter = 'latex';
                    nodePlot.HorizontalAlignment = 'left';
                    nodePlot.FontSize = 10;
                    nodePlot.Color = 'blue';
                    nodePlot.Tag = 'NodeNum';
                    visualization.nodeNumPlot(ii) = nodePlot;
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
                elPlot.Color = 'blue';
                visualization.deformedPlotLines(ii) = elPlot;
            end
            
        end
        
        function plotField(visualization,fieldType, step)
            if nargin == 2
                step = 'end';
            end

            visualization.clearContour;
            
            if ~isvalid(visualization.panel)
                visualization.panel = axes;
                visualization.panel.DataAspectRatio = [1 1 1]; % fix aspect ratio
            end

            hold on

            elements = visualization.model.getAllElements;
            nodes = visualization.model.getAllNodes;
            nElements = length(elements);
            [stressValue, element_connect] = computeElementStress(elements,nodes,step);
            
            disp_x = nodes.getDofValue('DISPLACEMENT_X',step);
            disp_y = nodes.getDofValue('DISPLACEMENT_Y',step);
            disp_absolute = sqrt(disp_x.^2+disp_y.^2);
            
            coords = zeros(length(nodes),2);
            scaling = visualization.scalingFactor;
            
            coords(:,1) = transpose(nodes.getX) + scaling .* disp_x;
            coords(:,2) = transpose(nodes.getY) + scaling .* disp_y;

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
                fpt = field(element_connect(ii,ord));

                fieldPlot = fill(xpt,ypt,fpt);
                fieldPlot.Parent = visualization.panel;
                visualization.contourPlot(ii) = fieldPlot;

            end
            hold off
            colormap parula(9)
            colorbar
        end
        
        function plotLoad(visualization, step)
            hold on
            
            visualization.clearLoad;
            
            nodes = visualization.model.getAllNodes;
            
            disp_x = nodes.getDofValue('DISPLACEMENT_X',step);
            disp_y = nodes.getDofValue('DISPLACEMENT_Y',step);
            scaling = visualization.scalingFactor;
            
            % scale load to fit model dimensions
            maxModelLength = max(abs([max(nodes.getX)-min(nodes.getX),max(nodes.getY)-min(nodes.getY)]));
            for i = 1:length(nodes)
                dofLoad(i,1) = nodes(i).getDof('DISPLACEMENT_X').getDofLoad();
                dofLoad(i,2) = nodes(i).getDof('DISPLACEMENT_Y').getDofLoad();
            end
            maxDofLoad = max(max(abs(dofLoad)));
            loadScaling = maxModelLength * 0.05 / maxDofLoad;
            
            x = nodes.getX' + scaling .* disp_x;
            y = nodes.getY' + scaling .* disp_y;

            for ii = 1:length(nodes)

                dofLoadX = nodes(ii).getDof('DISPLACEMENT_X').getDofLoad();
                dofLoadY = nodes(ii).getDof('DISPLACEMENT_Y').getDofLoad();
                plotX = x(ii) - dofLoadX * loadScaling;
                plotY = y(ii) - dofLoadY * loadScaling;
                
                if any([dofLoadX dofLoadY])
                    loadPlot = quiver(plotX,plotY,dofLoadX,dofLoadY);
                    loadPlot.Color = 'red';
                    loadPlot.MaxHeadSize = 0.6;
                    loadPlot.AutoScaleFactor = loadScaling;
                    loadPlot.LineWidth = 1;
                    loadPlot.Tag = 'load';
                    visualization.loadPlotLines(ii) = loadPlot;
                end
            end
            
            hold off
            
        end
        
        function plotConstrain(visualization, step)

            visualization.clearConstrain;
            
            nodes = visualization.model.getAllNodes;
            disp_x = nodes.getDofValue('DISPLACEMENT_X', step);
            disp_y = nodes.getDofValue('DISPLACEMENT_Y', step);
            scaling = visualization.scalingFactor;
            nodes = visualization.model.getAllNodes;
            
            % scale constrain to fit model dimension
            maxModelLength = max(abs([max(nodes.getX)-min(nodes.getX),max(nodes.getY)-min(nodes.getY)]));
            constrainScaling = maxModelLength * 0.01;
            
            x = nodes.getX' + scaling .* disp_x;
            y = nodes.getY' + scaling .* disp_y;
            
            for ii = 1:length(nodes)

                dofX = nodes(ii).getDof('DISPLACEMENT_X').isFixed();
                dofY = nodes(ii).getDof('DISPLACEMENT_Y').isFixed();
                
                if dofX == true
                    plotX = [x(ii);x(ii)-constrainScaling;x(ii)-constrainScaling;x(ii)];
                    plotY = [y(ii);y(ii)+constrainScaling;y(ii)-constrainScaling;y(ii)];
                    constrainPlot = line(plotX,plotY);
                    constrainPlot.Color = 'green';
                    constrainPlot.LineWidth = 1;
                    constrainPlot.Tag = 'constrain';
                    visualization.constrainPlotLines(ii,1) = constrainPlot;
                end
                    
                if dofY == true
                    plotX = [x(ii);x(ii)+constrainScaling;x(ii)-constrainScaling;x(ii)];
                    plotY = [y(ii);y(ii)-constrainScaling;y(ii)-constrainScaling;y(ii)];
                    constrainPlot = line(plotX,plotY);
                    constrainPlot.Color = 'green';
                    constrainPlot.LineWidth = 1;
                    constrainPlot.Tag = 'constrain';
                    visualization.constrainPlotLines(ii,2) = constrainPlot;
                end
            end
        end
        function plotLineData(visualization, selectedNodeIds, fieldType, step)
            if nargin == 3
                step = 'end';
            end
            
            hold on
            if nargin == 2
                step = 'end';
            end
            
            elements = visualization.model.getAllElements;
            nodes = visualization.model.getAllNodes;
            
            [stressValue, ~] = computeElementStress(elements,nodes,step);
            
            disp_x = nodes.getDofValue('DISPLACEMENT_X',step)';
            disp_y = nodes.getDofValue('DISPLACEMENT_Y',step)';
            disp_absolute = sqrt(disp_x.^2+disp_y.^2);
            
            if     strcmp(fieldType,'Select Field')
                
                warndlg('No Field selected. Please select Field.','Warning');
                field = 0;

            elseif     strcmp(fieldType,'displacement_x')
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
                
            else
                warning('fieldType not yet implemented.');
                field = 0;
            end
            
            if any(field)
                
                if length(selectedNodeIds) == 1
                    value = field(selectedNodeIds);
                    msg = [fieldType ' at node ' num2str(selectedNodeIds) ': ' num2str(value)];
                    msgbox(msg,'Result')
                else
                    for i = 1:length(selectedNodeIds)
                        selectedField(i) = field(selectedNodeIds(i));
                    end
                    
                    for i = 1:length(selectedNodeIds)
                        selectedNodes(i) = nodes(selectedNodeIds(i));
                    end
                    
                    if abs(selectedNodes(1).getX-selectedNodes(2).getX) < 1e-10 %ckeck if vertical
                        y = selectedNodes.getY';
                        value = selectedField';
                        %                 value = selectedNodes.getDofValue('DISPLACEMENT_X',step);
                        plot = sortrows([y value]);
                        
                        figure
                        line(plot(:,2),plot(:,1))
                        xlabel('\sigma_x');
                        ylabel('y');
                    else
                        x = selectedNodes.getX';
                        value = selectedField';
                        %                 value = selectedNodes.getDofValue('DISPLACEMENT_Y',step);
                        plot = sortrows([x value]);
                        figure
                        line(plot(:,1),plot(:,2))
                        xlabel('x');
                        ylabel(fieldType);
                        
                    end
                end
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
        
        function clearElemNumbering(v)
           delete(v.elemNumPlot);
           drawnow;
        end

        function clearNodeNumbering(v)
           delete(v.nodeNumPlot);
           drawnow;
        end
        
        function clearContour(v)
           delete(v.contourPlot);
           drawnow;
        end
        
        function clearLoad(v)
           delete(v.loadPlotLines);
           drawnow;
        end
        
        function clearConstrain(v)
           delete(v.constrainPlotLines);
           drawnow;
        end
        
    end
    
end

