classdef VisualizationGUI < handle
    %VISUALIZATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        model
        
        deformedPlotLines
        fieldPlotLines
        loadPlotLines
        constrainPlotLines
        numPlot
        
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
            axis equal
            hold off
            
        end
        function plotNumbering(visualization,numType,state)
            hold on
            
            try
                visualization.clearNumbering;
            end
            
            nodes = visualization.model.getAllNodes;
            disp_x = nodes.getDofValue('DISPLACEMENT_X');
            disp_y = nodes.getDofValue('DISPLACEMENT_Y');
            scaling = visualization.scalingFactor;
            nodes = visualization.model.getAllNodes;
            
            if strcmp(state,'undeformed')
                x = nodes.getX';
                y = nodes.getY';
            elseif strcmp(state,'deformed')
                x = nodes.getX' + scaling .* disp_x;
                y = nodes.getY' + scaling .* disp_y;
            else
                warning('Error plotting Load Boundary Condition.');
            end
            if strcmp(numType,'elements')
                for ii = 1:length(elements)
                    c = elements(ii).barycenter;
                    elemStr = strcat('\sffamily\fbox{',num2str(elements(ii).getId),'}');
                    elemPlot = text(c(1), c(2), elemStr);
                    elemPlot.Interpreter = 'latex';
                    elemPlot.HorizontalAlignment = 'center';
                    elemPlot.FontSize = 10;
                    elemPlot.Tag = 'ElemNum';
                end
            else strcmp(numType,'nodes')
                for ii = 1:length(nodes)
                    nodeStr = num2str(nodes(ii).getId());
                    nodePlot = text(x(ii), y(ii), nodeStr);
                    nodePlot.Interpreter = 'latex';
                    nodePlot.HorizontalAlignment = 'left';
                    nodePlot.FontSize = 10;
                    nodePlot.Color = 'blue';
                    nodePlot.Tag = 'NodeNum';
                    visualization.numPlot(ii) = nodePlot;
                end
            end   
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
        
        function plotField(visualization,fieldType,step)
            
            hold on
            if nargin == 2
                step = 'end';
            end
            
            try
            visualization.clearField;
            end
            
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
            
            if     strcmp(fieldType,'Select Field')

                for ii = 1:length(elements)
                    
                    ord = elements(ii).drawOrder();
                    
                    x = coords(element_connect(ii,ord),1);
                    y = coords(element_connect(ii,ord),2);
                    
                    defPlot = line(x,y);
                    defPlot.Tag = fieldType;
                    defPlot.Color = 'black';
                    visualization.fieldPlotLines(ii) = defPlot;
                end
                colorbar off
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
            end
            
            if ~strcmp(fieldType,'Select Field')
                for ii = 1:nElements

                    ord = elements(ii).drawOrder();

                    x = coords(element_connect(ii,ord),1);
                    y = coords(element_connect(ii,ord),2);
                    f = field(element_connect(ii,ord));

                    fieldPlot = fill(x,y,f);
                    fieldPlot.Tag = fieldType;
                    visualization.fieldPlotLines(ii) = fieldPlot;
                end
                colorbar
            end
            
            data = guidata(findobj('Tag','figure1'));
            data.scaling = visualization.scalingFactor;
            guidata(findobj('Tag','figure1'),data);
            
            hold off
        end
        
        function plotLineData(visualization, selectedNodeIds, fieldType, step)
            
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
                        xlabel(fieldType);
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

        function plotLoad(visualization,state)
            hold on
            
            try
                visualization.clearLoad;
            end
            
            nodes = visualization.model.getAllNodes;
            maxModelLength = max(abs([max(nodes.getX)-min(nodes.getX),max(nodes.getY)-min(nodes.getY)]));
            data = guidata(findobj('Tag','figure1'));
            maxDofLoad = data.maxDofLoad;
            loadScaling = maxModelLength * 0.05 / maxDofLoad;
            disp_x = nodes.getDofValue('DISPLACEMENT_X');
            disp_y = nodes.getDofValue('DISPLACEMENT_Y');
            scaling = visualization.scalingFactor;
            
            if strcmp(state,'undeformed')
                x = nodes.getX';
                y = nodes.getY';
            elseif strcmp(state,'deformed')
                x = nodes.getX' + scaling .* disp_x;
                y = nodes.getY' + scaling .* disp_y;
            else
                warning('Error plotting Load Boundary Condition.');
            end

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
        
        function plotConstrain(visualization, state)
            hold on
            
%             try
            delete(findobj(gcf,'Tag','constrain'));
%             if ~isempty(findobj(gcf,'Tag','constrain'))
%                 visualization.clearConstrain;
%             end
            
            nodes = visualization.model.getAllNodes;
            disp_x = nodes.getDofValue('DISPLACEMENT_X');
            disp_y = nodes.getDofValue('DISPLACEMENT_Y');
            scaling = visualization.scalingFactor;
            nodes = visualization.model.getAllNodes;
            maxModelLength = max(abs([max(nodes.getX)-min(nodes.getX),max(nodes.getY)-min(nodes.getY)]));
            constrainScaling = maxModelLength * 0.01;
            
            if strcmp(state,'undeformed')
                x = nodes.getX';
                y = nodes.getY';
            elseif strcmp(state,'deformed')
                x = nodes.getX' + scaling .* disp_x;
                y = nodes.getY' + scaling .* disp_y;
            else
                warning('Error plotting Load Boundary Condition.');
            end

            
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
        
        function clearLoad(v)
            delete(v.loadPlotLines);
            drawnow;
        end
        
        function clearConstrain(v)
            delete(v.constrainPlotLines);
            drawnow;
        end
        
        function a = getConstrain(v)
            a = v.constrainPlotLines;
            drawnow;
        end
        
        function a = getLoad(v)
            a = v.loadPlotLines;
            drawnow;
        end
        
        function clearField(v)
           delete(v.fieldPlotLines);
           drawnow;
        end
        
        function clearNumbering(v)
           delete(v.numPlot);
           drawnow;
        end
        
    end
    
end

