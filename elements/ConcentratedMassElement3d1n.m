classdef ConcentratedMassElement3d1n < Element
    %CONCENTRATEDMASSELEMENT3D1N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        % constructor
        function element = ConcentratedMassElement3d1n(id, nodeArray)
            requiredPropertyNames = cellstr(["ELEMENTAL_MASS", "VOLUME_ACCELERATION"]);
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 1 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class constructor
            element@Element(super_args{:});
            element.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
                        
        end
        
        function initialize(element)
            
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(element)
            stiffnessMatrix = sparse(3,3);
        end
        
        function forceVector = computeLocalForceVector(ele)
           mass = ele.getPropertyValue('ELEMENTAL_MASS');
           volume_acceleration = ele.getPropertyValue('VOLUME_ACCELERATION');
           forceVector = volume_acceleration .* mass;
        end
        
        function massMatrix = computeLocalMassMatrix(element)
           mass = element.getPropertyValue('ELEMENTAL_MASS');
           massMatrix = sparse([1 2 3],[1 2 3],[mass mass mass],3,3);        
        end
        
        function dampingMatrix = computeLocalDampingMatrix(element)
            dampingMatrix = sparse(3,3);
        end
        
        function dofs = getDofList(element)
            dofs(1) = element.nodeArray.getDof('DISPLACEMENT_X');
            dofs(2) = element.nodeArray.getDof('DISPLACEMENT_Y');
            dofs(3) = element.nodeArray.getDof('DISPLACEMENT_Z');
        end
        
        function vals = getValuesVector(element, step)
            vals = sparse(1,3);
            
            vals(1) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals(2) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
            vals(3) = element.nodeArray.getDofValue('DISPLACEMENT_Z',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = sparse(1,3);
            
            [~, vals(1), ~] = element.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, vals(2), ~] = element.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, vals(3), ~] = element.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = sparse(1,3);

            [~, ~, vals(1)] = element.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, ~, vals(2)] = element.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, ~, vals(3)] = element.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
        end
        
        
        
        function update(element)
            
        end
        
        function c = barycenter(element)
            nodes = element.getNodes;
            c = nodes(1).getCoords;
        end
        
        function plot = draw(element)
            plot = line(element.nodeArray.getX, element.nodeArray.getY, ...
                'color','b','Marker','d','MarkerFaceColor','b');
        end
        
        function plot = drawDeformed(element)
            plot = line(element.nodeArray.getX + ...
                element.nodeArray.getDofValue('DISPLACEMENT_X'), ...
                element.nodeArray.getY + ...
                element.nodeArray.getDofValue('DISPLACEMENT_Y'), ...
                'color','r','Marker','d','MarkerFaceColor','r');
        end
        
    end
    
end

