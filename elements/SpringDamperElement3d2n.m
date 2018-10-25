classdef SpringDamperElement3d2n < LinearElement
    %SPRINGDAMPERELEMENT3D2N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        length0 %initial length
    end
    
    methods
        % constructor
        function obj = SpringDamperElement3d2n(id, nodeArray)
            
            requiredPropertyNames = cellstr(["ELEMENTAL_STIFFNESS", "ELEMENTAL_DAMPING"]);
            
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 2 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class constructor
            obj@LinearElement(super_args{:});
            obj.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
                        
        end
        
        function initialize(obj)
            obj.localSystem = obj.computeLocalSystem();
            obj.length0 = computeLength(obj.nodeArray(1).getCoords, ...
                    obj.nodeArray(2).getCoords);
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            stiffness = obj.getPropertyValue('ELEMENTAL_STIFFNESS');
            
            stiffnessMatrix = sparse([1 1 4 4],[1 4 4 1],[stiffness -stiffness stiffness -stiffness],6,6);
            tMat = obj.getTransformationMatrix;
            stiffnessMatrix = tMat' * stiffnessMatrix * tMat;
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            massMatrix = sparse(6,6);
        end
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            damping = obj.getPropertyValue('ELEMENTAL_DAMPING');
            
            dampingMatrix = sparse([1 1 4 4],[1 4 4 1],[damping -damping damping -damping],6,6);
            tMat = obj.getTransformationMatrix;
            dampingMatrix = tMat' * dampingMatrix * tMat;
        end
        
        function forceVector = computeLocalForceVector(obj)
           forceVector = sparse(1,6);
%            stiffness = element.getPropertyValue('ELEMENTAL_STIFFNESS');
%            nodes = element.getNodes;
%            
%            tMat = element.getTransformationMatrix;
%            
%            disp = zeros(1,3);
%            disp(1) = nodes(2).getDofValue('DISPLACEMENT_X','end') - nodes(1).getDofValue('DISPLACEMENT_X','end');
%            disp(2) = nodes(2).getDofValue('DISPLACEMENT_Y','end') - nodes(1).getDofValue('DISPLACEMENT_Y','end');
%            disp(3) = nodes(2).getDofValue('DISPLACEMENT_Z','end') - nodes(1).getDofValue('DISPLACEMENT_Z','end');
%            localDisp = tMat(1:3,1:3)' * disp';
%            
%            forceVector(1) = - stiffness * localDisp(1);
%            forceVector(4) = stiffness * localDisp(1);
%            
%            forceVector = tMat' * forceVector';
%            forceVector = forceVector';
        end
        
        function dofs = getDofList(obj)
            dofs([1 4]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2 5]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
            dofs([3 6]) = obj.nodeArray.getDof('DISPLACEMENT_Z');
        end
        
        function vals = getValuesVector(obj, step)
            vals = zeros(1,6);
            
            vals([1 4]) = obj.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 5]) = obj.nodeArray.getDofValue('DISPLACEMENT_Y',step);
            vals([3 6]) = obj.nodeArray.getDofValue('DISPLACEMENT_Z',step);
        end
        
        function vals = getFirstDerivativesVector(obj, step)
            vals = zeros(1,6);
            
            [~, vals([1 4]), ~] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, vals([2 5]), ~] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, vals([3 6]), ~] = obj.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(obj, step)
            vals = zeros(1,6);

            [~, ~, vals([1 4])] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, ~, vals([2 5])] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, ~, vals([3 6])] = obj.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
        end
        
        function update(obj)
            obj.length0 = computeLength(obj.nodeArray(1).getCoords, ...
                    obj.nodeArray(2).getCoords);
        end
        
    end
    
end

