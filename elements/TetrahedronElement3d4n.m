classdef TetrahedronElement3d4n < Element  %Class Tetrahedron to be implemented
    %TETRAHEDRONELEMENT3D4N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        % constructor
        function obj = TetrahedronElement3d4n(id, nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO", "DENSITY"]);
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 4 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class constructor
            obj@Element(super_args{:});
            obj.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
            
        end
        
        % Compute Barycenter of a Tetrahedron
        function c = barycenter(obj)
            c = 1/4 * (obj.nodeArray(1).getCoords() ...
                + obj.nodeArray(2).getCoords() ...
                + obj.nodeArray(3).getCoords() ...
                + obj.nodeArray(4).getCoords());
        end
        
        %Initialization
        function initialize(obj)
            % check, if the element's nodes lie on the same plane
            x = obj.nodeArray.getX();
            y = obj.nodeArray.getY();
            z = obj.nodeArray.getZ();
            if all(abs(x - x(1)) < 100*eps) ...
                    || all(abs(y - y(1)) < 100*eps) ...
                    || all(abs(z - z(1)) < 100*eps)
                msg = ['TetrahedronElement3d4n: Element ', ...
                    num2str(obj.getId), 's nodes are on the same plane.'];
                e = MException('MATLAB:bm_mfem:elementInvalid',msg);
                throw(e);
            end
        end
        
        % member functions
        function [Jdet, NfDiff] = computeShapeFunction(obj)
            
            % Calculation of Determinant of Jacobian Matrix
            xDiff = zeros(4,4);
            yDiff = zeros(4,4);
            zDiff = zeros(4,4);
            
            for i = 1:4
                for j = 1:4
                    xDiff(i,j) = obj.nodeArray(i).getX - obj.nodeArray(j).getX;
                    yDiff(i,j) = obj.nodeArray(i).getY - obj.nodeArray(j).getY;
                    zDiff(i,j) = obj.nodeArray(i).getZ - obj.nodeArray(j).getZ;
                end
            end
            
            Jdet = xDiff(2,1)*(yDiff(2,3)*zDiff(3,4)-yDiff(3,4)*zDiff(2,3)) ...
                + xDiff(3,2)*(yDiff(3,4)*zDiff(1,2)-yDiff(1,2)*zDiff(3,4)) ...
                + xDiff(4,3)*(yDiff(1,2)*zDiff(2,3)-yDiff(2,3)*zDiff(1,2));
            
            if nargout == 2
                % Calculation of the derivatives of the Shape Functions
                NfDiffX = [yDiff(4,2)*zDiff(3,2)-yDiff(3,2)*zDiff(4,2) ...
                    yDiff(3,1)*zDiff(4,3)-yDiff(3,4)*zDiff(1,3) ...
                    yDiff(2,4)*zDiff(1,4)-yDiff(1,4)*zDiff(2,4) ...
                    yDiff(1,3)*zDiff(2,1)-yDiff(1,2)*zDiff(3,1)];  % Derivatives w.r.t. x
                NfDiffY = [xDiff(3,2)*zDiff(4,2)-xDiff(4,2)*zDiff(3,2) ...
                    xDiff(4,3)*zDiff(3,1)-xDiff(1,3)*zDiff(3,4) ...
                    xDiff(1,4)*zDiff(2,4)-xDiff(2,4)*zDiff(1,4) ...
                    xDiff(2,1)*zDiff(1,3)-xDiff(3,1)*zDiff(1,2)];  % Derivatives w.r.t. y
                NfDiffZ = [xDiff(4,2)*yDiff(3,2)-xDiff(3,2)*yDiff(4,2) ...
                    xDiff(3,1)*yDiff(4,3)-xDiff(3,4)*yDiff(1,3) ...
                    xDiff(2,4)*yDiff(1,4)-xDiff(1,4)*yDiff(2,4) ...
                    xDiff(1,3)*yDiff(2,1)-xDiff(1,2)*yDiff(3,1)];  % Derivatives w.r.t. z

                NfDiff = [NfDiffX; NfDiffY; NfDiffZ];
            end
            
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            
            [Jdet, NfDiff] = computeShapeFunction(obj);
            
            % Compute B-Matrix
            Be = [NfDiff(1,1),0,0,            NfDiff(1,2),0,0,            NfDiff(1,3),0,0,            NfDiff(1,4),0,0;...
                0,NfDiff(2,1),0,            0,NfDiff(2,2),0,            0,NfDiff(2,3),0,            0,NfDiff(2,4),0;...
                0,0,NfDiff(3,1),            0,0,NfDiff(3,2),            0,0,NfDiff(3,3),            0,0,NfDiff(3,4);...
                NfDiff(2,1),NfDiff(1,1),0,  NfDiff(2,2),NfDiff(1,2),0,  NfDiff(2,3),NfDiff(1,3),0,  NfDiff(2,4),NfDiff(1,4),0;...
                0,NfDiff(3,1),NfDiff(2,1),  0,NfDiff(3,2),NfDiff(2,2),  0,NfDiff(3,3),NfDiff(2,3),  0,NfDiff(3,4),NfDiff(2,4);...
                NfDiff(3,1),0,NfDiff(1,1),   NfDiff(3,2),0,NfDiff(1,2), NfDiff(3,3),0,NfDiff(1,3),  NfDiff(3,4),0,NfDiff(1,4)];
            
            % Compute Emat-Matrix (for isotropic material only)
            Emodul = obj.getPropertyValue('YOUNGS_MODULUS');
            PoissonRatio = obj.getPropertyValue('POISSON_RATIO');
            
            a=Emodul*(1-PoissonRatio)/( (1-2*PoissonRatio)*(1+PoissonRatio) );
            b=Emodul*PoissonRatio/( (1-2*PoissonRatio)*(1+PoissonRatio) );
            c=Emodul/(2*(1+PoissonRatio));
            
            Emat=[a,b,b,0,0,0; ...
                b,a,b,0,0,0; ...
                b,b,a,0,0,0; ...
                0,0,0,c,0,0; ...
                0,0,0,0,c,0; ...
                0,0,0,0,0,c];
            
            % Compute Element Stiffness Matrix
            stiffnessMatrix=(1/(6*Jdet))*(transpose(Be)*(Emat*Be));
            
        end
        
        
        function massMatrix = computeLocalMassMatrix(obj)
            
            rho = obj.getPropertyValue('DENSITY');
            Jdet = computeShapeFunction(obj);
            V = 1/6*Jdet;
            
            massMatrix = rho*V/20 * ...
                [2 0 0 1 0 0 1 0 0 1 0 0;...
                0 2 0 0 1 0 0 1 0 0 1 0;...
                0 0 2 0 0 1 0 0 1 0 0 1;...
                1 0 0 2 0 0 1 0 0 1 0 0;...
                0 1 0 0 2 0 0 1 0 0 1 0;...
                0 0 1 0 0 2 0 0 1 0 0 1;...
                1 0 0 1 0 0 2 0 0 1 0 0;...
                0 1 0 0 1 0 0 2 0 0 1 0;...
                0 0 1 0 0 1 0 0 2 0 0 1;...
                1 0 0 1 0 0 1 0 0 2 0 0;...
                0 1 0 0 1 0 0 1 0 0 2 0;...
                0 0 1 0 0 1 0 0 1 0 0 2];
            
        end
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            % up to now only Rayleigh-Damping
            eProperties = obj.getProperties;
            dampingMatrix = sparse(12,12);
            
            if (eProperties.hasValue('RAYLEIGH_ALPHA'))
                alpha = eProperties.getValue('RAYLEIGH_ALPHA');
                dampingMatrix = dampingMatrix + alpha * element.computeLocalMassMatrix;
            end
            
            if (eProperties.hasValue('RAYLEIGH_BETA'))
                beta = eProperties.getValue('RAYLEIGH_BETA');
                dampingMatrix = dampingMatrix + beta * element.computeLocalStiffnessMatrix;
            end
            
        end
        
        function dofs = getDofList(obj)
            dofs([1 4 7 10]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2 5 8 11]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
            dofs([3 6 9 12]) = obj.nodeArray.getDof('DISPLACEMENT_Z');
        end
        
        function vals = getValuesVector(obj, step)
            vals = zeros(1,12);
            
            vals([1 4 7 10]) = obj.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 5 8 11]) = obj.nodeArray.getDofValue('DISPLACEMENT_Y',step);
            vals([3 6 9 12]) = obj.nodeArray.getDofValue('DISPLACEMENT_Z',step);
        end
        
        function vals = getFirstDerivativesVector(obj, step)
            vals = zeros(1,12);
            
            [~, vals([1 4 7 10]), ~] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, vals([2 5 8 11]), ~] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, vals([3 6 9 12]), ~] = obj.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(obj, step)
            vals = zeros(1,12);
            
            [~, ~, vals([1 4 7 10])] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, ~, vals([2 5 8 11])] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, ~, vals([3 6 9 12])] = obj.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
        end
        
        function pl = draw(obj)
            x = [obj.nodeArray(1).getX, obj.nodeArray(2).getX, ...
                obj.nodeArray(3).getX, obj.nodeArray(1).getX,...
                obj.nodeArray(4).getX, obj.nodeArray(2).getX, ...
                obj.nodeArray(3).getX, obj.nodeArray(4).getX];
            
            y = [obj.nodeArray(1).getY, obj.nodeArray(2).getY, ...
                obj.nodeArray(3).getY, obj.nodeArray(1).getY,...
                obj.nodeArray(4).getY, obj.nodeArray(2).getY, ...
                obj.nodeArray(3).getY, obj.nodeArray(4).getY];
            
            z = [obj.nodeArray(1).getZ, obj.nodeArray(2).getZ, ...
                obj.nodeArray(3).getZ, obj.nodeArray(1).getZ,...
                obj.nodeArray(4).getZ, obj.nodeArray(2).getZ, ...
                obj.nodeArray(3).getZ, obj.nodeArray(4).getZ];
            
            pl = line(x,y,z);
        end
        
        function pl = drawDeformed(obj, step, scaling)
            x = [obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_X', step),...
                obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_X', step)];
            
            y = [obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Y', step),...
                obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Y', step)];
            
            z = [obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(2).getZ + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(3).getZ + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(4).getZ + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Z', step),...
                obj.nodeArray(2).getZ + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(3).getZ + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(4).getZ + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Z', step)];
            
            pl = line(x,y,z);
            
        end
        
        function update(obj)
        end
        
        function F = computeLocalForceVector(obj)
            F = zeros(1,12);
        end
        
        
    end
end

