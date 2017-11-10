classdef BeamElement3d2n < LinearElement
    %BEAMELEMENT3D2N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        length0 %initial length
    end
    
    methods
        % constructor
        function beamElement = BeamElement3d2n(id, nodeArray)
            
            requiredPropertyNames = cellstr(["IY", "IZ", "IT", ...
                "YOUNGS_MODULUS", "SHEAR_MODULUS", "CROSS_SECTION"]);
            
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 2 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class constructor
            beamElement@LinearElement(super_args{:});
            beamElement.dofNames = cellstr(["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]);
                       
        end
        
        function initialize(element)
            element.localSystem = element.computeLocalSystem();
            element.length0 = computeLength(element.nodeArray(1).getCoords, ...
                element.nodeArray(2).getCoords);
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(element)
           E = element.getPropertyValue('YOUNGS_MODULUS');
           G = element.getPropertyValue('SHEAR_MODULUS');
           A = element.getPropertyValue('CROSS_SECTION');
           Iy = element.getPropertyValue('IY');
           Iz = element.getPropertyValue('IZ');
           It = element.getPropertyValue('IT');
           L = element.length0;
           L2 = L*L;
           L3 = L2*L;
           
           stiffnessMatrix = sparse(12,12);
           
           stiffnessMatrix(1,1) = E*A / L;
           stiffnessMatrix(1,7) = - stiffnessMatrix(1,1);
           stiffnessMatrix(7,7) = stiffnessMatrix(1,1);
           stiffnessMatrix(7,1) = - stiffnessMatrix(1,1);
           
           EIz = E*Iz;
           stiffnessMatrix(2,2) = 12 * EIz / L3;
           stiffnessMatrix(2,8) = - stiffnessMatrix(2,2);
           stiffnessMatrix(8,8) = stiffnessMatrix(2,2);
           stiffnessMatrix(8,2) = - stiffnessMatrix(2,2);
           stiffnessMatrix(6,6) = 4 * EIz / L;
           stiffnessMatrix(12,12) = stiffnessMatrix(6,6);
           stiffnessMatrix(2,6) = 6 * EIz / L2;
           stiffnessMatrix(2,12) = stiffnessMatrix(2,6);
           stiffnessMatrix(6,2) = stiffnessMatrix(2,6);
           stiffnessMatrix(12,2) = stiffnessMatrix(2,6);
           stiffnessMatrix(6,8) = - stiffnessMatrix(2,6);
           stiffnessMatrix(8,6) = - stiffnessMatrix(2,6);
           stiffnessMatrix(8,12) = - stiffnessMatrix(2,6);
           stiffnessMatrix(12,8) = - stiffnessMatrix(2,6);
           stiffnessMatrix(6,12) = 2 * EIz / L;
           stiffnessMatrix(12,6) = stiffnessMatrix(6,12);
           
           EIy = E*Iy;
           stiffnessMatrix(3,3) = 12 * EIy / L3;
           stiffnessMatrix(3,9) = - stiffnessMatrix(3,3);
           stiffnessMatrix(9,9) = stiffnessMatrix(3,3);
           stiffnessMatrix(9,3) = - stiffnessMatrix(3,3);
           stiffnessMatrix(5,5) = 4 * EIy / L;
           stiffnessMatrix(11,11) = stiffnessMatrix(5,5);
           stiffnessMatrix(3,5) = - 6 * EIy / L2;
           stiffnessMatrix(3,11) = stiffnessMatrix(3,5);
           stiffnessMatrix(5,3) = stiffnessMatrix(3,5);
           stiffnessMatrix(11,3) = stiffnessMatrix(3,5);
           stiffnessMatrix(5,9) = - stiffnessMatrix(3,5);
           stiffnessMatrix(9,5) = - stiffnessMatrix(3,5);
           stiffnessMatrix(9,11) = - stiffnessMatrix(3,5);
           stiffnessMatrix(11,9) = - stiffnessMatrix(3,5);
           stiffnessMatrix(5,11) = 2 * EIy / L;
           stiffnessMatrix(11,5) = stiffnessMatrix(5,11);
           
           stiffnessMatrix(4,4) = G*It / L;
           stiffnessMatrix(10,10) = stiffnessMatrix(4,4);
           stiffnessMatrix(10,4) = - stiffnessMatrix(4,4);
           stiffnessMatrix(4,10) = - stiffnessMatrix(4,4);
           
           tMat = element.getTransformationMatrix;
           stiffnessMatrix = tMat' * stiffnessMatrix * tMat;
        end
        
        function massMatrix = computeLocalMassMatrix(element)
            massMatrix = zeros(12);
        end
        
        function dampingMatrix = computeLocalDampingMatrix(element)
            dampingMatrix = zeros(12);
%            damping = element.getPropertyValue('ELEMENTAL_DAMPING');
%            dampingMatrix = zeros(6);
%            dampingMatrix(1,1) = damping;
%            dampingMatrix(1,4) = - damping;
%            dampingMatrix(4,4) = damping;
%            dampingMatrix(4,1) = - damping;
%            
%            tMat = element.getTransformationMatrix;
%            dampingMatrix = tMat' * dampingMatrix * tMat;           
        end
        
        function forceVector = computeLocalForceVector(element)
           forceVector = zeros(1,12);
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
        
        function update(springDamperElement)
            springDamperElement.length0 = computeLength(springDamperElement.nodeArray(1).getCoords, ...
                    springDamperElement.nodeArray(2).getCoords);
        end
        
    end
    
    methods (Access = private)
        
    end
    
end

