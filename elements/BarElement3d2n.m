classdef BarElement3d2n < BarElement
    %BARELEMENT3D2N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    properties (Access = private, Constant = true)
       dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
    end
    
    methods
        % constructor
        function barElement3d2n = BarElement3d2n(id, nodeArray, material, crossSectionArea)
           barElement3d2n@BarElement(id, material, crossSectionArea);
           
           if (length(nodeArray) == 2 && isa(nodeArray,'Node'))
               barElement3d2n.nodeArray = nodeArray;
           else
               error('problem with the nodes in element %d', id);
           end
           
           barElement3d2n.addDofs(barElement3d2n.dofNames);
           
           barElement3d2n.length = computeLength(barElement3d2n.nodeArray(1).getCoords, ...
               barElement3d2n.nodeArray(2).getCoords);
           
        end
        
        % getter functions
        
        % member functions
        function stiffnessMatrix = computeLocalStiffnessMatrix(barElement)
            dist = barElement.nodeArray(2).getCoords - barElement.nodeArray(1).getCoords;
            x21 = dist(1);
            y21 = dist(2);
            z21 = dist(3);
            stiffnessMatrix = [x21*x21 x21*y21 x21*z21 -x21*x21 -x21*y21 -x21*z21; ...
                x21*y21 y21*y21 y21*z21 -x21*y21 -y21*y21 -y21*z21; ...
                x21*z21 y21*z21 z21*z21 -x21*z21 -y21*z21 -z21*z21; ...
                -x21*x21 -x21*y21 -x21*z21 x21*x21 x21*y21 x21*z21; ...
                -x21*y21 -y21*y21 -y21*z21 x21*y21 y21*y21 y21*z21; ...
                -x21*z21 -y21*z21 -z21*z21 x21*z21 y21*z21 z21*z21];
            factor = (barElement.getMaterial().getParameterValue('YOUNGS_MODULUS') ...
                * barElement.crossSectionArea) ...
                / (barElement.length^3);
            stiffnessMatrix = factor * stiffnessMatrix;
        end
        
    end
    
end

