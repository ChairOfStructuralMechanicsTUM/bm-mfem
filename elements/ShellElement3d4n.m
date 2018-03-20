classdef ShellElement3d4n < QuadrilateralElement 
    %SHELLELEMENT3D4N A quadrilateral shell element
    %   Detailed explanation goes here
    
    properties (Access = private)
    end 
    
    methods
        % Constructor
        function shellElement3d4n = ShellElement3d4n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO", ...
                                             "THICKNESS", "NUMBER_GAUSS_POINT", ...
                                             "DENSITY", "SHEAR_CORRECTION_FACTOR"]);
                                         
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 4 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class contructor
            shellElement3d4n@QuadrilateralElement(super_args{:});
            shellElement3d4n.dofNames = cellstr([ ...
                "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]);
        end
        
        %Initialization
        function initialize(shellElement3d4n)
            shellElement3d4n.lengthX = computeLength(shellElement3d4n.nodeArray(1).getCoords, ...
                shellElement3d4n.nodeArray(2).getCoords);
            
            shellElement3d4n.lengthY = computeLength(shellElement3d4n.nodeArray(1).getCoords, ...
                shellElement3d4n.nodeArray(4).getCoords);
            
            checkConvexity(shellElement3d4n);
        end
        
        function responseDoF = getResponseDofArray(shellElement3d4n, step)
           
            responseDoF = zeros(24,1);
            for itNodes = 1:1:6
                nodalDof = shellElement3d4n.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';
                
                for itDof = 3:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
        end
        
        function [] = computeShapeFunction(shellElement3d4n, xi, eta)
            
            M = [(1-xi)/2    (1+xi)/2    ;   (1-eta)/2     (1+eta)/2];
            N = [1/2 - 3/4*xi + xi^3/4     ,  1/4 - xi/4 - xi^2/4 + xi^3/4       , 1/2 + 3/4 * xi - xi^3/4      ,  -1/4 - xi/4 +xi^2/4 + xi^3/4      
                    1/2 - 3/4*eta + eta^3/4 , 1/4 - eta/4 - eta^2/4 + eta^3/4 , 1/2 + 3/4 * eta - eta^3/4 ,  -1/4 - eta/4 +eta^2/4 + eta^3/4]
            
            MN = sparse(2,16); 
            
            MN(1,1)   =  M(1,1) * N(2,1); 
            MN(2,2)   =  M(2,1) * N(1,1);
            MN(1,3)   = -M(1,1) * N(2,2); 
            MN(2,4)   =  M(2,1) * N(1,2);
            MN(1,5)   =  M(1,2) * N(2,1);
            MN(2,6)   =  M(2,1) * N(1,3);
            MN(1,7)   = -M(1,2) * N(2,2);
            MN(2,8)   =  M(2,1) * N(1,3); 
            MN(1,9)   =  M(1,2) * N(2,3); 
            MN(2,10) =  M(2,2) * N(1,3); 
            MN(1,11) = -M(1,2) * N(2,4); 
            MN(2,12) =  M(2,2) * N(1,4); 
            MN(1,13) =  M(1,1) * N(2,3); 
            MN(2,14) =  M(2,2) * N(1,1); 
            MN(1,15) = -M(1,1) * N(2,4); 
            MN(2,16) =  M(2,2) * N(1,2); 
            
            
            Tr = sparse(16,12);
            
            Tr(1,1) = 1; 
            Tr(2,2) = 1; 
            Tr(3,3) = 
            
            
            
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(shellElement3d4n)
            EModul = reissnerMindlinElement3d4n.getPropertyValue('YOUNGS_MODULUS');
            prxy = reissnerMindlinElement3d4n.getPropertyValue('POISSON_RATIO');
            nr_gauss_points = reissnerMindlinElement3d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            thickness = reissnerMindlinElement3d4n.getPropertyValue('THICKNESS');
            alpha_shear = reissnerMindlinElement3d4n.getPropertyValue('SHEAR_CORRECTION_FACTOR');     % shear correction factor
        
        end
        
    end
end




