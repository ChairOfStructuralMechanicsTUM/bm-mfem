classdef ShellElement3d4n < QuadrilateralElement 
    %SHELLELEMENT  Summary of this class goes here
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
        
        
        
        
    end
end
