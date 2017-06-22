classdef SimpleSolvingStrategy < Solver
    %SIMPLESOLVINGSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        femModel
        isInitialized
    end
    
    
    
    methods
        
        %constructor
        function simpleSolver = SimpleSolvingStrategy(femModel)
            if (nargin > 0)
               simpleSolver.femModel = femModel;
               simpleSolver.isInitialized = false;
            end
        end
        
        function x = solve(simpleSolver)
            if ~ simpleSolver.isInitialized
                simpleSolver.initialize();
            end
            
            Kred = SimpleAssembler(simpleSolver.femModel).reducedStiffnessMatrix;
            f = SimpleAssembler(simpleSolver.femModel).reducedForceVector;
            
            x = linsolve(Kred, f.');
            
            SimpleAssembler.assignResultsToDofs(simpleSolver.femModel, x);
        end
        
        function initialize(simpleSolver)
            elements = simpleSolver.femModel.getAllElements();
            for ii = 1:length(elements)
                elements(ii).check();
            end
        end
        
        function nodalForces = getNodalForces(simpleSolver)
            
            dofs = simpleSolver.femModel.getDofArray;
            nDofs = length(dofs);
            nodalForces = SimpleAssembler(simpleSolver.femModel).forceVector;
            fixedDofs = [];
            
            for itDof = 1:nDofs
                if (dofs(itDof).isFixed)
                    fixedDofs = [fixedDofs itDof];                          % array of fixed dofs and their location
                else
                    nodalForces(itDof) = dofs(itDof).getDofLoad;
                end
            end
            
            nodalForces(fixedDofs) = SimpleAssembler(simpleSolver.femModel).stiffnessMatrix(fixedDofs, :) ...
                * getValue(simpleSolver.femModel.getDofArray)';
            
        end
        
    end
end

