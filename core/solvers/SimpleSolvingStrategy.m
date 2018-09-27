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
            
            [~, Kred] = SimpleAssembler.assembleGlobalStiffnessMatrix(simpleSolver.femModel);
            [~, fred] = SimpleAssembler.applyExternalForces(simpleSolver.femModel);
            
%             x = linsolve(Kred, fred.');
            x = Kred\fred.';
            
            SimpleAssembler.assignResultsToDofs(simpleSolver.femModel, x);
        end
        
        function initialize(simpleSolver)
            simpleSolver.femModel.initialize;
            simpleSolver.isInitialized = true;
        end
        
        function nodalForces = getNodalForces(simpleSolver, step)
            nodalForces = SimpleAssembler.applyExternalForces(simpleSolver.femModel);
            [~, fixedDofs] = simpleSolver.femModel.getDofConstraints();
            
            K = SimpleAssembler.assembleGlobalStiffnessMatrix(simpleSolver.femModel);
            nodalForces(fixedDofs.getId) = K(fixedDofs.getId, :) ...
                * simpleSolver.femModel.getDofArray.getValue(step);
            
        end
        
    end
end

