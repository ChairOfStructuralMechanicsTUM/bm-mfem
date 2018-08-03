classdef SimpleHarmonicSolvingStrategy < Solver
    %SIMPLESOLVINGSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        femModel
        omega
        isInitialized
    end
        
    methods
        
        %constructor
        function simpleHarmonicSolver = SimpleHarmonicSolvingStrategy(femModel,omega)
            if (nargin > 0)
               simpleHarmonicSolver.femModel = femModel;
               simpleHarmonicSolver.omega = omega;
               simpleHarmonicSolver.isInitialized = false;
            end
        end
        
        function x = solve(simpleHarmonicSolver)
            if ~ simpleHarmonicSolver.isInitialized
                simpleHarmonicSolver.initialize();
            end
            
            [~, Kred] = SimpleAssembler.assembleGlobalStiffnessMatrix(simpleHarmonicSolver.femModel);
            [~, Mred] = SimpleAssembler.assembleGlobalMassMatrix(simpleHarmonicSolver.femModel);
            [~, fred] = SimpleAssembler.applyExternalForces(simpleHarmonicSolver.femModel);
            
            Kdyn = Kred - simpleHarmonicSolver.omega^2 * Mred;
            x = linsolve(Kdyn, fred.');
            
            SimpleAssembler.assignResultsToDofs(simpleHarmonicSolver.femModel, x);
        end
        
        function initialize(simpleHarmonicSolver)
            simpleHarmonicSolver.femModel.initialize;
            simpleHarmonicSolver.isInitialized = true;
        end
        
        function nodalForces = getNodalForces(simpleHarmonicSolver, step)
            nodalForces = SimpleAssembler.applyExternalForces(simpleHarmonicSolver.femModel);
            [~, fixedDofs] = simpleHarmonicSolver.femModel.getDofConstraints();
            
            K = SimpleAssembler.assembleGlobalStiffnessMatrix(simpleHarmonicSolver.femModel);
            nodalForces(fixedDofs.getId) = K(fixedDofs.getId, :) ...
                * simpleHarmonicSolver.femModel.getDofArray.getValue(step);
            
        end
        
    end
end

