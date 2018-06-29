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
        
        function x = solve(simpleHarmoniSolver)
            if ~ simpleHarmoniSolver.isInitialized
                simpleHarmoniSolver.initialize();
            end
            
            [~, Kred] = SimpleAssembler.assembleGlobalStiffnessMatrix(simpleHarmoniSolver.femModel);
            [~, Mred] = SimpleAssembler.assembleGlobalMassMatrix(simpleHarmoniSolver.femModel);
            [~, fred] = SimpleAssembler.applyExternalForces(simpleHarmoniSolver.femModel);
            
            Kdyn = Kred - simpleHarmoniSolver.omega^2 * Mred;
            x = linsolve(Kdyn, fred.');
            
            SimpleAssembler.assignResultsToDofs(simpleHarmoniSolver.femModel, x);
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

