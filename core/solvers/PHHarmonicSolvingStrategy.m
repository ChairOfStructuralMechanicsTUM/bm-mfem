classdef PHHarmonicSolvingStrategy < Solver
    %PHHSOLVINGSTRATEGY works as SimpleHarmonicSolvingStrategy, but is specified
    %for the application on mixed porous and homogen models. Hence, it
    %calls the specified assembler Assembler_Porous_Homogen
    
    properties (Access = private)
        femModel
        omega
        isInitialized
    end
        
    methods
        
        %constructor
        function PHHarmonicSolver = PHHarmonicSolvingStrategy(femModel,omega)
            if (nargin > 0)
               PHHarmonicSolver.femModel = femModel;
               PHHarmonicSolver.omega = omega;
               PHHarmonicSolver.isInitialized = false;
            end
        end
        
        function x = solve(PHHarmonicSolver)
            if ~ PHHarmonicSolver.isInitialized
                PHHarmonicSolver.initialize();
            end
            
            [~, Kred] = Assembler_Porous_Homogen.assembleGlobalStiffnessMatrix(PHHarmonicSolver.femModel);
            [~, Mred] = Assembler_Porous_Homogen.assembleGlobalMassMatrix(PHHarmonicSolver.femModel);
            [~, fred,delete] = Assembler_Porous_Homogen.applyExternalForces(PHHarmonicSolver.femModel);
            
            PHHarmonicSolver.femModel.deleteDof(delete);
            
            Kdyn = Kred - PHHarmonicSolver.omega^2 * Mred;
            x = linsolve(Kdyn, fred.');
            
            Assembler_Porous_Homogen.assignResultsToDofs(PHHarmonicSolver.femModel, x);
        end
        
        function initialize(PHHarmonicSolver)
            PHHarmonicSolver.femModel.initialize;
            PHHarmonicSolver.isInitialized = true;
        end
        
        function nodalForces = getNodalForces(PHHarmonicSolver, step)
            nodalForces = Assembler_Porous_Homogen.applyExternalForces(PHHarmonicSolver.femModel);
            [~, fixedDofs] = PHHarmonicSolver.femModel.getDofConstraints();
            
            K = Assembler_Porous_Homogen.assembleGlobalStiffnessMatrix(PHHarmonicSolver.femModel);
            nodalForces(fixedDofs.getId) = K(fixedDofs.getId, :) ...
                * PHHarmonicSolver.femModel.getDofArray.getValue(step);
            
        end
        
    end
end

