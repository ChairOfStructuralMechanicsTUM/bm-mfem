classdef MORStrategy < Solver
    %MORSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        femModel
        assembler
        isInitialized
        
        %stiffnessMatrix
        %dampingMatrix
        %massMatrix
        %basis
        
        stiffnessMatrixR
        massMatrixR
        basisR
        forceR
        
        MORmethod = 'krylov-galerkin-proj'
    end
    
    methods
        
        function strategy = MORStrategy(femModel)
            if nargin > 0
                strategy.femModel = femModel;
                strategy.assembler = SimpleAssembler(femModel);
                strategy.isInitialized = false;
            else
                error("Error (EigensolverStratey): no fem model defined!")
            end
        end
        
        function solve(s, samplingPoints)
            if ~ s.isInitialized
                error('MOR strategy not initialized!')
            end
            
            for i = samplingPoints
                disp(i)
                uR = (-i^2*s.massMatrixR + s.stiffnessMatrixR) \ s.forceR;
                
                s.assembler.appendValuesToDofs(s.femModel, uR);
            end
                
        end
        
        function initialize(s, sampling)
            if ~ s.femModel.isInitialized()
                s.femModel.initialize;
            end
            s.femModel.initialize;
            
            % assemble and reduce matrices
            massMatrix = s.assembler.assembleGlobalMassMatrix(s.femModel);
            %s.dampingMatrix = s.assembler.assembleGlobalDampingMatrix(s.femModel);
            stiffnessMatrix = s.assembler.assembleGlobalStiffnessMatrix(s.femModel);
            
            [~, fixedDofs] = s.femModel.getDofConstraints();
            if ~ isempty(fixedDofs)
                fixedDofIds = fixedDofs.getId();
                massMatrix = applyMatrixBoundaryConditions(massMatrix, fixedDofIds);
                %s.dampingMatrix = applyMatrixBoundaryConditions(s.dampingMatrix, fixedDofIds);
                stiffnessMatrix = applyMatrixBoundaryConditions(stiffnessMatrix, fixedDofIds);                
            end
            [~, force] = s.assembler.applyExternalForces(s.femModel);
            
            switch s.MORmethod
                case 'krylov-galerkin-proj'
                    
                    K_dyn = stiffnessMatrix - (sampling)^2 * massMatrix;
                    b = K_dyn \ force;
                    Ab = K_dyn \ (massMatrix * b);
                    AAb = K_dyn\ (massMatrix * Ab);
                    basis = [b,Ab,AAb];
                    
                    [s.basisR, ~] = qr(basis, 0);
                    s.massMatrixR = s.basisR.' * s.massMatrix * s.basisR;
                    s.stiffnessMatrixR = s.basisR.' * s.stiffnessMatrix * s.basisR;
                    s.forceR = s.basisR.' * force;
                    
                otherwise
                    error('unknown MOR method!')
            end
            
            s.isInitialized = true;
        end
        
    end
    
end

