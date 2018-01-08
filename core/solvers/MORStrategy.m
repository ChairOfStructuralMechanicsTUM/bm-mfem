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
        
        MORmethod
    end
    
    methods
        
        function strategy = MORStrategy(femModel, method)
            if nargin > 0
                strategy.femModel = femModel;
                strategy.assembler = SimpleAssembler(femModel);
                strategy.MORmethod = method;
                strategy.isInitialized = false;
            else
                error("Error (MORStratey): no fem model defined!")
            end
        end
        
        function solve(s, samplingPoints)
            if ~ s.isInitialized
                error('MOR strategy not initialized!')
            end
            
            for ii = samplingPoints
                uR = (-ii^2*s.massMatrixR + s.stiffnessMatrixR) \ s.forceR;
                u = s.basisR * uR;
                s.assembler.appendValuesToDofs(s.femModel, u);
            end
                
        end
        
        function initialize(s, sampling)
            if ~ s.femModel.isInitialized()
                s.femModel.initialize;
            end
            
            % assemble and reduce matrices
            massMatrix = s.assembler.assembleGlobalMassMatrix(s.femModel);
            stiffnessMatrix = s.assembler.assembleGlobalStiffnessMatrix(s.femModel);
            
            [~, fixedDofs] = s.femModel.getDofConstraints();
            if ~ isempty(fixedDofs)
                fixedDofIds = fixedDofs.getId();
                massMatrix = applyMatrixBoundaryConditions(massMatrix, fixedDofIds);
                stiffnessMatrix = applyMatrixBoundaryConditions(stiffnessMatrix, fixedDofIds);                
            end
            [~, force] = s.assembler.applyExternalForces(s.femModel);
            
            switch s.MORmethod
                case 'krylov-galerkin-proj'
                    for ii = 1:length(sampling)
                        K_dyn = stiffnessMatrix - (sampling(ii))^2 * massMatrix;
                        b = K_dyn \ force';
                        Ab = K_dyn \ (massMatrix * b);
                        AAb = K_dyn\ (massMatrix * Ab);
                        basis(:,(ii-1)*3+1:(ii-1)*3+3) = [b,Ab,AAb];
                    end
                    [s.basisR, ~] = qr(basis, 0);
                    s.massMatrixR = s.basisR.' * massMatrix * s.basisR;
                    s.stiffnessMatrixR = s.basisR.' * stiffnessMatrix * s.basisR;
                    s.forceR = s.basisR.' * force';
                    
                case 'derivative-based-galerkin-proj'
                    ddS = -2 * massMatrix; %system matrix second derivative
                    for ii = 1:length(sampling)
                        sp = sampling(ii);
                        S = -sp^2 * massMatrix + stiffnessMatrix; %system matrix
                        dS = -2 * sp * massMatrix; %system matrix first derivative
                        u = S \ force';
                        dudw = S \ (-dS * u);
                        dduddw = S \ (-2*dS*dudw - ddS*u);
                        basis(:,(ii-1)*3+1:(ii-1)*3+3) = [u,dudw,dduddw];
                    end
                    [s.basisR, ~] = qr(basis, 0);
                    s.massMatrixR = s.basisR.' * massMatrix * s.basisR;
                    s.stiffnessMatrixR = s.basisR.' * stiffnessMatrix * s.basisR;
                    s.forceR = s.basisR.' * force';
                    
                otherwise
                    error('unknown MOR method!')
            end
            
            s.isInitialized = true;
        end
        
        function setMethod(s, method)
            s.method = method;
        end
        
    end
    
end

