classdef EigensolverStrategy < Solver
    %EIGENSOLVERSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        femModel
        assembler
        isInitialized
        
        stiffnessMatrix
        dampingMatrix
        massMatrix
        
        eigenfrequencies
        modalMatrix
        spectralMatrix
        
        normalizedModalMatrix
        modalSuperpositionIsInitialized
        rayleighAlpha
        rayleighBeta
        modalDampingRatio
    end
    
    methods
        
        function eigensolver = EigensolverStrategy(femModel)
            if nargin > 0
                eigensolver.femModel = femModel;
                eigensolver.assembler = SimpleAssembler(femModel);
                eigensolver.isInitialized = false;
                eigensolver.modalSuperpositionIsInitialized = false;
            end
        end
        
        function solve(eigensolver)
            if ~ eigensolver.isInitialized
                eigensolver.initialize();
            end
            
            % solve
            [eigensolver.modalMatrix, eigensolver.spectralMatrix] = ...
                eig(eigensolver.stiffnessMatrix, eigensolver.massMatrix);
            
            % get eigenfrequencies
            eigensolver.eigenfrequencies = diag(eigensolver.spectralMatrix);
            eigensolver.eigenfrequencies = sqrt(eigensolver.eigenfrequencies) ./ (2*pi);
        end
        
        function assignModeShapes(eigensolver)
            nModes = size(eigensolver.modalMatrix, 1);
            for itEv = 1:nModes
                eigensolver.assembler.appendValuesToDofs(eigensolver.femModel, eigensolver.modalMatrix(:,itEv));
            end
        end
        
        function solveModalSuperposition(solver)
            if isempty(solver.eigenfrequencies)
                solver.solve();
            end
            
            if ~solver.modalSuperpositionIsInitialized
                solver.initializeModalSuperposition();
            end
            
            
            %modalMassMatrix = eye(nEigenvectors); %= normModalMatrix.' * solver.massMatrix * normModalMatrix
            %omega^2: diag(solver.spectralMatrix)
            
            %(y_n°°) + (2*zeta_n*omega_n*y_n°) + (omega_n^2*y_n) = p_n/M_n
            %zeta_n = alpha/2/omega_n + beta/2*omega_n
            %M_n is identity -> M_n = 1
            %->(y_n°°) + ((alpha+beta*omega_n^2)*y_n°) + (omega_n^2*y_n) = p_n
            
            
            
        end
        
        function solveUndampedModalSuperposition(solver)
            if isempty(solver.eigenfrequencies)
                solver.solve();
            end
            solver.normalizeModalMatrix();
            
            solver.femModel.getAllNodes.addNewValue(["VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z", ...
                "ACCELERATION_X", "ACCELERATION_Y", "ACCELERATION_Z"]);
            [~, fixedDofs] = solver.femModel.getDofConstraints();
            fixedDofIds = fixedDofs.getId();
            
            disp = solver.assembler.assemble3dDofVector(solver.femModel, 'DISPLACEMENT');
            disp0 = applyVectorBoundaryConditions(disp, fixedDofIds)';
            vel = solver.assembler.assemble3dValueVector(solver.femModel, 'VELOCITY');
            vel0 = applyVectorBoundaryConditions(vel, fixedDofIds)';
            [~, force0] = solver.assembler.applyExternalForces(solver.femModel);
            
            t = 0;
            endTime = 2000;
            dt = 3;
            
            nModes = length(solver.eigenfrequencies);
            
            mc = zeros(nModes,1);  %modal coordinate eta
            mcd = zeros(nModes,1); %derivative of modal coordinate eta
            pfo = zeros(nModes,1); %old modal participation factor phi (load)
            pfn = zeros(nModes,1); %new modal participation factor phi (load)
            for m = 1:nModes
                pfo(m) = solver.modalMatrix(:,m).' * solver.massMatrix * force0';
                mc(m) = solver.modalMatrix(:,m).' * solver.massMatrix * disp0;
                mcd(m) = solver.modalMatrix(:,m).' * solver.massMatrix * vel0;
            end
            
            while t < endTime
                for m = 1:nModes
                    eigenvalues = diag(solver.spectralMatrix);
                    alpha = sqrt(eigenvalues(m)) * dt;
                    [~, force] = solver.assembler.applyExternalForces(solver.femModel);
                    pfn(m) = solver.modalMatrix(:,m).' * solver.massMatrix * force';
                    res = [cos(alpha) sin(alpha)/alpha; -alpha*sin(alpha) cos(alpha)] ...
                        * [mc(m); dt*mcd(m)] ...
                        + dt^2/alpha^2 ...
                        * [(sin(alpha)/alpha-cos(alpha)) (1-sin(alpha)/alpha); (alpha*sin(alpha)+cos(alpha)-1) (1-cos(alpha))] ...
                        * [pfo(m); pfn(m)];
                    
                    %new modal coordinates
                    mc(m) = res(1);
                    mcd(m) = res(2) / dt;
                    
                    pfo(m) = pfn(m);
                end
                
                %transform back
                disp = 0;
                vel = 0;
                for m = 1:nModes
                    disp = disp + mc(m) * solver.modalMatrix(:,m);
                    vel = vel + mcd(m) * solver.modalMatrix(:,m);
                end
                solver.assembler.appendValuesToDofs(solver.femModel, disp);
                solver.assembler.appendValuesToNodes(solver.femModel, 'VELOCITY', vel);
                
                t = t + dt;
            end
            
        end
        
        function harmonicAnalysis(solver, excitations)
            if isempty(solver.eigenfrequencies)
                solver.solve();
            end
            
%             solver.femModel.getAllNodes.addNewValue(["VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z", ...
%                 "ACCELERATION_X", "ACCELERATION_Y", "ACCELERATION_Z"]);
            
            [~, force] = solver.assembler.applyExternalForces(solver.femModel);
            nModes = length(solver.eigenfrequencies);
            
            eigenvalues = diag(solver.spectralMatrix);
            
            for e=1:length(excitations)
                excitation = excitations(e);
                result = zeros;
                for n=1:nModes
                    if solver.modalDampingRatio ~= 0
                        dampingRatio = solver.modalDampingRatio;
                    elseif (solver.rayleighAlpha ~= 0) || (solver.rayleighBeta ~= 0)
                        dampingRatio = solver.rayleighAlpha / (2 * sqrt(eigenvalues(n))) + solver.rayleighBeta * sqrt(eigenvalues(n)) / 2;
                    else
                        dampingRatio = 0.0;
                    end
                    factor = (eigenvalues(n) - excitation^2) + 2i * dampingRatio * sqrt(eigenvalues(n)) * excitation;
                    result = result + 1/factor .* ((solver.modalMatrix(:,n) * solver.modalMatrix(:,n)') * force');
                end
                solver.assembler.appendValuesToDofs(solver.femModel, result);
%                 solver.assembler.appendValuesToNodes(solver.femModel, 'VELOCITY', (1i * excitation) .* result)
%                 solver.assembler.appendValuesToNodes(solver.femModel, 'ACCELERATION', (- power(excitation, 2)) .* result)
            end    
            solver.femModel.getDofArray.removeValue(1);
            
        end
        
        function normalizeModalMatrix(solver)
            nEigenvectors = size(solver.modalMatrix,1);
            solver.normalizedModalMatrix = zeros(nEigenvectors);
            for ii = 1:nEigenvectors
               solver.normalizedModalMatrix(:,ii) = ...
                   1 / sqrt(solver.modalMatrix(:,ii).' * solver.massMatrix * solver.modalMatrix(:,ii))...
                   .* solver.modalMatrix(:,ii);
            end
        end
        
%         function modalSuperposition(eigensolver)
%             generalizedMassMatrix = eigensolver.modalMatrix.' ...
%                 * eigensolver.massMatrix ...
%                 * eigensolver.modalMatrix;
%             init_disp = [1 0];
%             
% %             generalizedStiffnessMatrix = eigensolver.modalMatrix.' ...
% %                 * eigensolver.stiffnessMatrix ...
% %                 * eigensolver.modalMatrix;
% 
%             xmu = eigensolver.modalMatrix * eigensolver.massMatrix * init_disp';
%             init = 1./diag(generalizedMassMatrix) .* xmu;
% 
% %             init = 1/diag(generalizedMassMatrix) * eigensolver.modalMatrix * eigensolver.massMatrix * init_disp;
%         end
        
        function initialize(eigensolver)
            eigensolver.femModel.initialize;
            
            % assemble and reduce matrices
            [~, fixedDofs] = eigensolver.femModel.getDofConstraints();
            fixedDofIds = fixedDofs.getId();
            
            mass = eigensolver.assembler.assembleGlobalMassMatrix(eigensolver.femModel);
            eigensolver.massMatrix = applyMatrixBoundaryConditions(mass, fixedDofIds);
            damping = eigensolver.assembler.assembleGlobalDampingMatrix(eigensolver.femModel);
            eigensolver.dampingMatrix = applyMatrixBoundaryConditions(damping, fixedDofIds);
            stiffness = eigensolver.assembler.assembleGlobalStiffnessMatrix(eigensolver.femModel);
            eigensolver.stiffnessMatrix = applyMatrixBoundaryConditions(stiffness, fixedDofIds);
            
            %get the damping coefficients
            if (eigensolver.femModel.getElement(1).getMaterial.hasValue('RAYLEIGH_ALPHA')) ...
                    && (eigensolver.femModel.getElement(1).getMaterial.hasValue('RAYLEIGH_BETA'))
                eigensolver.rayleighAlpha = eigensolver.femModel.getElement(1).getMaterial.getValue('RAYLEIGH_ALPHA');
                eigensolver.rayleighBeta = eigensolver.femModel.getElement(1).getMaterial.getValue('RAYLEIGH_BETA');
            else
                eigensolver.rayleighAlpha = 0.0;
                eigensolver.rayleighBeta = 0.0;
            end
            
            if (eigensolver.femModel.getElement(1).getMaterial.hasValue('MODAL_DAMPING_RATIO'))
                eigensolver.modalDampingRatio = eigensolver.femModel.getElement(1).getMaterial.getValue('MODAL_DAMPING_RATIO');
            else
                eigensolver.modalDampingRatio = 0.0;
            end
            
            eigensolver.isInitialized = true;
        end
        
        function initializeModalSuperposition(solver)
            % normalize eigenvectors
            solver.normalizeModalMatrix();
            
            % initial acceleration FALSCH
%             solver.femModel.getAllNodes.addNewValue(["VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z", ...
%                 "ACCELERATION_X", "ACCELERATION_Y", "ACCELERATION_Z"]);
%             
%             [~, force0] = solver.assembler.applyExternalForces(solver.femModel);
%             disp = solver.assembler.assemble3dDofVector(solver.femModel, 'DISPLACEMENT');
%             disp0 = applyVectorBoundaryConditions(disp, fixedDofIds)';
%             vel = solver.assembler.assemble3dValueVector(solver.femModel, 'VELOCITY');
%             vel0 = applyVectorBoundaryConditions(vel, fixedDofIds)';
%             
%             acc0 = (solver.massMatrix) \ (force0' - solver.stiffnessMatrix * disp0 - solver.dampingMatrix * vel0);
%             solver.assembler.assignValuesToNodes(solver.femModel, 'ACCELERATION', acc0, 1);
            
            % get rayleigh alpha and beta
            %this assumes, that all elements have the same alpha and beta!
            if (solver.femModel.getElement(1).getMaterial.hasValue('RAYLEIGH_ALPHA')) ...
                    && (solver.femModel.getElement(1).getMaterial.hasValue('RAYLEIGH_BETA'))
                solver.rayleighAlpha = solver.femModel.getElement(1).getMaterial.getValue('RAYLEIGH_ALPHA');
                solver.rayleighBeta = solver.femModel.getElement(1).getMaterial.getValue('RAYLEIGH_BETA');
            else
                solver.rayleighAlpha = 0.0;
                solver.rayleighBeta = 0.0;
            end
            
            % set flag to true
            solver.modalSuperpositionIsInitialized = true;
        end
        
        function ef = getEigenfrequencies(eigensolver)
            ef = eigensolver.eigenfrequencies;
        end
        
        function mm = getModalMatrix(eigensolver)
            mm = eigensolver.modalMatrix;
        end
        
    end
    
end

