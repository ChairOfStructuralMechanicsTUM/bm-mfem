classdef EigensolverStrategy < Solver
    %EIGENSOLVERSTRATEGY Solving strategy to obtain an eigenvalue analysis
    %   Contains also modal superposition and harmonic analysis
    
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
            else
                error("Error (EigensolverStratey): no fem model defined!")
            end
            
        end
        
        function solve(eigensolver, nModes)
            if ~ eigensolver.isInitialized
                eigensolver.initialize();
            end
            
            if size(eigensolver.stiffnessMatrix,1) < nModes
                nModes = size(eigensolver.stiffnessMatrix,1);
                fprintf("Warning (EigensolverStrategy): setting number of modes to %d\n",nModes);
            end
            
            % solve
            [eigensolver.modalMatrix, eigensolver.spectralMatrix] = ...
                eigs(eigensolver.stiffnessMatrix, eigensolver.massMatrix, ...
                nModes, 'sm');
            
            % get eigenfrequencies
            eigensolver.eigenfrequencies = diag(eigensolver.spectralMatrix);
            eigensolver.eigenfrequencies = sqrt(eigensolver.eigenfrequencies) ./ (2*pi);
        end
        
        function assignModeShapes(eigensolver)
            nModes = size(eigensolver.modalMatrix, 2);
            for itEv = 1:nModes
                eigensolver.assembler.appendValuesToDofs(eigensolver.femModel, eigensolver.modalMatrix(:,itEv));
            end
        end
        
        function solveUndampedModalSuperposition(solver)
            %SOLVEUNDAMPEDMODALSUPERPOSITION using a time integration
            %scheme for arbitrary excitations
            
            if isempty(solver.eigenfrequencies)
                solver.solve();
            end
            solver.normalizeModalMatrix();
            
            solver.femModel.getAllNodes.addNewValue(["VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z", ...
                "ACCELERATION_X", "ACCELERATION_Y", "ACCELERATION_Z"]);
            
            disp = solver.assembler.assemble3dDofVector(solver.femModel, 'DISPLACEMENT');
            vel = solver.assembler.assemble3dValueVector(solver.femModel, 'VELOCITY');
            [~, force0] = solver.assembler.applyExternalForces(solver.femModel);
            
            [~, fixedDofs] = solver.femModel.getDofConstraints();
            if ~ isempty(fixedDofs)
                fixedDofIds = fixedDofs.getId();
                vel0 = applyVectorBoundaryConditions(vel, fixedDofIds)';
                disp0 = applyVectorBoundaryConditions(disp, fixedDofIds)';
            end
            
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
        
        function harmonicAnalysis(solver, excitations, nModes)
            %HARMONIC ANALYSIS performs an analysis in the frequency domain
            %using a harmonic excitation
            %EXCITATIONS: vector of all excitation frequencies in rad/s
            %NMODES: number of modes
            if isempty(solver.eigenfrequencies)
                solver.solve(nModes);
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
        
        function initialize(eigensolver)
            if ~ eigensolver.femModel.isInitialized()
                eigensolver.femModel.initialize;
            end
            
            % assemble and reduce matrices
            eigensolver.massMatrix = eigensolver.assembler.assembleGlobalMassMatrix(eigensolver.femModel);
            eigensolver.dampingMatrix = eigensolver.assembler.assembleGlobalDampingMatrix(eigensolver.femModel);
            eigensolver.stiffnessMatrix = eigensolver.assembler.assembleGlobalStiffnessMatrix(eigensolver.femModel);
            
            [~, fixedDofs] = eigensolver.femModel.getDofConstraints();
            if ~ isempty(fixedDofs)
                fixedDofIds = fixedDofs.getId();
                eigensolver.massMatrix = applyMatrixBoundaryConditions(eigensolver.massMatrix, fixedDofIds);
                eigensolver.dampingMatrix = applyMatrixBoundaryConditions(eigensolver.dampingMatrix, fixedDofIds);
                eigensolver.stiffnessMatrix = applyMatrixBoundaryConditions(eigensolver.stiffnessMatrix, fixedDofIds);                
            end
            
            %get the damping coefficients
            if (eigensolver.femModel.getElement(1).getProperties.hasValue('RAYLEIGH_ALPHA')) ...
                    && (eigensolver.femModel.getElement(1).getProperties.hasValue('RAYLEIGH_BETA'))
                eigensolver.rayleighAlpha = eigensolver.femModel.getElement(1).getProperties.getValue('RAYLEIGH_ALPHA');
                eigensolver.rayleighBeta = eigensolver.femModel.getElement(1).getProperties.getValue('RAYLEIGH_BETA');
            else
                eigensolver.rayleighAlpha = 0.0;
                eigensolver.rayleighBeta = 0.0;
            end
            
            if (eigensolver.femModel.getElement(1).getProperties.hasValue('MODAL_DAMPING_RATIO'))
                eigensolver.modalDampingRatio = eigensolver.femModel.getElement(1).getProperties.getValue('MODAL_DAMPING_RATIO');
            else
                eigensolver.modalDampingRatio = 0.0;
            end
            
            eigensolver.isInitialized = true;
        end
        
        function ef = getEigenfrequencies(eigensolver)
            ef = eigensolver.eigenfrequencies;
        end
        
        function mm = getModalMatrix(eigensolver)
            mm = eigensolver.modalMatrix;
        end
        
    end
    
end

