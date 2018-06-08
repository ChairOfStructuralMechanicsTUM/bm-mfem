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
        
        normalizedModalMatrix
        rayleighAlpha
        rayleighBeta
        modalDampingRatio
    end
    
    methods
        
        function obj = EigensolverStrategy(femModel)
            if nargin > 0
                obj.femModel = femModel;
                obj.assembler = SimpleAssembler(femModel);
                obj.isInitialized = false;
            else
                error("Error (EigensolverStratey): no fem model defined!")
            end
            
        end
        
        function solve(obj, nModes)
            %EIGENSOLVERSTRATEGY.SOLVE computes the eigenfrequencies and mode shapes of the
            %given system
            if ~ obj.isInitialized
                obj.initialize();
            end
            
            if size(obj.stiffnessMatrix,1) < nModes
                nModes = size(obj.stiffnessMatrix,1);
                warning('EigensolverStrategy: setting number of modes to %d\n',nModes);
            end
            
            % solve
            [obj.modalMatrix, spectralMatrix] = ...
                eigs(obj.stiffnessMatrix, obj.massMatrix, ...
                nModes, 'sm');
            
            % save eigenfrequencies in rad/s
            obj.eigenfrequencies = sqrt(diag(spectralMatrix));
            
            % sort eigenfrequencies and mode shapes
            [obj.eigenfrequencies, order] = sort(obj.eigenfrequencies);
            obj.modalMatrix = obj.modalMatrix(:,order);
        end
        
        function assignModeShapes(obj)
            %EIGENSOLVERSTRATEGY.ASSIGNMODESHAPES assigns the modal displacement to all dofs.
            %   step 1 visualizes the first eigenvector etc.
            nModes = size(obj.modalMatrix, 2);
            obj.assembler.assignResultsToDofs(obj.femModel, obj.modalMatrix(:,1));
            for itEv = 2:nModes
                obj.assembler.appendValuesToDofs(obj.femModel, obj.modalMatrix(:,itEv));
            end
        end
        
        function solveUndampedModalSuperposition(obj)
            %SOLVEUNDAMPEDMODALSUPERPOSITION using a time integration
            %scheme for arbitrary excitations
            
            if isempty(obj.eigenfrequencies)
                obj.solve();
            end
            obj.normalizeModalMatrix();
            
            obj.femModel.getAllNodes.addNewValue(["VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z", ...
                "ACCELERATION_X", "ACCELERATION_Y", "ACCELERATION_Z"]);
            
            disp = obj.assembler.assemble3dDofVector(obj.femModel, 'DISPLACEMENT');
            vel = obj.assembler.assemble3dValueVector(obj.femModel, 'VELOCITY');
            [~, force0] = obj.assembler.applyExternalForces(obj.femModel);
            
            [~, fixedDofs] = obj.femModel.getDofConstraints();
            if ~ isempty(fixedDofs)
                fixedDofIds = fixedDofs.getId();
                vel0 = applyVectorBoundaryConditions(vel, fixedDofIds)';
                disp0 = applyVectorBoundaryConditions(disp, fixedDofIds)';
            end
            
            t = 0;
            endTime = 2000;
            dt = 3;
            
            nModes = length(obj.eigenfrequencies);
            
            mc = zeros(nModes,1);  %modal coordinate eta
            mcd = zeros(nModes,1); %derivative of modal coordinate eta
            pfo = zeros(nModes,1); %old modal participation factor phi (load)
            pfn = zeros(nModes,1); %new modal participation factor phi (load)
            for m = 1:nModes
                pfo(m) = obj.modalMatrix(:,m).' * obj.massMatrix * force0';
                mc(m) = obj.modalMatrix(:,m).' * obj.massMatrix * disp0;
                mcd(m) = obj.modalMatrix(:,m).' * obj.massMatrix * vel0;
            end
            
            while t < endTime
                for m = 1:nModes
                    eigenvalues = diag(obj.spectralMatrix);
                    alpha = sqrt(eigenvalues(m)) * dt;
                    [~, force] = obj.assembler.applyExternalForces(obj.femModel);
                    pfn(m) = obj.modalMatrix(:,m).' * obj.massMatrix * force';
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
                    disp = disp + mc(m) * obj.modalMatrix(:,m);
                    vel = vel + mcd(m) * obj.modalMatrix(:,m);
                end
                obj.assembler.appendValuesToDofs(obj.femModel, disp);
                obj.assembler.appendValuesToNodes(obj.femModel, 'VELOCITY', vel);
                
                t = t + dt;
            end
            
        end
        
        function harmonicAnalysis(obj, excitations, nModes)
            %HARMONIC ANALYSIS performs an analysis in the frequency domain
            %using a harmonic excitation
            %EXCITATIONS: vector of all excitation frequencies in rad/s
            %NMODES: number of modes
            if isempty(obj.eigenfrequencies)
                obj.solve(nModes);
            end
            
%             solver.femModel.getAllNodes.addNewValue(["VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z", ...
%                 "ACCELERATION_X", "ACCELERATION_Y", "ACCELERATION_Z"]);
            
            [~, force] = obj.assembler.applyExternalForces(obj.femModel);
            nModes = length(obj.eigenfrequencies);
            
            for e=1:length(excitations)
                excitation = excitations(e);
                result = zeros;
                for n=1:nModes
                    if obj.modalDampingRatio ~= 0
                        dampingRatio = obj.modalDampingRatio;
                    elseif (obj.rayleighAlpha ~= 0) || (obj.rayleighBeta ~= 0)
                        dampingRatio = obj.rayleighAlpha / (2 * obj.eigenfrequencies(n) ) ...
                            + obj.rayleighBeta * obj.eigenfrequencies(n) / 2;
                    else
                        dampingRatio = 0.0;
                    end
                    factor = (obj.eigenfrequencies(n)^2 - excitation^2) + 2i * dampingRatio * obj.eigenfrequencies(n) * excitation;
                    result = result + 1/factor .* ((obj.modalMatrix(:,n) * obj.modalMatrix(:,n)') * force');
                end
                obj.assembler.appendValuesToDofs(obj.femModel, result);
%                 solver.assembler.appendValuesToNodes(solver.femModel, 'VELOCITY', (1i * excitation) .* result)
%                 solver.assembler.appendValuesToNodes(solver.femModel, 'ACCELERATION', (- power(excitation, 2)) .* result)
            end    
            obj.femModel.getDofArray.removeValue(1);
            
        end
        
        function normalizeModalMatrix(obj)
            nEigenvectors = size(obj.modalMatrix,1);
            obj.normalizedModalMatrix = zeros(nEigenvectors);
            for ii = 1:nEigenvectors
               obj.normalizedModalMatrix(:,ii) = ...
                   1 / sqrt(obj.modalMatrix(:,ii).' * obj.massMatrix * obj.modalMatrix(:,ii))...
                   .* obj.modalMatrix(:,ii);
            end
        end
        
        function initialize(obj)
            if ~ obj.femModel.isInitialized()
                obj.femModel.initialize;
            end
            
            % assemble and reduce matrices
            obj.massMatrix = obj.assembler.assembleGlobalMassMatrix(obj.femModel);
            obj.dampingMatrix = obj.assembler.assembleGlobalDampingMatrix(obj.femModel);
            obj.stiffnessMatrix = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);
            
            [~, fixedDofs] = obj.femModel.getDofConstraints();
            if ~ isempty(fixedDofs)
                fixedDofIds = fixedDofs.getId();
                obj.massMatrix = applyMatrixBoundaryConditions(obj.massMatrix, fixedDofIds);
                obj.dampingMatrix = applyMatrixBoundaryConditions(obj.dampingMatrix, fixedDofIds);
                obj.stiffnessMatrix = applyMatrixBoundaryConditions(obj.stiffnessMatrix, fixedDofIds);                
            end
            
            %get the damping coefficients
            if (obj.femModel.getElement(1).getProperties.hasValue('RAYLEIGH_ALPHA')) ...
                    && (obj.femModel.getElement(1).getProperties.hasValue('RAYLEIGH_BETA'))
                obj.rayleighAlpha = obj.femModel.getElement(1).getProperties.getValue('RAYLEIGH_ALPHA');
                obj.rayleighBeta = obj.femModel.getElement(1).getProperties.getValue('RAYLEIGH_BETA');
            else
                obj.rayleighAlpha = 0.0;
                obj.rayleighBeta = 0.0;
            end
            
            if (obj.femModel.getElement(1).getProperties.hasValue('MODAL_DAMPING_RATIO'))
                obj.modalDampingRatio = obj.femModel.getElement(1).getProperties.getValue('MODAL_DAMPING_RATIO');
            else
                obj.modalDampingRatio = 0.0;
            end
            
            obj.isInitialized = true;
        end
        
        function ef = getEigenfrequencies(obj)
            %GETEIGENFREQUENCIES returns the eigenfrequencies in rad/s
            %sorted from low to high
            ef = obj.eigenfrequencies;
        end
        
        function mm = getModalMatrix(obj)
            %GETMODALMATRIX returns the modal matrix of the system sorted
            %according to the eigenfrequencies
            mm = obj.modalMatrix;
        end
        
    end
    
end

