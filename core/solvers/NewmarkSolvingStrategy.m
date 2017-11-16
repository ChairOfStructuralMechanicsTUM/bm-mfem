classdef NewmarkSolvingStrategy < Solver
    %NEWMARKSOLVINGSTRATEGY A time integration method based on the Newmark
    %method
    %   Newmark's method as outlined in Géradin, Michel, and Daniel J. 
    %   Rixen. Mechanical vibrations: theory and application to structural
    %   dynamics. John Wiley & Sons, 2014.
    
    properties (Access = private)
        femModel
        assembler
        isInitialized
        
        dt
        alpha
        beta
        gamma
        newmarkCoefficients
        
        massMatrix
        dampingMatrix
        stiffnessMatrix
        
        lhs
    end
    
    methods
        
        %constructor
        function solver = NewmarkSolvingStrategy(femModel, dt, varargin)
            solver.femModel = femModel;
            solver.dt = dt;
            solver.assembler = SimpleAssembler(femModel);
            solver.isInitialized = false;
            
            numvarargs = length(varargin);
            optargs = { 0.25 0.5 0 };
            optargs(1:numvarargs) = varargin;
            [solver.beta, solver.gamma, solver.alpha] = optargs{:};
        end
        
        function solve(solver)
            if ~solver.isInitialized
                solver.initialize;
            end
            [~, fixedDofs] = solver.femModel.getDofConstraints();
            fixedDofIds = fixedDofs.getId();
            
            [~, force] = solver.assembler.applyExternalForces(solver.femModel);
            dispOld = solver.assembler.assemble3dDofVector(solver.femModel, 'DISPLACEMENT');
            dispOld = applyVectorBoundaryConditions(dispOld, fixedDofIds)';
            velOld = solver.assembler.assemble3dValueVector(solver.femModel, 'VELOCITY');
            velOld = applyVectorBoundaryConditions(velOld, fixedDofIds)';
            accOld = solver.assembler.assemble3dValueVector(solver.femModel, 'ACCELERATION');
            accOld = applyVectorBoundaryConditions(accOld, fixedDofIds)';
            
            if solver.alpha == 0
                rhs = force' ...
                    - solver.dampingMatrix * (velOld + ((1 - solver.gamma) .* solver.dt) * accOld) ...
                    - solver.stiffnessMatrix * (dispOld + solver.dt .* velOld + ((0.5 - solver.beta) * power(solver.dt,2)) .* accOld);
                
                accNew = linsolve(solver.lhs, rhs);
                velNew = velOld + ((1 - solver.gamma) * solver.dt) .* accOld + (solver.gamma * solver.dt) .* accNew;
                dispNew = dispOld + solver.dt .* velOld + (power(solver.dt, 2) * (0.5 - solver.beta)) .* accOld + (power(solver.dt, 2) * solver.beta) .* accNew;
                
                
            else
                c = solver.newmarkCoefficients;
                rhs = force' ...
                    + solver.massMatrix * ((1 - solver.alpha) .* (dispOld .* c(1) + velOld .* c(3) + accOld .* c(4)) - solver.alpha .* accOld) ...
                    + solver.dampingMatrix * (dispOld .* c(2) + velOld .* c(5) + accOld .* c(6));
                
                dispNew = linsolve(solver.lhs, rhs);
                accNew = dispNew .* c(1) - dispOld .* c(1) - velOld .* c(3) - accOld .* c(4);
                velNew = velOld + accOld .* (solver.dt * (1 - solver.gamma)) + accNew .* (solver.dt * solver.gamma);
            end
            
            %write the values back to the dofs / nodes
            solver.assembler.appendValuesToDofs(solver.femModel, dispNew);
            solver.assembler.appendValuesToNodes(solver.femModel, 'VELOCITY', velNew);
            solver.assembler.appendValuesToNodes(solver.femModel, 'ACCELERATION', accNew);
            
        end
        
        function initialize(solver)
            solver.femModel.initialize;
            solver.femModel.getAllNodes.addNewValue(["VELOCITY", "ACCELERATION"]);
            
            % assemble and reduce matrices
            [~, fixedDofs] = solver.femModel.getDofConstraints();
            fixedDofIds = fixedDofs.getId();
            
            mass = solver.assembler.assembleGlobalMassMatrix(solver.femModel);
            solver.massMatrix = applyMatrixBoundaryConditions(mass, fixedDofIds);
            damping = solver.assembler.assembleGlobalDampingMatrix(solver.femModel);
            solver.dampingMatrix = applyMatrixBoundaryConditions(damping, fixedDofIds);
            stiffness = solver.assembler.assembleGlobalStiffnessMatrix(solver.femModel);
            solver.stiffnessMatrix = applyMatrixBoundaryConditions(stiffness, fixedDofIds);
            
            % initial acceleration
            [~, force0] = solver.assembler.applyExternalForces(solver.femModel);
            disp = solver.assembler.assemble3dDofVector(solver.femModel, 'DISPLACEMENT');
            disp0 = applyVectorBoundaryConditions(disp, fixedDofIds)';
            vel = solver.assembler.assemble3dValueVector(solver.femModel, 'VELOCITY');
            vel0 = applyVectorBoundaryConditions(vel, fixedDofIds)';
            
            acc0 = (solver.massMatrix) \ (force0' - solver.stiffnessMatrix * disp0 - solver.dampingMatrix * vel0);
            solver.assembler.assignValuesToNodes(solver.femModel, 'ACCELERATION', acc0, 1);
            solver.isInitialized = true;
            
            if solver.alpha == 0
                %set up left hand side of the les
                solver.lhs = solver.massMatrix ...
                    + (solver.gamma * solver.dt) .* solver.dampingMatrix ...
                    + (solver.beta * power(solver.dt,2)) .* solver.stiffnessMatrix;
            else
                %coefficients for bossak/newmark
                solver.beta = power((1.0 - solver.alpha),2) * solver.beta;
                solver.gamma = solver.gamma - solver.alpha;
                c = zeros(6,1);
                c(1) = 1 / (solver.beta * power(solver.dt,2));
                c(2) = solver.gamma / (solver.beta * solver.dt);
                c(3) = 1 / (solver.beta * solver.dt);
                c(4) = 1 / (2 * solver.beta) - 1;
                c(5) = solver.gamma / solver.beta - 1;
                c(6) = solver.dt * 0.5 * (solver.gamma / solver.beta - 2);
                solver.newmarkCoefficients = c;
                
                solver.lhs = solver.stiffnessMatrix ...
                    + ((1 - solver.alpha) * c(1)) .* solver.massMatrix ...
                    + c(2) .* solver.dampingMatrix;
            end
        end
    end
    
end

