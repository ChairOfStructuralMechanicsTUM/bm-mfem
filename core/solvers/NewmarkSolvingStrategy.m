classdef NewmarkSolvingStrategy < Solver
    %NEWMARKSOLVINGSTRATEGY A time integration method based on the Newmark
    %method
    %   Newmark's method as outlined in Géradin, Michel, and Daniel J. 
    %   Rixen. Mechanical vibrations: theory and application to structural
    %   dynamics. John Wiley & Sons, 2014.
    
    properties (Access = private)
        dt
        beta
        gamma
        femModel
        isInitialized = false
        assembler
        
        massMatrix
        dampingMatrix
        stiffnessMatrix
    end
    
    methods
        
        %constructor
        function solver = NewmarkSolvingStrategy(femModel, dt, varargin)
            solver.femModel = femModel;
            solver.dt = dt;
            solver.assembler = SimpleAssembler(femModel);
            
            numvarargs = length(varargin);
            optargs = { 0.25 0.5 };
            optargs(1:numvarargs) = varargin;
            [solver.beta, solver.gamma] = optargs{:};
        end
        
        function solve(solver)
            if ~solver.isInitialized
                solver.initialize;
            end
            [freeDofs, fixedDofs] = solver.femModel.getDofConstraints();
            fixedDofIds = fixedDofs.getId();
            
            [forceFull, force] = solver.assembler.applyExternalForces(solver.femModel);
            dispOld = solver.assembler.assemble3dDofVector(solver.femModel, 'DISPLACEMENT');
            dispOld = applyVectorBoundaryConditions(dispOld, fixedDofIds)';
            velOld = solver.assembler.assemble3dValueVector(solver.femModel, 'VELOCITY');
            velOld = applyVectorBoundaryConditions(velOld, fixedDofIds)';
            accOld = solver.assembler.assemble3dValueVector(solver.femModel, 'ACCELERATION');
            accOld = applyVectorBoundaryConditions(accOld, fixedDofIds)';
            
            %set up left hand side of the les
            lhs = solver.massMatrix ...
                + (solver.gamma * solver.dt) .* solver.dampingMatrix ...
                + (solver.beta * power(solver.dt,2)) .* solver.stiffnessMatrix;
            rhs = force' ...
                - solver.dampingMatrix * (velOld + ((1 - solver.gamma) .* solver.dt) * accOld) ...
                - solver.stiffnessMatrix * (dispOld + solver.dt .* velOld + ((0.5 - solver.beta) * power(solver.dt,2)) .* accOld);
            
            accNew = linsolve(lhs, rhs);
            velNew = velOld + ((1 - solver.gamma) * solver.dt) .* accOld + (solver.gamma * solver.dt) .* accNew;
            dispNew = dispOld + solver.dt .* velOld + (power(solver.dt, 2) * (0.5 - solver.beta)) .* accOld + (power(solver.dt, 2) * solver.beta) .* accNew;
            
            %write the values back to the dofs / nodes
            solver.assembler.appendValuesToDofs(solver.femModel, dispNew);
            solver.assembler.appendValuesToNodes(solver.femModel, 'VELOCITY', velNew);
            solver.assembler.appendValuesToNodes(solver.femModel, 'ACCELERATION', accNew);
            
        end
        
        function initialize(solver)
            solver.femModel.initialize;
            solver.femModel.getAllNodes.addNewValue(["VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z", ...
                "ACCELERATION_X", "ACCELERATION_Y", "ACCELERATION_Z"]);
            
            % assemble and reduce matrices
            [freeDofs, fixedDofs] = solver.femModel.getDofConstraints();
            fixedDofIds = fixedDofs.getId();
            
            mass = solver.assembler.assembleGlobalMassMatrix(solver.femModel);
            solver.massMatrix = applyMatrixBoundaryConditions(mass, fixedDofIds);
            damping = solver.assembler.assembleGlobalDampingMatrix(solver.femModel);
            solver.dampingMatrix = applyMatrixBoundaryConditions(damping, fixedDofIds);
            stiffness = solver.assembler.assembleGlobalStiffnessMatrix(solver.femModel);
            solver.stiffnessMatrix = applyMatrixBoundaryConditions(stiffness, fixedDofIds);
            
            % initial acceleration
            [forceFull, force0] = solver.assembler.applyExternalForces(solver.femModel);
            disp = solver.assembler.assemble3dDofVector(solver.femModel, 'DISPLACEMENT');
            disp0 = applyVectorBoundaryConditions(disp, fixedDofIds)';
            vel = solver.assembler.assemble3dValueVector(solver.femModel, 'VELOCITY');
            vel0 = applyVectorBoundaryConditions(vel, fixedDofIds)';
            
            acc0 = (solver.massMatrix) \ (force0' - solver.stiffnessMatrix * disp0 - solver.dampingMatrix * vel0);
            solver.assembler.assignValuesToNodes(solver.femModel, 'ACCELERATION', acc0, 1);
            solver.isInitialized = true;
        end
    end
    
end

