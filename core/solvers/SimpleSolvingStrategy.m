classdef SimpleSolvingStrategy < Solver
    %SIMPLESOLVINGSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
%         nodalDisplacments
    end
    


    %     methods
        
        %     constructor
%         function nodalDisplacements = SimpleSolvingStrategy(femModel)
%             if (nargin > 0)
%                 nodalDisplacements.res = SimpleSolvingStrategy.solve(femModel);
%             end
%         end
%     end
    
    
    
methods (Static)
    function x = solve(femModel_01, femModel_02)
        
        %Feti solution
        if nargin == 2
            disp('FETI');
            K_01 = SimpleAssembler(femModel_01).reducedStiffnessMatrix;
            f_01 = SimpleAssembler(femModel_01).reducedForceVector;
            
            K_02 = SimpleAssembler(femModel_02).reducedStiffnessMatrix;
            f_02 = SimpleAssembler(femModel_02).reducedForceVector;
            
            x = FetiSolver.solveFeti(K_01, K_02, f_01, f_02);
            
        %FEM solution    
        else
            disp('FEM');
            Kred = SimpleAssembler(femModel_01).reducedStiffnessMatrix;
            f = SimpleAssembler(femModel_01).reducedForceVector;
            x = Kred \ f';
        end
        
        SimpleAssembler.assignResultsToDofs(femModel_01, x);
    end
 
    function nodalForces = getNodalForces(femModel)
        
        dofs = femModel.getDofArray;
        nDofs = length(dofs);
        nodalForces = SimpleAssembler(femModel).forceVector;
        fixedDofs = [];
        
        for itDof = 1:nDofs
            if (dofs(itDof).isFixed)
                fixedDofs = [fixedDofs itDof];                          % array of fixed dofs and their location
            else
                nodalForces(itDof) = dofs(itDof).getValue;
            end
        end
        
        nodalForces(fixedDofs) = SimpleAssembler(femModel).stiffnessMatrix(fixedDofs, :) * getValue(femModel.dofArray)';
        
    end
    
end
end

