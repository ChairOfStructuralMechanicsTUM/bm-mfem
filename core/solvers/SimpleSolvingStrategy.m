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
    function [u01, u02] = solve(femModel01, femModel02)
        
        %Feti solution
        if nargin == 2
            disp('FETI');
            K_01 = SimpleAssembler(femModel01).reducedStiffnessMatrix;
            f_01 = SimpleAssembler(femModel01).reducedForceVector;
            
            K_02 = SimpleAssembler(femModel02).reducedStiffnessMatrix;
            f_02 = SimpleAssembler(femModel02).reducedForceVector;
            
            [u01, u02] = FetiSolver.solveFeti(K_01, K_02, f_01, f_02, femModel01, femModel02);
            
            SimpleAssembler.assignResultsToDofs(femModel01, u01);
            SimpleAssembler.assignResultsToDofs(femModel02, u02);
            
        %FEM solution    
        else
            disp('FEM');
            Kred = SimpleAssembler(femModel01).reducedStiffnessMatrix;
            f = SimpleAssembler(femModel01).reducedForceVector;
            u01 = Kred \ f';
            SimpleAssembler.assignResultsToDofs(femModel01, u01);
        end
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

