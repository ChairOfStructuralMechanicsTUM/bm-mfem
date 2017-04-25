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
    function x = solve(femModel)
        Kred = SimpleAssembler(femModel).reducedStiffnessMatrix;
        f = SimpleAssembler(femModel).reducedForceVector;
        
        x = Kred \ f';
        
        SimpleAssembler.assignResultsToDofs(femModel, x);
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

