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
    function [u] = solve(substructures)
        
        %distinguish between FETI and FEM by number of subdomains
        if length(substructures) == 1
            %FEM
            disp('FEM');
            Kred = SimpleAssembler(substructures).reducedStiffnessMatrix;
            f = SimpleAssembler(substructures).reducedForceVector;
          
            u = Kred \ f';
            SimpleAssembler.assignResultsToDofs(substructures, u);
        else
            %FETI
            disp('FETI');
            K = cell(1,length(substructures));
            f = cell(1,length(substructures));
            %calculate stiffness Matrix and Force Vector for all
            %substructues and save them in a cell array
            for ii = 1:length(substructures)
                K{1,ii} = [SimpleAssembler(substructures(ii)).reducedStiffnessMatrix];
                f{1,ii} = [SimpleAssembler(substructures(ii)).reducedForceVector];
            end

            u = FetiPreparer.solveFeti(K, f, substructures);
            
            %SimpleAssembler.assignResultsToDofs(substructures, u);
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
                %nodalForces(itDof) = dofs(itDof).getValue;
                nodalForces(itDof) = dofs(itDof).getDofLoad;
            end
        end
        
        nodalForces(fixedDofs) = SimpleAssembler(femModel).stiffnessMatrix(fixedDofs, :) * getValue(femModel.dofArray)';
        
    end
    
end
end

