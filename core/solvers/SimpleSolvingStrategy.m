classdef SimpleSolvingStrategy < Solver
    %SIMPLESOLVINGSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function solve(femModel)
            K = SimpleAssembler.assembleGlobalStiffnessMatrix(femModel);
            f = SimpleAssembler.applyExternalForces(femModel);
            
            x = K \ f';
            
            SimpleAssembler.assignResultsToDofs(femModel, x);
        end
    end
    
end

