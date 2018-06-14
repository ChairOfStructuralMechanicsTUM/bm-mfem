classdef NewtonRaphsoSolvingStrategy < Solver
    %SIMPLESOLVINGSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        femModel
        isInitialized
        
        max_iterations
        tolerance
    end
    
    
    
    methods
        
        %constructor
        function newtonRaphsonSolver = NewtonRaphsonSolvingStrategy(femModel)
            if (nargin > 0)
               newtonRaphsonSolver.femModel = femModel;
               newtonRaphsonSolver.isInitialized = false;
            end
        end
        
        function x = solve(newtonRaphsonSolver)
            if ~ newtonRaphsonSolver.isInitialized
                newtonRaphsonSolver.initialize();
            end
            
            nIter = 0;
            while (nIter < newtonRaphsonSolver.max_iterations && residual_norm > newtonRaphsonSolver.tolerance)
                
               
                nIter = nIter + 1; 
            end
            

                
                
                
            
                
            end

            
             
%             [~, Kred] = SimpleAssembler.assembleGlobalStiffnessMatrix(simpleSolver.femModel);
%             [~, fred] = SimpleAssembler.applyExternalForces(simpleSolver.femModel);
%             
%             x = linsolve(Kred, fred.');
%             
%             SimpleAssembler.assignResultsToDofs(simpleSolver.femModel, x);
        end
        
        function initialize(newtonRaphsonSolver)
            newtonRaphsonSolver.femModel.initialize;
            newtonRaphsonSolver.isInitialized = true;
        end
        

end

