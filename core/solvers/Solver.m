classdef Solver < handle
    %SOLVER Base class for solver implementations
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Abstract)
        initialize(solver)
        solve(solver)
    end
    
end

