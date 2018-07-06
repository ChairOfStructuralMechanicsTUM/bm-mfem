classdef ModelIO < handle
    %MODELIO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        file
    end
    
    methods
        function modelIO = ModelIO(file)
            if nargin > 0
                modelIO.file = file;
            end
        end
    end
    
end

