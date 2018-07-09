classdef ModelIO < handle
    %MODELIO Base class for all external input and output functionalities
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

