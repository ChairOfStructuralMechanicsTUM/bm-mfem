classdef HexahedronElement < Element %Class HexahedronElement to be implemented if necessary
    %HexahedronElement Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
    end
    
    methods
        % constructor
        function hexahedron = HexahedronElement(id, nodeArray, requiredProperties)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id, nodeArray, requiredProperties};
            end
            
            hexahedron@Element(super_args{:});
        end
    end
    
end