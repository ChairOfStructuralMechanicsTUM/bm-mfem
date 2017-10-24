classdef (Abstract) Hexahedron < Element
    %HEXAHEDRON Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        area
    end
    
    methods
        % constructor
        function hexahedron = Hexahedron(id, material)
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                super_args = {id; material};
            end
            
            hexahedron@Element(super_args{:});
            
            % For hexahedron in comparison to truss no definition of crosssection 
            
        end
        
        % getter functions
        
          % Implement function for calculation of volume here
        
        % setter functions
        
        % member functions
        
          % Implement functions for update, calculation of center and
          % drawing here
        
    end
    
    methods (Access = protected)
       function cp = copyElement(obj)
           cp = copyElement@Element(obj);
           % Extension might be necessary at this point
        end 
    end
    
end

