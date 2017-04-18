classdef Node < handle %handle class
    %NODE The node class
    %   Parameters:
    %       id: unique identifier
    %       x, y(, z): coordinates
    
    properties (Access = private)
        id
        x
        y
        z
    end
    
    methods
        % constructor
        function node = Node(id, x, y, z)
            switch nargin
                case 3
                    node.id = id;
                    node.x = x;
                    node.y = y;
                case 4
                    node.id = id;
                    node.x = x;
                    node.y = y;
                    node.z = z;
                otherwise
                    error('Wrong number of arguments')
            end
        end
        
        % getter functions
        function id = getId(node)
            id = node.id;
        end
        
        function coords = getCoords(node)
           coords = [node.x node.y node.z]; 
        end
        
        
        
    end
    
    methods (Access = private)
        
    end
    
end

