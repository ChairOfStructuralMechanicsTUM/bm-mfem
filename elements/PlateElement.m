classdef (Abstract) PlateElement < Element
    
    %PLATEELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        lengthX
        lengthY
    end
    
    methods
        % Constructor
        function plateElement = PlateElement(id)
            if nargin == 0
                super_args = {};
            elseif nargin == 1
                super_args = {id};
            end
            
            plateElement@Element(super_args{:});
        end
        
        %getter functions
        function lenX = getLengthX(plateElement)
        if plateElement.lengthX == 0 
            plateElement.lengthX = computeLength(plateElement.nodeArray(1).getCoords,...
                plateElement.nodeArray(2).getCoords);
        end
        lenX = plateElement.lengthX;
        end
        
        function lenY = getLengthY(plateElement)
        if plateElement.lengthY == 0 
            plateElement.lengthY = computeLength(plateElement.nodeArray(1).getCoords,...
                plateElement.nodeArray(4).getCoords);
        end
        lenY = plateElement.lengthY;
        end
        
        % Check wether Geometry is a convex Quadrilateral
        function checkGeometry(plateElement)
            diag1X = [plateElement.nodeArray(1).getX() plateElement.nodeArray(3).getX()];
            diag1Y = [plateElement.nodeArray(1).getY() plateElement.nodeArray(3).getY()];
            
            diag2X = [plateElement.nodeArray(2).getX() plateElement.nodeArray(4).getX()];
            diag2Y = [plateElement.nodeArray(2).getY() plateElement.nodeArray(4).getY()];
            
            intersection = polyxpoly(diag1X, diag1Y, diag2X, diag2Y);
            if  ~isempty(intersection)
                disp("Geometry is convex!");
            else
                disp("Geometry is not convex!");
            end
        end

        function update(plateElement)
        end
        
        function computeLocalForceVector
        end

        function barycenter
        end

    end
    
        
        
        
    
    
end
