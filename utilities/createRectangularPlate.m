function [ plate ] = createRectangularPlate( lx, ly, nx, ny, varargin )
%CREATERECTANGULARPLATE Creates the mesh of a rectangular plate
%   lx: length in x
%   ly: length in y
%   nx: number of elements in x
%   ny: number of elements in y
%   'elementType': used element type

p = inputParser();
p.addParameter('elementType','ReissnerMindlinElement3d4n');
p.parse(varargin{:});
elementType = p.Results.elementType;

plate = FemModel();

%create nodes
idn = 1;
for idy = 1:(ny+1)
    y = (idy-1) * ly / ny;
    for idx = 1:(nx+1)
        x = (idx-1) * lx / nx;
        plate.addNewNode(idn,x,y,0);
        idn = idn + 1;
    end
end

%create elements
ide = 1;
for idy = 1:ny
    for idx = 1:nx
        plate.addNewElement(elementType,ide, ...
            [idx+(idy-1)*(nx+1) idx+(idy-1)*(nx+1)+1 idx+idy*(nx+1)+1 idx+idy*(nx+1)]);
        ide = ide + 1;
    end
end

end

