function [ plate, x0, xl, y0, yl ] = createRectangularPlate( lx, ly, nx, ny, varargin )
%CREATERECTANGULARPLATE Creates the mesh of a rectangular plate
%   inputs:
%       lx: length in x
%       ly: length in y
%       nx: number of elements in x
%       ny: number of elements in y
%       'elementType': used element type
%    outputs:
%        plate: FemModel with all nodes and elements
%        x0: nodes on x=0
%        xl: nodes on x=ly
%        y0: nodes on y=0
%        yl: nodes on y=lx

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

nNodes = (nx+1)*(ny+1);
x0 = plate.getNodes(1:nx+1:nNodes);
y0 = plate.getNodes(1:nx+1);
xl = plate.getNodes(nx+1:nx+1:nNodes);
yl = plate.getNodes(nNodes-nx:nNodes);

end

