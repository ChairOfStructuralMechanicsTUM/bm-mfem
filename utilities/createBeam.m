function [ beam ] = createBeam( l, n, varargin )
%CREATEBEAM Creates the mesh of a beam
%   inputs:
%       l: length
%       n: number of elements
%       'elementType': used element type
%       'IY'
%       'IZ'
%       'IT'
%       'YOUNGS_MODULUS'
%       'POISSON_RATIO'
%       'CROSS_SECTION'
%       'DENSITY'
%   outputs:
%       beam: FemModel with all nodes, elements, and properties

p = inputParser();
p.addParameter('elementType','BeamElement3d2n');
p.addParameter('IY',0.00001,@isnumeric);
p.addParameter('IZ',0.00001,@isnumeric);
p.addParameter('IT',0.001,@isnumeric);
p.addParameter('YOUNGS_MODULUS',2.1e11,@isnumeric);
p.addParameter('POISSON_RATIO',0.3,@isnumeric);
p.addParameter('CROSS_SECTION',.001,@isnumeric);
p.addParameter('DENSITY',7850,@isnumeric);
p.parse(varargin{:});
elementType = p.Results.elementType;
iy = p.Results.IY;
iz = p.Results.IZ;
it = p.Results.IT;
young_modulus = p.Results.YOUNGS_MODULUS;
poisson_ratio = p.Results.POISSON_RATIO;
cross_section = p.Results.CROSS_SECTION;
density = p.Results.DENSITY;

beam = FemModel();

%create beam
beam.addNewNode(1,0,0,0);
for ii = 1:n
    beam.addNewNode(ii+1,l*ii/n,0,0);
    beam.addNewElement(elementType,ii,[ii ii+1]);
end

beam.getAllElements.setPropertyValue('IY',iy);
beam.getAllElements.setPropertyValue('IZ',iz);
beam.getAllElements.setPropertyValue('IT',it);
beam.getAllElements.setPropertyValue('YOUNGS_MODULUS',young_modulus);
beam.getAllElements.setPropertyValue('POISSON_RATIO',poisson_ratio);
beam.getAllElements.setPropertyValue('CROSS_SECTION',cross_section);
beam.getAllElements.setPropertyValue('DENSITY',density);

end

