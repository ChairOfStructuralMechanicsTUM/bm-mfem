clc;
clear all;

model = FemModel();

% Dimension of the structure
Lx=10;
Ly=2;

% Number of elements in specific directions
nx=20;
ny=4;

% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;

numnodes=(nx + 1)*(ny + 1);
numele=nx*ny;

% Generation of nodes
id=0;

for j=1:(ny + 1)
    for i=1:(nx + 1)
        id=id+1;
        model.addNewNode(id,(i-1)*dx,(j-1)*dy);
    end
end

% Assignment DOFs
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

% Generation of elements
id = 0;

for j=1:ny
    for i=1:nx
        id=id+1;
        a = i + (j-1)*(nx+1);
        model.addNewElement('QuadrilateralElement2d4n',id,[a, a+1, a+1+(nx+1), a+(nx+1)]);
    end
end

% assignment of material properties
model.getAllElements.setPropertyValue('YOUNGS_MODULUS',30e6);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.0);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',3);
model.getAllElements.setPropertyValue('DENSITY',7860);

% Definition of BCs
for i=1:(ny+1)
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_X');
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_Y');
end

% Definition of loading
addPointLoad(model.getNode(numnodes),300,[0 -1]);

% Determination of global matrices
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);           
massMatrix = assembling.assembleGlobalMassMatrix(model);

% Solving
solver = SimpleSolvingStrategy(model);
x = solver.solve();

step = 1;
VerschiebungDofs = model.getDofArray.getValue(step);
nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(1000);
v.plotUndeformed
v.plotDeformed
    
