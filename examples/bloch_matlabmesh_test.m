
clear all;
model = FemModel();

% Dimension of the structure
Lx=0.1;
Ly=0.02;

% Number of elements in specific directions
nx=3;
ny=2;

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
model.getAllElements.setPropertyValue('YOUNGS_MODULUS',78*10^6);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',70000);

solver = BlochInverse1D(model);
assembling = SimpleAssembler(model);

stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);            
massMatrix = assembling.assembleGlobalMassMatrix(model);

initialize(solver)
[Ksorted,Msorted] = sortKandM(solver,stiffnessMatrix,massMatrix);
numberOfPhases = 10;
[Kred,Mred] = reducedStiffnesAndMass (solver,Ksorted,Msorted,numberOfPhases);  
omega = cell(numberOfPhases,1);
nob = 3;
[kx,miu] = propConst(solver,numberOfPhases); 

    
for j = 1:nob
    for i = 1:numberOfPhases
        omega{i,1} = solver.calcOmega(Kred{i,1},Mred{i,1},nob);
        f(j,i) = omega{i,1}(j,1)/(2*pi);
    end
    figure(1)
    plot(kx,f(j,:),'r')
    title('Dispersion curves')
    xlabel('Phase k')
    ylabel('frequenzy f')
    xlim([0 pi])
    legend(['bandnumber: ' num2str(j)],'Location','EastOutside')
    hold on
end