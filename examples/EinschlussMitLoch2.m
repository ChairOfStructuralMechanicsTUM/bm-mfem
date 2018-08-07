clear all;
model = FemModel();

% % ELEMENT1

Lx=0.2;
Ly=0.055;
% Number of elements in specific directions
nx=20;
ny=11;
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
        model.addNewElement('QuadrilateralElement2d4n',id,[a, a+1, a+1+(nx+1), a+(nx+1)]); %letzter Teil: KnotenIDS
    end
end





% % % ELEMENT2

Lx=0.07;
Ly=0.03;
% Number of elements in specific directions
nx=7;
ny=6;
% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;
numnodes=(nx + 1)*(ny + 1);
numele=nx*ny;
% Generation of nodes
id=252;
for j=1:(ny + 1)
    for i=1:(nx + 1)
        id=id+1;
        model.addNewNode(id,(i-1)*dx,(j-1)*dy+0.06);
    end
end

% Assignment DOFs
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
% Generation of elements



id = 220;
for j=1:ny
    for i=1:nx
        id=id+1;
        a = i + (j-1)*(nx+1)+252;
        model.addNewElement('QuadrilateralElement2d4n',id,[a, a+1, a+1+(nx+1), a+(nx+1)]);
    end
end


% % % % %ELEMENT3

Lx=0.07;
Ly=0.03;
% Number of elements in specific directions
nx=7;
ny=6;
% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;
numnodes=(nx + 1)*(ny + 1);
numele=nx*ny;
% Generation of nodes
id=308;
for j=1:(ny + 1)
    for i=1:(nx + 1)
        id=id+1;
        model.addNewNode(id,(i-1)*dx+0.13,(j-1)*dy+0.06);
    end
end

% Assignment DOFs
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
% Generation of elements

id = 262;

for j=1:ny
    for i=1:nx
        id=id+1;
        a = i + (j-1)*(nx+1)+308;
        model.addNewElement('QuadrilateralElement2d4n',id,[a, a+1, a+1+(nx+1), a+(nx+1)]);
    end
end


%ELEMENT4

Lx=0.2;
Ly=0.055;
% Number of elements in specific directions
nx=20;
ny=11;
% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;
numnodes=(nx + 1)*(ny + 1);
numele=nx*ny;
% Generation of nodes
id=364;
for j=1:(ny + 1)
    for i=1:(nx + 1)
        id=id+1;
        model.addNewNode(id,(i-1)*dx,(j-1)*dy+0.095);
    end
end

% Assignment DOFs
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
% Generation of elements

id = 304;
for j=1:ny
    for i=1:nx
        id=id+1;
        a = i + (j-1)*(nx+1)+364;
        model.addNewElement('QuadrilateralElement2d4n',id,[a, a+1, a+1+(nx+1), a+(nx+1)]);
    end
end



%ELEMENT 1-2

Lx=0.07;
Ly=0.005;
% Number of elements in specific directions
nx=7;
ny=1;
% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;
numnodes=(nx + 1)*(ny + 1);
numele=nx*ny;
a=0;
id = 524;
for i=1:nx
    id=id+1;
    a = i + 231;
    model.addNewElement('QuadrilateralElement2d4n',id,[a, a+1, a+22, a+21]);
end



%ELEMENT 1-3

Lx=0.07;
Ly=0.005;
% Number of elements in specific directions
nx=7;
ny=1;
% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;
numnodes=(nx + 1)*(ny + 1);
numele=nx*ny;
a=0;
id = 531;
for i=1:nx
    id=id+1;
    a = i + 244;
    model.addNewElement('QuadrilateralElement2d4n',id,[a, a+1, a+65, a+64]);
end


% % %ELEMENT 2-4

Lx=0.07;
Ly=0.005;
% Number of elements in specific directions
nx=7;
ny=1;
% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;
numnodes=(nx + 1)*(ny + 1);
numele=nx*ny;
a=0;
id = 538;
for i=1:nx
    id=id+1;
    a = i + 300;
    model.addNewElement('QuadrilateralElement2d4n',id,[a, a+1, a+65, a+64]);
end

% % % Element 3-4

Lx=0.07;
Ly=0.005;
% Number of elements in specific directions
nx=7;
ny=1;
% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;
numnodes=(nx + 1)*(ny + 1);
numele=nx*ny;
a=0;
id = 545;
for i=1:nx
    id=id+1;
    a = i + 356;
    model.addNewElement('QuadrilateralElement2d4n',id,[a, a+1, a+22, a+21]);
end



model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7*10^10);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);



v = Visualization(model);
v.setScaling(1);
v.plotUndeformed

solver = BlochInverse1D(model);
assembling = SimpleAssembler(model);

stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);


initialize(solver)


[Ksorted,Msorted] = sortKandM(solver,stiffnessMatrix,massMatrix);

numberOfPhases = 20;

[Kred,Mred] = reducedStiffnesAndMass (solver,Ksorted,Msorted,numberOfPhases);  

omega = cell(numberOfPhases,1);

nob = 10;
[kx,miu] = propConst(solver,numberOfPhases); %already used in reducedStiffnesAndMass(..)

    
for j = 1:nob
%     f(j,1) = zeros(1,numberOfPhases); 
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
    hold on
    

end




