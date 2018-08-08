clear;

node01 = Node(1,0,0);
node02 = Node(2,2,0);
node03 = Node(3,4,0);
node04 = Node(4,0,1);
node05 = Node(5,2,1);
node06 = Node(6,4,1);
node07 = Node(7,0,2);
node08 = Node(8,2,2);
node09 = Node(9,4,2);

nodeArray = [node01 node02 node03 node04 node05 ...
    node06 node07 node08 node09];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

ele01 = QuadrilateralElement2d4n(1,[node01 node02 node05 node04]);
ele02 = QuadrilateralElement2d4n(2,[node02 node03 node06 node05]);
ele03 = QuadrilateralElement2d4n(3,[node04 node05 node08 node07]);
ele04 = QuadrilateralElement2d4n(4,[node05 node06 node09 node08]);

elementArray = [ele01 ele02 ele03 ele04];

elementArray.setPropertyValue('YOUNGS_MODULUS',96);
elementArray.setPropertyValue('POISSON_RATIO',1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);
elementArray.setPropertyValue('DENSITY',7860);

model = FemModel(nodeArray, elementArray);

model.getNode(5).fixDof('DISPLACEMENT_Y');
model.getNode(5).fixDof('DISPLACEMENT_X');

solver = BlochInverse1D(model);
assembling = SimpleAssembler(model);

[stiffnessMatrix,Kred1] = assembling.assembleGlobalStiffnessMatrix(model);
            
[massMatrix,Mred1] = assembling.assembleGlobalMassMatrix(model);

% v = Visualization(model);
% v.setScaling(1);
% v.plotUndeformed

initialize(solver)

%  [Ksorted,Msorted] = sortKandM(solver,stiffnessMatrix,massMatrix);
 [Ksorted,Msorted] = sortKandM(solver,Kred1,Mred1);

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

    figure(2)
    plot(kx,f(j,:),'r')
    title('Dispersion curves')
    xlabel('Phase k')
    ylabel('frequenzy f')
    xlim([0 pi])
    hold on
    

end