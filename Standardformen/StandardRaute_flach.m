%%%Raute mit A=20cm (Diagonalenlängen d1=15 und d2=8/3cm)
    
clear

io=MdpaInput('SpringHorizontal_raute_flach.mdpa'); %specify input file
model = io.readModel(); %read the model

model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
a=model.getAllModelParts
leftNodes = model.getModelPart('GENERIC_leftNodes').getNodes();
rightNodes = model.getModelPart('GENERIC_rightNodes').getNodes();
springNodes = model.getModelPart('GENERIC_springNodes').getNodes()
%     innerMaterial = model.getModelPart('GENERIC_innerMaterial').getNodes()

% coords = getCoords(node)


model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);


v=Visualization(model); %set up visualization
v.plotUndeformed()  %visualize


solver = BlochInverse1D(model);
assembling = SimpleAssembler(model);

[stiffnessMatrix,Kred1] = assembling.assembleGlobalStiffnessMatrix(model);
[massMatrix,Mred1] = assembling.assembleGlobalMassMatrix(model);

initialize(solver)
[Ksorted,Msorted] = sortKandM(solver,Kred1,Mred1);

numberOfPhases = 20;

[Kred,Mred] = reducedStiffnesAndMass (solver,Ksorted,Msorted,numberOfPhases);  

omega = cell(numberOfPhases,1);

nob = 10;
[kx,miu] = propConst(solver,numberOfPhases);
figure(2)

for j = 1:nob

    for i = 1:numberOfPhases
        omega{i,1} = solver.calcOmega(Kred{i,1},Mred{i,1},nob);
        f(j,i) = omega{i,1}(j,1)/(2*pi);

    end
    
    plot(kx,f(j,:),'r')
    title('Dispersion curves - S-Raute, flach')
    xlabel('Phase k')
    ylabel('frequenzy f')
    xlim([0 pi])
    legend(['bandnumbers: ' num2str(j)],'Location','EastOutside')
    hold on
end
