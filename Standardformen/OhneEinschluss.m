%%% Ages=300cm^2 b=20cm h=15cm

io=MdpaInput('OhneEinschluss.mdpa'); %specify input file
model = io.readModel(); %read the model

model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
a=model.getAllModelParts;


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

Dofs = model.getDofArray();

initialize(solver)
[Ksorted,Msorted] = sortKandM(solver,Kred1,Mred1);
numberOfPhases = 50;
[Kred,Mred] = reducedStiffnesAndMass (solver,Ksorted,Msorted,numberOfPhases);  
omega = cell(numberOfPhases,1);
nob = 10;
[kx,miu] = propConst(solver,numberOfPhases); %already used in reducedStiffnesAndMass(..)
    
for j = 1:nob
    for i = 1:numberOfPhases
        omega{i,1} = solver.calcOmega(Kred{i,1},Mred{i,1},nob);
        f(j,i) = omega{i,1}(j,1)/(2*pi);
    end
    figure(2)
    plot(kx,f(j,:),'r')
    title('Dispersion curves - OhneEinschluss')
    xlabel('Phase k')
    ylabel('frequenzy f')
    xlim([0 pi])
    legend(['bandnumber: ' num2str(j)],'Location','EastOutside')
    hold on
end


