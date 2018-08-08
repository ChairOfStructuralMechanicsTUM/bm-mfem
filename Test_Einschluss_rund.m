
clear

io=MdpaInput('Test_Einschluss_rund.mdpa'); %specify input file
model = io.readModel(); %read the model

model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
a=model.getAllModelParts;
innerCircle = model.getModelPart('GENERIC_innerCircle').getNodes();
id = getId(innerCircle) ;
coords = getCoords(innerCircle);

x = getX(innerCircle);
leftNode = model.getNode(66);
rightNode = model.getNode(130);
leftcoords = getCoords(leftNode);
rightcoords = getCoords(rightNode);

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);


model.addNewNode(239,0.1,0.075,0);
model.addNewElement('SpringDamperElement3d2n',199,[66 239]);
model.addNewElement('SpringDamperElement3d2n',200,[239 130]);
model.addNewElement('ConcentratedMassElement3d1n',201, 239);

model.getNode(66).addDof("DISPLACEMENT_Z");
model.getNode(130).addDof("DISPLACEMENT_Z");
model.getNode(239).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"]);

model.getNode(66).fixDof('DISPLACEMENT_Z');
model.getNode(130).fixDof('DISPLACEMENT_Z');
model.getNode(239).fixDof('DISPLACEMENT_Y');
model.getNode(239).fixDof('DISPLACEMENT_Z');

model.getElement(199).setPropertyValue('ELEMENTAL_STIFFNESS',7e5);
model.getElement(199).setPropertyValue('ELEMENTAL_DAMPING',0);
model.getElement(200).setPropertyValue('ELEMENTAL_STIFFNESS',7e5);
model.getElement(200).setPropertyValue('ELEMENTAL_DAMPING',0);

model.getElement(201).setPropertyValue('ELEMENTAL_MASS',7e2);
model.getElement(201).setPropertyValue('VOLUME_ACCELERATION',10);

v=Visualization(model); %set up visualization
v.plotUndeformed()  %visualize

solver = BlochInverse1D(model);
assembling = SimpleAssembler(model);

[stiffnessMatrix,Kred1] = assembling.assembleGlobalStiffnessMatrix(model);            
[massMatrix,Mred1] = assembling.assembleGlobalMassMatrix(model);

Dofs = model.getDofArray();

initialize(solver)
% [Ksorted,Msorted] = sortKandM(solver,stiffnessMatrix,massMatrix);
[Ksorted,Msorted] = sortKandM(solver,Kred1,Mred1);
numberOfPhases = 20;
[Kred,Mred] = reducedStiffnesAndMass (solver,Ksorted,Msorted,numberOfPhases);  
omega = cell(numberOfPhases,1);
nob = 10;
[kx,miu] = propConst(solver,numberOfPhases); %already used in reducedStiffnesAndMass(..)
    
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


