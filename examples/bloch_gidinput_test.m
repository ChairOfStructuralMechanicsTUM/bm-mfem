clear
close all

io=MdpaInput('cell_circle.mdpa'); %specify input file
model = io.readModel(); %read the model

model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);

innerNodes = model.getModelPart('GENERIC_Inner_circle').getNodes();
id = getId(innerNodes); 
coords = getCoords(innerNodes);
% leftNodes2 = model.getModelPart('GENERIC_leftNodes2').getNodes()
% rightNodes = model.getModelPart('GENERIC_rightNodes').getNodes()
% rightNodes2 = model.getModelPart('GENERIC_rightNodes2').getNodes()
% innerNodes = model.getModelPart('GENERIC_innerNodes').getNodes()
% innerNodes2 = model.getModelPart('GENERIC_innerNodes2').getNodes()

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);

% model.addNewNode(434,0.015,0.01,0);
% model.addNewElement('SpringDamperElement3d2n',367,[254 434]);
% model.addNewElement('SpringDamperElement3d2n',368,[434 115]);
% model.addNewElement('ConcentratedMassElement3d1n',369, 434);
% 
% model.getNode(254).addDof("DISPLACEMENT_Z");
% model.getNode(115).addDof("DISPLACEMENT_Z");
% model.getNode(434).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"]);
% 
% model.getNode(254).fixDof('DISPLACEMENT_Z');
% model.getNode(115).fixDof('DISPLACEMENT_Z');
% model.getNode(434).fixDof('DISPLACEMENT_Y');
% model.getNode(434).fixDof('DISPLACEMENT_Z');
% 
% model.getElement(367).setPropertyValue('ELEMENTAL_STIFFNESS',7e5);
% model.getElement(367).setPropertyValue('ELEMENTAL_DAMPING',0);
% model.getElement(368).setPropertyValue('ELEMENTAL_STIFFNESS',7e5);
% model.getElement(368).setPropertyValue('ELEMENTAL_DAMPING',0);
% 
% model.getElement(369).setPropertyValue('ELEMENTAL_MASS',7e2);
% model.getElement(369).setPropertyValue('VOLUME_ACCELERATION',10);
% 
solver = BlochInverse1D(model);
assembling = SimpleAssembler(model);

[stiffnessMatrix,Kred1] = assembling.assembleGlobalStiffnessMatrix(model);            
[massMatrix,Mred1] = assembling.assembleGlobalMassMatrix(model);
% stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);            
% massMatrix = assembling.assembleGlobalMassMatrix(model);

initialize(solver)
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

% fig=figure
% v=Visualization(model); 
% v.plotUndeformed()  