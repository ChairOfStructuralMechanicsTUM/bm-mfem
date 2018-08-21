%%%% A = 20 cm^2 b*h = 2cm*10cm 
    
clear

io=MdpaInput('SpringHorizontal_rechteckig_hoch.mdpa'); %specify input file
model = io.readModel(); %read the model

model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
a=model.getAllModelParts
leftNodes = model.getModelPart('GENERIC_leftNodes').getNodes();
rightNodes = model.getModelPart('GENERIC_rightNodes').getNodes();
springNodes = model.getModelPart('GENERIC_springNodes').getNodes()
%     innerMaterial = model.getModelPart('GENERIC_innerMaterial').getNodes()

% coords = getCoords(node)
% springNodeIds = arrayfun(@(node) node.getId, springNodes);
% springNodeCoords = arrayfun(@(node) node.getCoords, springNodes,'UniformOutput',false);
% 
% leftSpringNode = springNodes(1,7);
% leftSNCoords = getCoords(leftSpringNode);
% leftSNId = getId(leftSpringNode);
% rightSpringNode = springNodes(1,2);
% rightSNCoords = getCoords(rightSpringNode);
% rightSNId = getId(rightSpringNode);
% horizontalSpringNodes = [leftSpringNode rightSpringNode];

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);

%%% 2 springs, 1 mass
% % model.addNewNode(366,0.1,0.075,0);
% % model.addNewElement('SpringDamperElement3d2n',308,[112 366]);
% % model.addNewElement('SpringDamperElement3d2n',309,[366 198]);
% % model.addNewElement('ConcentratedMassElement3d1n',310, 366);
% % 
% % model.getNode(112).addDof("DISPLACEMENT_Z");
% % model.getNode(198).addDof("DISPLACEMENT_Z");
% % model.getNode(366).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"]);
% % 
% % 
% % model.getNode(112).fixDof('DISPLACEMENT_Z');
% % model.getNode(198).fixDof('DISPLACEMENT_Z');
% % model.getNode(366).fixDof('DISPLACEMENT_Y');
% % model.getNode(366).fixDof('DISPLACEMENT_Z');
% % 
% % 
% % model.getElement(308).setPropertyValue('ELEMENTAL_STIFFNESS',7e5);
% % model.getElement(308).setPropertyValue('ELEMENTAL_DAMPING',0);
% % model.getElement(309).setPropertyValue('ELEMENTAL_STIFFNESS',7e5);
% % model.getElement(309).setPropertyValue('ELEMENTAL_DAMPING',0);
% % 
% % model.getElement(310).setPropertyValue('ELEMENTAL_MASS',7e2);
% % model.getElement(310).setPropertyValue('VOLUME_ACCELERATION',10);


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

for j = 1:nob

    for i = 1:numberOfPhases
        omega{i,1} = solver.calcOmega(Kred{i,1},Mred{i,1},nob);
        f(j,i) = omega{i,1}(j,1)/(2*pi);

    end

    figure(2)
    plot(kx,f(j,:),'r')
    title('Dispersion curves - StandardRechteck,hoch')
    xlabel('Phase k')
    ylabel('frequenzy f')
    xlim([0 pi])
    legend(['bandnumbers: ' num2str(j)],'Location','EastOutside')
    hold on
end
