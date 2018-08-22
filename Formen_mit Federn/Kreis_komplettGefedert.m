%%%% Kreis, Masse_Feder_system 90Grad
    
clear

io=MdpaInput('Federtests_Kreis.mdpa'); %specify input file
model = io.readModel(); %read the model
model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);

% % 
a=model.getAllModelParts;

allNodes = model.getAllNodes();
massNodeID = length(allNodes)+1;

allElements = model.getAllElements();
spring1ID = length(allElements)+1;
spring2ID = length(allElements)+2;
massID = length(allElements)+3;

spring3ID = length(allElements)+3;
spring4ID = length(allElements)+4;
spring5ID = length(allElements)+5;
spring6ID = length(allElements)+6;
spring7ID = length(allElements)+7;
spring8ID = length(allElements)+8;
spring9ID = length(allElements)+9;
spring10ID = length(allElements)+10;
spring11ID = length(allElements)+11;
spring12ID = length(allElements)+12;
spring13ID = length(allElements)+13;
spring14ID = length(allElements)+14;
spring15ID = length(allElements)+15;
spring16ID = length(allElements)+16;



springNodes = model.getModelPart('GENERIC_0Grad').getNodes();

leftSpringNode = springNodes(1,1);
% % leftSNCoords = getCoords(leftSpringNode);
leftSNId = getId(leftSpringNode);
rightSpringNode = springNodes(1,2);
% % rightSNCoords = getCoords(rightSpringNode);
rightSNId = getId(rightSpringNode);


% %% 2 springs, 1 mass
model.addNewNode(massNodeID,0.1,0.075,0);
model.addNewElement('SpringDamperElement3d2n',spring1ID,[leftSNId massNodeID]);
model.addNewElement('SpringDamperElement3d2n',spring2ID,[massNodeID rightSNId]);
model.addNewElement('ConcentratedMassElement3d1n',massID, massNodeID);

model.getNode(leftSNId).addDof("DISPLACEMENT_Z");
model.getNode(rightSNId).addDof("DISPLACEMENT_Z");
model.getNode(massNodeID).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"]);


model.getNode(leftSNId).fixDof('DISPLACEMENT_Z');
model.getNode(rightSNId).fixDof('DISPLACEMENT_Z');
model.getNode(massNodeID).fixDof('DISPLACEMENT_Y');
model.getNode(massNodeID).fixDof('DISPLACEMENT_Z');


model.getElement(spring1ID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
model.getElement(spring1ID).setPropertyValue('ELEMENTAL_DAMPING',0);
model.getElement(spring2ID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
model.getElement(spring2ID).setPropertyValue('ELEMENTAL_DAMPING',0);

model.getElement(massID).setPropertyValue('ELEMENTAL_MASS',7e0);
model.getElement(massID).setPropertyValue('VOLUME_ACCELERATION',10);



%%%%%%%%Spring 3 und 4 ; 22.5Grad

springNodes = model.getModelPart('GENERIC_22.5Grad').getNodes();

leftSpringNode = springNodes(1,1);
% % leftSNCoords = getCoords(leftSpringNode);
leftSNId = getId(leftSpringNode);
rightSpringNode = springNodes(1,2);
% % rightSNCoords = getCoords(rightSpringNode);
rightSNId = getId(rightSpringNode);

model.addNewElement('SpringDamperElement3d2n',spring3ID,[leftSNId massNodeID]);
model.addNewElement('SpringDamperElement3d2n',spring4ID,[massNodeID rightSNId]);

model.getNode(leftSNId).addDof("DISPLACEMENT_Z");
model.getNode(rightSNId).addDof("DISPLACEMENT_Z");

model.getNode(leftSNId).fixDof('DISPLACEMENT_Z');
model.getNode(rightSNId).fixDof('DISPLACEMENT_Z');

model.getElement(spring3ID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
model.getElement(spring3ID).setPropertyValue('ELEMENTAL_DAMPING',0);
model.getElement(spring4ID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
model.getElement(spring4ID).setPropertyValue('ELEMENTAL_DAMPING',0);



%%%%%%%%Spring 5 und 6 ; 45Grad

springNodes = model.getModelPart('GENERIC_45Grad').getNodes();
leftSpringNode = springNodes(1,1);
leftSNId = getId(leftSpringNode);
rightSpringNode = springNodes(1,2);
rightSNId = getId(rightSpringNode);

model.addNewElement('SpringDamperElement3d2n',spring5ID,[leftSNId massNodeID]);
model.addNewElement('SpringDamperElement3d2n',spring6ID,[massNodeID rightSNId]);

model.getNode(leftSNId).addDof("DISPLACEMENT_Z");
model.getNode(rightSNId).addDof("DISPLACEMENT_Z");

model.getNode(leftSNId).fixDof('DISPLACEMENT_Z');
model.getNode(rightSNId).fixDof('DISPLACEMENT_Z');

model.getElement(spring5ID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
model.getElement(spring5ID).setPropertyValue('ELEMENTAL_DAMPING',0);
model.getElement(spring6ID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
model.getElement(spring6ID).setPropertyValue('ELEMENTAL_DAMPING',0);















v=Visualization(model); %set up visualization
v.plotUndeformed()  %visualize


% solver = BlochInverse1D(model);
% assembling = SimpleAssembler(model);
% 
% [stiffnessMatrix,Kred1] = assembling.assembleGlobalStiffnessMatrix(model);
% [massMatrix,Mred1] = assembling.assembleGlobalMassMatrix(model);
% 
% initialize(solver)
% [Ksorted,Msorted] = sortKandM(solver,Kred1,Mred1);
% 
% numberOfPhases = 20;
% 
% [Kred,Mred] = reducedStiffnesAndMass (solver,Ksorted,Msorted,numberOfPhases);  
% 
% omega = cell(numberOfPhases,1);
% 
% nob = 10;
% [kx,miu] = propConst(solver,numberOfPhases);
% 
% 
% figure(2)
% title('Dispersion curves, Kreis 90Grad')
% xlabel('Phase k')
% ylabel('frequenzy f')
% xlim([0 pi])
% 
% ylim([0 2e4])
% 
% hold on
% for j = 1:nob
%     
%     for i = 1:numberOfPhases
%         omega{i,1} = solver.calcOmega(Kred{i,1},Mred{i,1},nob);
%         f(j,i) = omega{i,1}(j,1)/(2*pi);
% 
%     end
%     plot(kx,f(j,:),'r')       
%     legend(['bandnumbers: ' num2str(j)],'Location','EastOutside')
% end
% 
% 
% m = model.getElement(massID).getPropertyValue('ELEMENTAL_MASS');
% k1 = model.getElement(spring1ID).getPropertyValue('ELEMENTAL_STIFFNESS');
% k2 = model.getElement(spring2ID).getPropertyValue('ELEMENTAL_STIFFNESS');
% 
% fe_SpringMass = sqrt((k1+k2)/m)/(2*pi());
% fprintf('Erwartetes Stopband bei Eigenfrequenz des MF-Systems: %s \n',fe_SpringMass)
