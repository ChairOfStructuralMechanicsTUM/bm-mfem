    
clear


io=MdpaInput('Federtests_Raute.mdpa'); %specify input file   
model = io.readModel(); %read the model
model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);

% model.getModelPart('GENERIC_fixedNodes').getNodes.fixDof('DISPLACEMENT_Y');
% model.getModelPart('GENERIC_fixedNodes').getNodes.fixDof('DISPLACEMENT_X');
% % 
a=model.getAllModelParts;

allNodes = model.getAllNodes();
massNodeID = length(allNodes)+1;

allElements = model.getAllElements();
massID = length(allElements)+1;
spring1ID = massID+1;
spring2ID = massID+2;
spring3ID = massID+3;
spring4ID = massID+4;


springNodes = model.getModelPart('GENERIC_45Grad').getNodes();
leftSpringNode = springNodes(1,1);
leftSNId = getId(leftSpringNode);
rightSpringNode = springNodes(1,2);
rightSNId = getId(rightSpringNode);

springNodes2 = model.getModelPart('GENERIC_135Grad').getNodes();
upperSpringNode = springNodes2(1,1);
upperSNId = getId(upperSpringNode);
lowerSpringNode = springNodes2(1,2);
lowerSNId = getId(lowerSpringNode);


% %% 2 springs, 1 mass
model.addNewNode(massNodeID,0.1,0.075,0);
model.addNewElement('SpringDamperElement3d2n',spring1ID,[leftSNId massNodeID]);
model.addNewElement('SpringDamperElement3d2n',spring2ID,[massNodeID rightSNId]);
model.addNewElement('SpringDamperElement3d2n',spring3ID,[massNodeID upperSNId]);
model.addNewElement('SpringDamperElement3d2n',spring4ID,[massNodeID lowerSNId]);
model.addNewElement('ConcentratedMassElement3d1n',massID, massNodeID);

model.getNode(leftSNId).addDof("DISPLACEMENT_Z");
model.getNode(rightSNId).addDof("DISPLACEMENT_Z");
model.getNode(upperSNId).addDof("DISPLACEMENT_Z");
model.getNode(lowerSNId).addDof("DISPLACEMENT_Z");
model.getNode(massNodeID).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"]);


model.getNode(leftSNId).fixDof('DISPLACEMENT_Z');
model.getNode(rightSNId).fixDof('DISPLACEMENT_Z');
model.getNode(upperSNId).fixDof('DISPLACEMENT_Z');
model.getNode(lowerSNId).fixDof('DISPLACEMENT_Z');

% % % model.getNode(massNodeID).fixDof('DISPLACEMENT_Y');
model.getNode(massNodeID).fixDof('DISPLACEMENT_Z');


model.getElement(spring1ID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
model.getElement(spring1ID).setPropertyValue('ELEMENTAL_DAMPING',0);
model.getElement(spring2ID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
model.getElement(spring2ID).setPropertyValue('ELEMENTAL_DAMPING',0);
model.getElement(spring3ID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
model.getElement(spring3ID).setPropertyValue('ELEMENTAL_DAMPING',0);
model.getElement(spring4ID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
model.getElement(spring4ID).setPropertyValue('ELEMENTAL_DAMPING',0);

model.getElement(massID).setPropertyValue('ELEMENTAL_MASS',7e0);
model.getElement(massID).setPropertyValue('VOLUME_ACCELERATION',10);

% % % Fixing Displacement Y

% % fixingNode = model.getModelPart('GENERIC_fixedNodes').getNodes();
% % fixingNodeID = getId(fixingNode);
% % fixingSpringID = massID + 3;
% % 
% % model.getNode(fixingNodeID).addDof("DISPLACEMENT_Z");
% % 
% % model.getNode(fixingNodeID).fixDof('DISPLACEMENT_X');
% % model.getNode(fixingNodeID).fixDof('DISPLACEMENT_Y');
% % model.getNode(fixingNodeID).fixDof('DISPLACEMENT_Z');
% % 
% % model.addNewElement('SpringDamperElement3d2n',fixingSpringID,[fixingNodeID massNodeID]);
% % 
% % model.getElement(fixingSpringID).setPropertyValue('ELEMENTAL_STIFFNESS',10e20);
% % model.getElement(fixingSpringID).setPropertyValue('ELEMENTAL_DAMPING',0);

v=Visualization(model); %set up visualization
v.plotUndeformed()  %visualize
 
assembling = SimpleAssembler(model);

%% Up to here everything is done as it would be set up of Standard FEM

% Create Bloch Solver
solver = BlochInverse1D_mm(model);

% define number of phases and number of bands
numberOfPhases = 50;
numberOfBands = 10;

% call the solve function of the solver
% phases: contains the discret values of the phase
% frequencies: contains the solution, each row of represents one band
[phases,frequencies]=solver.solve(numberOfPhases,numberOfBands);



%% Visualize results
figure()

for i=1:numberOfBands
    plot(phases,frequencies(i,:),'r')
    hold on
end

title('Raute 45und135Grad - Doppelfeder')
%legend(['bandnumbers: ' numberOfBands],'Location','EastOutside')
xlabel('Phase k')
ylabel('frequenzy f')
xlim([0 pi])
ylim([0 2e4])


% 
% 
% m = model.getElement(massID).getPropertyValue('ELEMENTAL_MASS');
% k1 = model.getElement(spring1ID).getPropertyValue('ELEMENTAL_STIFFNESS');
% k2 = model.getElement(spring2ID).getPropertyValue('ELEMENTAL_STIFFNESS');
% 
% fe_SpringMass = sqrt((k1+k2)/m)/(2*pi());
% fprintf('Erwartetes Stopband bei Eigenfrequenz des MF-Systems: %s \n',fe_SpringMass)