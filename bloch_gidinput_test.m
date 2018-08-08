
clear

io=MdpaInput('Cell_Input.mdpa'); %specify input file
model = io.readModel(); %read the model

model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);

innerNodes = model.getModelPart('GENERIC_Inner_Points').getNodes();
id = getId(innerNodes) ;
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

% model.addNewNode(244,0.1,0.075,0);
% model.addNewElement('SpringDamperElement3d2n',188,[59 244]);
% model.addNewElement('SpringDamperElement3d2n',189,[244 159]);
% model.addNewElement('ConcentratedMassElement3d1n',190, 244);
% 
% model.getNode(59).addDof("DISPLACEMENT_Z");
% model.getNode(159).addDof("DISPLACEMENT_Z");
% model.getNode(244).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"]);
% 
% model.getNode(59).fixDof('DISPLACEMENT_Z');
% model.getNode(159).fixDof('DISPLACEMENT_Z');
% model.getNode(244).fixDof('DISPLACEMENT_Y');
% model.getNode(244).fixDof('DISPLACEMENT_Z');
% 
% model.getElement(188).setPropertyValue('ELEMENTAL_STIFFNESS',7e5);
% model.getElement(188).setPropertyValue('ELEMENTAL_DAMPING',0);
% model.getElement(189).setPropertyValue('ELEMENTAL_STIFFNESS',7e5);
% model.getElement(189).setPropertyValue('ELEMENTAL_DAMPING',0);
% 
% model.getElement(190).setPropertyValue('ELEMENTAL_MASS',7e2);
% model.getElement(190).setPropertyValue('VOLUME_ACCELERATION',10);

% solver = BlochInverse1D(model);
% assembling = SimpleAssembler(model);
% 
% [stiffnessMatrix,Kred1] = assembling.assembleGlobalStiffnessMatrix(model);            
% [massMatrix,Mred1] = assembling.assembleGlobalMassMatrix(model);
% 
% Dofs = model.getDofArray();
% 
% 
% 
% initialize(solver)
% % [Ksorted,Msorted] = sortKandM(solver,stiffnessMatrix,massMatrix);
% [Ksorted,Msorted] = sortKandM(solver,Kred1,Mred1);
% numberOfPhases = 20;
% [Kred,Mred] = reducedStiffnesAndMass (solver,Ksorted,Msorted,numberOfPhases);  
% omega = cell(numberOfPhases,1);
% nob = 10;
% [kx,miu] = propConst(solver,numberOfPhases); %already used in reducedStiffnesAndMass(..)
%     
% for j = 1:nob
%     for i = 1:numberOfPhases
%         omega{i,1} = solver.calcOmega(Kred{i,1},Mred{i,1},nob);
%         f(j,i) = omega{i,1}(j,1)/(2*pi);
%     end
%     figure(1)
%     plot(kx,f(j,:),'r')
%     title('Dispersion curves')
%     xlabel('Phase k')
%     ylabel('frequenzy f')
%     xlim([0 pi])
%     legend(['bandnumber: ' num2str(j)],'Location','EastOutside')
%     hold on
% end

v=Visualization(model); %set up visualization
v.plotUndeformed()  %visualize