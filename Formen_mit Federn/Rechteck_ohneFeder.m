%%%% Kreis ohne Feder
    
clear

io=MdpaInput('Federtests_Rechteck.mdpa'); %specify input file
model = io.readModel(); %read the model
model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);

% % 
% % a=model.getAllModelParts;
% % leftNodes = model.getModelPart('GENERIC_leftNodes').getNodes();
% % rightNodes = model.getModelPart('GENERIC_rightNodes').getNodes();
% % springNodes = model.getModelPart('GENERIC_springNodes').getNodes();
%     innerMaterial = model.getModelPart('GENERIC_innerMaterial').getNodes()
% 
% coords = getCoords(node)
% % springNodeIds = arrayfun(@(node) node.getId, springNodes);
% % springNodeCoords = arrayfun(@(node) node.getCoords, springNodes,'UniformOutput',false);
% % 
% % leftSpringNode = springNodes(1,2);
% % leftSNCoords = getCoords(leftSpringNode);
% % leftSNId = getId(leftSpringNode);
% % rightSpringNode = springNodes(1,7);
% % rightSNCoords = getCoords(rightSpringNode);
% % rightSNId = getId(rightSpringNode);
% % horizontalSpringNodes = [leftSpringNode rightSpringNode];
% % 
% % 
% % 
% % % %% 2 springs, 1 mass
% % model.addNewNode(229,0.1,0.075,0);
% % model.addNewElement('SpringDamperElement3d2n',199,[58 229]);
% % model.addNewElement('SpringDamperElement3d2n',200,[229 114]);
% % model.addNewElement('ConcentratedMassElement3d1n',201, 229);
% % 
% % model.getNode(58).addDof("DISPLACEMENT_Z");
% % model.getNode(114).addDof("DISPLACEMENT_Z");
% % model.getNode(229).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"]);
% % 
% % 
% % model.getNode(58).fixDof('DISPLACEMENT_Z');
% % model.getNode(114).fixDof('DISPLACEMENT_Z');
% % model.getNode(229).fixDof('DISPLACEMENT_Y');
% % % % model.getNode(229).fixDof('DISPLACEMENT_Z');
% % 
% % 
% % model.getElement(199).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
% % model.getElement(199).setPropertyValue('ELEMENTAL_DAMPING',0);
% % model.getElement(200).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
% % model.getElement(200).setPropertyValue('ELEMENTAL_DAMPING',0);
% % 
% % model.getElement(201).setPropertyValue('ELEMENTAL_MASS',7e0);
% % model.getElement(201).setPropertyValue('VOLUME_ACCELERATION',10);

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
title('Dispersion curves - Rechteck oF')
xlabel('Phase k')
ylabel('frequenzy f')
xlim([0 pi])

ylim([0 2e4])

hold on
for j = 1:nob
    
    for i = 1:numberOfPhases
        omega{i,1} = solver.calcOmega(Kred{i,1},Mred{i,1},nob);
        f(j,i) = omega{i,1}(j,1)/(2*pi);

    end
    plot(kx,f(j,:),'r')       
    legend(['bandnumbers: ' num2str(j)],'Location','EastOutside')
end


% m = model.getElement(201).getPropertyValue('ELEMENTAL_MASS');
% k1 = model.getElement(199).getPropertyValue('ELEMENTAL_STIFFNESS');
% k2 = model.getElement(200).getPropertyValue('ELEMENTAL_STIFFNESS');
% 
% fe_SpringMass = sqrt((k1+k2)/m)/(2*pi());
% fprintf('Erwartetes Stopband bei Eigenfrequenz des MF-Systems: %s \n',fe_SpringMass)
