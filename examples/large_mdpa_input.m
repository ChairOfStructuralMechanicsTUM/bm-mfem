% LARGE_MDPA_INPUT Read a large mdpa model and plot it
clear
close all

io=MdpaInput('large_mdpa_input.mdpa'); %specify input file
model = io.readModel(); %read the model

model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y", "DISPLACEMENT_Z", ...
    "ROTATION_X", "ROTATION_Y",  "ROTATION_Z"]);


%set boundary conditions for all elements in 'left_support' and 'right_support':
model.getModelPart('GENERIC_left_support').getNodes().fixAllDofs();
model.getModelPart('GENERIC_right_support').getNodes().fixAllDofs();

ln = model.getModelPart('GENERIC_left_support').getNodes()
rn = model.getModelPart('GENERIC_right_support').getNodes()
%set load for all elements in 'inner_circle':
% model.getModelPart('GENERIC_inner_circle').getNodes().setDofLoad('DISPLACEMENT_X',10000);


v=Visualization(model); %set up visualization
v.plotUndeformed()  %visualize


model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);


solver = BlochInverse1D(model);
assembling = SimpleAssembler(model);

stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);


initialize(solver)


[Ksorted,Msorted] = sortKandM(solver,stiffnessMatrix,massMatrix);

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

    figure(1)
    plot(kx,f(j,:),'r')
    title('Dispersion curves')
    xlabel('Phase k')
    ylabel('frequenzy f')
    xlim([0 pi])
    legend(['bandnumber: ' num2str(j)],'Location','EastOutside')
    hold on
    

end