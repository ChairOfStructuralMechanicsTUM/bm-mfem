    
clear


io=MdpaInput('Federtests_Raute.mdpa'); %specify input file   
model = io.readModel(); %read the model
model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);
% 

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

title('Raute')
%legend(['bandnumbers: ' numberOfBands],'Location','EastOutside')
xlabel('Phase Im(k)')
ylabel('Frequenz f')
xlim([0 pi])
ylim([0 1.5e4])


%