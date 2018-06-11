clear;

node01 = Node(1,-1,-1,-1);
node02 = Node(2,1,-1,-1);
node03 = Node(3,1,1,-1);
node04 = Node(4,-1,1,-1);
node05 = Node(5,-1,-1,1);
node06 = Node(6,1,-1,1);
node07 = Node(7,1,1,1);
node08 = Node(8,-1,1,1);
node09 = Node(9,-1,-1,2);
node10 = Node(10,1,-1,2);
node11 = Node(11,1,1,2);
node12 = Node(12,-1,1,2);

nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08 node09 node10 node11 node12];

nodeArray.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_SOLID_Z','DISPLACEMENT_FLUID_X', 'DISPLACEMENT_FLUID_Y', 'DISPLACEMENT_FLUID_Z'});

ele01 = PorousElement3d8n(1,[node01 node02 node03 node04 node05 node06 node07 node08]);
ele02 = PorousElement3d8n(2,[node05 node06 node07 node08 node09 node10 node11 node12]);

elementArray = [ele01 ele02];

elementArray.setPropertyValue('DENSITY_SOLID',30);
elementArray.setPropertyValue('LAMBDA_SOLID',905357);
elementArray.setPropertyValue('MUE_SOLID',264062);
elementArray.setPropertyValue('DAMPING_SOLID',0);

elementArray.setPropertyValue('DENSITY_FLUID',1.21);
elementArray.setPropertyValue('VISCOSITY_FLUID',1.84e-5);
elementArray.setPropertyValue('STANDARD_PRESSURE_FLUID',101);
elementArray.setPropertyValue('HEAT_CAPACITY_FLUID',1.4);
elementArray.setPropertyValue('PRANDTL_NUMBER_FLUID',0.71);

elementArray.setPropertyValue('POROSITY',0.96);
elementArray.setPropertyValue('TORTUOSITY',1.7);
elementArray.setPropertyValue('FLOW_RESISTIVITY',32e3);
elementArray.setPropertyValue('VISCOUS_LENGHT',90);
elementArray.setPropertyValue('THERMAL_LENGTH',165);

elementArray.setPropertyValue('FREQUENCY',9);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);

node01.fixDof('DISPLACEMENT_SOLID_X');
node01.fixDof('DISPLACEMENT_SOLID_Y');
node01.fixDof('DISPLACEMENT_SOLID_Z');
node04.fixDof('DISPLACEMENT_SOLID_X');
node04.fixDof('DISPLACEMENT_SOLID_Y');
node04.fixDof('DISPLACEMENT_SOLID_Z');
node05.fixDof('DISPLACEMENT_SOLID_X');
node05.fixDof('DISPLACEMENT_SOLID_Y');
node05.fixDof('DISPLACEMENT_SOLID_Z');
node08.fixDof('DISPLACEMENT_SOLID_X');
node08.fixDof('DISPLACEMENT_SOLID_Y');
node08.fixDof('DISPLACEMENT_SOLID_Z');

node01.fixDof('DISPLACEMENT_FLUID_X');
node01.fixDof('DISPLACEMENT_FLUID_Y');
node01.fixDof('DISPLACEMENT_FLUID_Z');
node04.fixDof('DISPLACEMENT_FLUID_X');
node04.fixDof('DISPLACEMENT_FLUID_Y');
node04.fixDof('DISPLACEMENT_FLUID_Z');
node05.fixDof('DISPLACEMENT_FLUID_X');
node05.fixDof('DISPLACEMENT_FLUID_Y');
node05.fixDof('DISPLACEMENT_FLUID_Z');
node08.fixDof('DISPLACEMENT_FLUID_X');
node08.fixDof('DISPLACEMENT_FLUID_Y');
node08.fixDof('DISPLACEMENT_FLUID_Z');

addPointLoadPorous(node02,1,[1 0 0]);
addPointLoadPorous(node03,2,[0 -1 0]);
addPointLoadPorous(node07,1,[0 0 -1]);

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);

% solver = SimpleSolvingStrategy(model);
% x = solver.solve();
% 
% step = 1;
% 
% VerschiebungDofs = model.getDofArray.getValue(step);
% 
% nodalForces = solver.getNodalForces(step);
% 
v = Visualization(model);
v.setScaling(5000);
v.plotUndeformed
% v.plotDeformed

%% Dynamic Test

% mid = 13;
% 
% Solver
% model.getNode(mid).setDofLoad('DISPLACEMENT_Z', 1);
% 
% dt = .01;
% time = 0;
% endTime = .2;
% solver = NewmarkSolvingStrategy(model, dt);
% 
% while time < endTime
%     solver.solve();
%     time = time + dt;
% end

dt = .05;
time = 0;
endTime = 1.5;
ls = linspace(time,endTime,endTime/dt+1);
% disp = zeros(1,endTime/dt+1);
solver = NewmarkSolvingStrategy(model, dt);
% disp = 0;
% while time < endTime
%     solver.solve();
%     time = time + dt;
%         hold on
%         v.plotDeformed();
%         hold off
%         pause(0.1)
% %         plot(model.getAllNodes.getDofValue('DISPLACEMENT_SOLID_Z','end'))
% %     
% %         disp(model.getProperties.getValue('STEP')) = node05.getDofValue('DISPLACEMENT_SOLID_Z','end');
%     
% end
% plot(ls,endnode.getDofValue('DISPLACEMENT_Y','all'))
% solver = SimpleSolvingStrategy(model);
% solver.solve();

%% Movie

for j = 1:30
    solver.solve();
    time = j*0.05;
    hold on
    v.plotDeformed();
    view(3)
    axis([-1.1 1.1 -1.1 1.1 -1.1 2.1])
    F(j) = getframe(gcf);
end

fig = figure;
movie(fig,F,1,2)
