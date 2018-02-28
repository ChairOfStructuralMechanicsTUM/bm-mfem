clear;

node01 = Node(1,0,0);
node02 = Node(2,2,0);
node03 = Node(3,2,1);
node04 = Node(4,0,1);

nodeArray = [node01 node02 node03 node04];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

ele01 = QuadrilateralElement2d4n(1,[node01 node02 node03 node04]);

elementArray = [ele01];

elementArray.setPropertyValue('YOUNGS_MODULUS',96000000);
elementArray.setPropertyValue('POISSON_RATIO',1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);
elementArray.setPropertyValue('DENSITY',7860);

elementIds = elementArray.getId;

node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');
node04.fixDof('DISPLACEMENT_X');
node04.fixDof('DISPLACEMENT_Y');

addPointLoad(node03,1,[0 -1]);

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);

for i=1:8
    for j=1:8
        if stiffnessMatrix(i,j)<1e-6
            stiffnessMatrix(i,j)=0;
        end
    end
end
            
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleSolvingStrategy(model);
x = solver.solve();

step = 1;

% stressVector = computeElementStress(elementArray, step)';


VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);

%x=0:0.25:1.75;
fig=figure;
plot(0:1:1.75,VerschiebungDofs(2:2:4))
    
    
