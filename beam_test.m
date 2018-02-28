clear;

node01 = Node(1,0,0);
node02 = Node(2,0.25,0);
node03 = Node(3,0.5,0);
node04 = Node(4,0.75,0);
node05 = Node(5,1,0);
node06 = Node(6,1.25,0);
node07 = Node(7,1.5,0);
node08 = Node(8,1.75,0);
node09 = Node(9,0,0.01);
node10 = Node(10,0.25,0.01);
node11 = Node(11,0.5,0.01);
node12 = Node(12,0.75,0.01);
node13 = Node(13,1,0.01);
node14 = Node(14,1.25,0.01);
node15 = Node(15,1.5,0.01);
node16 = Node(16,1.75,0.01);


nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08 node09 node10 node11 node12 node13 node14 node15 node16];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

ele01 = QuadrilateralElement2d4n(1,[node01 node02 node10 node09]);
ele02 = QuadrilateralElement2d4n(2,[node02 node03 node11 node10]);
ele03 = QuadrilateralElement2d4n(3,[node03 node04 node12 node11]);
ele04 = QuadrilateralElement2d4n(4,[node04 node05 node13 node12]);
ele05 = QuadrilateralElement2d4n(5,[node05 node06 node14 node13]);
ele06 = QuadrilateralElement2d4n(6,[node06 node07 node15 node14]);
ele07 = QuadrilateralElement2d4n(7,[node07 node08 node16 node15]);

elementArray = [ele01 ele02 ele03 ele04 ele05 ele06 ele07];

elementArray.setPropertyValue('YOUNGS_MODULUS',210e9);
elementArray.setPropertyValue('POISSON_RATIO',0.3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);
elementArray.setPropertyValue('DENSITY',7860);

elementIds = elementArray.getId;

node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');
node09.fixDof('DISPLACEMENT_X');
node09.fixDof('DISPLACEMENT_Y');

addPointLoad(node08,1,[0 -1]);

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);

for i=1:32
    for j=1:32
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
plot(0:0.25:1.75,VerschiebungDofs(2:2:16))
    
    
