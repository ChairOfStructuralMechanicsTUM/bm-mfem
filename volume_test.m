clear;

node01 = Node(1,-1,-1,-1);
node02 = Node(2,1,-1,-1);
node03 = Node(3,1,1,-1);
node04 = Node(4,-1,1,-1);
node05 = Node(5,-1,-1,1);
node06 = Node(6,1,-1,1);
node07 = Node(7,1,1,1);
node08 = Node(8,-1,1,1);


nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});

ele01 = Hexahedron3d8n(1,[node01 node02 node03 node04 node05 node06 node07 node08]);

elementArray = [ele01];

elementArray.setPropertyValue('YOUNGS_MODULUS',32000);
elementArray.setPropertyValue('POISSON_RATIO',1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',3);
elementArray.setPropertyValue('DENSITY',7860);

elementIds = elementArray.getId;

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);

massMatrix = assembling.assembleGlobalMassMatrix(model);

a=[16,6,6,-8,2,2,-6,-6,1,4,-2,3,4,3,-2,-6,1,-6,-4,-3,-3,0,-1,-1;
    6,16,6,-2,4,3,-6,-6,1,2,-8,2,3,4,-2,-1,0,-1,-3,-4,-3,1,-6,-6;
    6,6,16,-2,3,4,-1,-1,0,3,-2,4,2,2,-8,-6,1,-6,-3,-3,-4,1,-6,-6;
    -8,-2,-2,16,-6,-6,4,2,-3,-6,6,-1,-6,-1,6,4,-3,2,0,1,1,-4,3,3;
    2,4,3,-6,16,6,-2,-8,2,6,-6,1,1,0,-1,-3,4,-2,-1,-6,-6,3,-4,-3;
    2,3,4,-6,6,16,-3,-2,4,1,-1,0,6,1,-6,-2,2,-8,-1,-6,-6,3,-3,-4;
    -6,-6,-1,4,-2,-3,16,6,-6,-8,2,-2,-4,-3,3,0,-1,1,4,3,2,-6,1,6;
    -6,-6,-1,2,-8,-2,6,16,-6,-2,4,-3,-3,-4,3,1,-6,6,3,4,2,-1,0,1;
    1,1,0,-3,2,4,-6,-6,16,2,-3,4,3,3,-4,-1,6,-6,-2,-2,-8,6,-1,-6;
    4,2,3,-6,6,1,-8,-2,2,16,-6,6,0,1,-1,-4,3,-3,-6,-1,-6,4,-3,-2;
    -2,-8,-2,6,-6,-1,2,4,-3,-6,16,-6,-1,-6,6,3,-4,3,1,0,1,-3,4,2;
    3,2,4,-1,1,0,-2,-3,4,6,-6,16,1,6,-6,-3,3,-4,-6,-1,-6,2,-2,-8;
    4,3,2,-6,1,6,-4,-3,3,0,-1,1,16,6,-6,-8,2,-2,-6,-6,-1,4,-2,-3;
    3,4,2,-1,0,1,-3,-4,3,1,-6,6,6,16,-6,-2,4,-3,-6,-6,-1,2,-8,-2;
    -2,-2,-8,6,-1,-6,3,3,-4,-1,6,-6,-6,-6,16,2,-3,4,1,1,0,-3,2,4;
    -6,-1,-6,4,-3,-2,0,1,-1,-4,3,-3,-8,-2,2,16,-6,6,4,2,3,-6,6,1;
    1,0,1,-3,4,2,-1,-6,6,3,-4,3,2,4,-3,-6,16,-6,-2,-8,-2,6,-6,-1;
    -6,-1,-6,2,-2,-8,1,6,-6,-3,3,-4,-2,-3,4,6,-6,16,3,2,4,-1,1,0;
    -4,-3,-3,0,-1,-1,4,3,-2,-6,1,-6,-6,-6,1,4,-2,3,16,6,6,-8,2,2;
    -3,-4,-3,1,-6,-6,3,4,-2,-1,0,-1,-6,-6,1,2,-8,2,6,16,6,-2,4,3;
    -3,-3,-4,1,-6,-6,2,2,-8,-6,1,-6,-1,-1,0,3,-2,4,6,6,16,-2,3,4;
    0,1,1,-4,3,3,-6,-1,6,4,-3,2,4,2,-3,-6,6,-1,-8,-2,-2,16,-6,-6;
    -1,-6,-6,3,-4,-3,1,0,-1,-3,4,-2,-2,-8,2,6,-6,1,2,4,3,-6,16,6;
    -1,-6,-6,3,-3,-4,6,1,-6,-2,2,-8,-3,-2,4,1,-1,0,2,3,4,-6,6,16];

Fehler=stiffnessMatrix-a;

for i=1:24
    for j=1:24
        if abs(Fehler(i,j))<1e-5
            Fehler(i,j)=0;
        end
    end
end
    
    
