% generated model
clear

%NODES:

n=10;   % n= #Number of Nodes each row/col

   %z=0;
   x=0:2:2*n-2;  % n Coordinates each row
   y=0:2:2*n-2;  % n Coordinates each colum
   
nodeArray=[];

id=1;
for row=1:length(y)
    for col=1:length(x)
        node=Node(id,x(col),y(row));
        nodeArray=[nodeArray node];
        id=id+1;
    end
end

% mat = Material('test');
% mat.addParameter('YOUNGS_MODULUS', 1000);
nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
elementArray=[];
id=1;

for i=1:(n^2-n)
    
    if mod(i,n)~= 0

% Horizontal
         element=BarElement2d2n(id,[nodeArray(i) nodeArray(i+1)]);
         elementArray=[elementArray element];
         id=id+1;
    
% Diagonal  
        element=BarElement2d2n(id,[nodeArray(i) nodeArray(i+(n+1))]);
        id=id+1;
        elementArray=[elementArray element];
    end
    
% Vertical
   element=BarElement2d2n(id,[nodeArray(i) nodeArray(i+n)]);
   elementArray=[elementArray element];
   id=id+1;
end

id=length(elementArray);
% add top line of horizonatl elements
for i=(n^2-n+1):(n^2-1)   
     id=id+1;
     element=BarElement2d2n(id,[nodeArray(i) nodeArray(i+1)]);
     elementArray=[elementArray element];
end
elementArray.setPropertyValue('CROSS_SECTION',1);
elementArray.setPropertyValue('YOUNGS_MODULUS',1000);

%arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

addPointLoad(nodeArray(3),10,[-1 0]);

nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(n^2-n+1).fixDof('DISPLACEMENT_X');
nodeArray(n^2-n+1).fixDof('DISPLACEMENT_Y');

model = FemModel(nodeArray, elementArray);
