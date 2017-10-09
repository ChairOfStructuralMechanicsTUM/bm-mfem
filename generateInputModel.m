% generated model
clear

%NODES:

n=10;   % n= #Number of Nodes each row/col

   z=0;
   x=0:2:2*n-2;  % n Coordinates each row
   y=0:2:2*n-2;  % n Coordinates each colum
   
nodeArray=[];

id=1;
for row=1:length(y)
    for col=1:length(x)
        node=Node(id,x(col),y(row),z);
        nodeArray=[nodeArray node];
        id=id+1;
    end
end

mat = Material('test');
mat.addParameter('YOUNGS_MODULUS', 1000);

elementArray=[];
id=1;

for i=1:(n^2-n)
    
    if mod(i,n)~= 0

% Horizontal
         element=BarElement3d2n(id,[nodeArray(i) nodeArray(i+1)], mat, 5);
         elementArray=[elementArray element];
         id=id+1;
    
% Diagonal  
        element=BarElement3d2n(id,[nodeArray(i) nodeArray(i+(n+1))], mat, 5);
        id=id+1;
        elementArray=[elementArray element];
    end
    
% Vertical
   element=BarElement3d2n(id,[nodeArray(i) nodeArray(i+n)], mat, 5);
   elementArray=[elementArray element];
   id=id+1;
end

id=length(elementArray);
% add top line of horizonatl elements
for i=(n^2-n+1):(n^2-1)   
     id=id+1;
     element=BarElement3d2n(id,[nodeArray(i) nodeArray(i+1)], mat, 5);
     elementArray=[elementArray element];
end


arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

addPointLoad(nodeArray(3:3:60),10,[-1 0 0]);

nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(n^2-n+1).fixDof('DISPLACEMENT_X');
nodeArray(n^2-n+1).fixDof('DISPLACEMENT_Y');

model = FemModel(nodeArray, elementArray);

% Compute:
[Substructure]=createSubstructure(model,1,n);  
 Substructure01=Substructure(1);
 Substructure02=Substructure(2);
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), Substructure02.nodeArray);
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), Substructure02.nodeArray)

tic
[u]=CallPCG(Substructure01,Substructure02);
FetiTime=toc;

tic
 x = SimpleSolvingStrategy.solve(model);
FemTime=toc; 

 