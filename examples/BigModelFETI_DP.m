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

addPointLoad(nodeArray(100),100,[0 -1]);

nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(n^2-n+1).fixDof('DISPLACEMENT_X');
nodeArray(n^2-n+1).fixDof('DISPLACEMENT_Y');

model = FemModel(nodeArray, elementArray);
assembling = SimpleAssembler2(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
[forceVector,reducedForceVector] = assembling.applyExternalForces(model);
%% testcase substructuring
dim=[n,n];
Ns=16;
v=4;
hz=4;
nodematrixtest=substructureFETI_DP.setupNodeMatrix(model,dim);
[K,bc,br,in,gbc,gbr,gin]=substructureFETI_DP.substructureNodeMatrix(nodematrixtest,Ns,v,hz,dim);
%[DoubleNodes]=substructureFETI_DP.getDoubleNodes(gbr);
[sNodeIdArray]=substructureFETI_DP.getSubstructureNodeIdArray(K,v,hz);
[sNodeArray] = substructureFETI_DP.getSubstructureNodeArray(model,sNodeIdArray,v,hz);
[sElementArray,id,nodes] = substructureFETI_DP.getSubstructureElementArray(model,sNodeIdArray,v,hz);
[sDofArray]= substructureFETI_DP.getSubstrucureDofArray(sNodeArray,v,hz);
[gstiffnessMatrix, greducedStiffnessMatrix] = substructureFETI_DP.assembleSubstructureStiffnessMatrix(model,sElementArray,sDofArray,gbr,v,hz);
[Ksort,Krr,Kcc,Krc,Kcr,suDofId,srDofId,suDofIdLoc,srDofIdLoc]=substructureFETI_DP.splitMatrix(model,gstiffnessMatrix,sDofArray,v,hz,in,br,bc);
[sForceVector,ubcId]=substructureFETI_DP.getSubstructureForceVector(model,assembling,suDofId,gbc,gbr,v,hz);
[gfr,gfbc]=substructureFETI_DP.sortSubstructureForceVector(sForceVector,srDofId,v,hz);
[Bbr,ur,ur2,sinDofId,sbrDofId]=substructureFETI_DP.getInterfaceBooleanMatrix(model,in,sDofArray,srDofId,v,hz);
[Bc,bcgl,bcdofId]=substructureFETI_DP.getCornerBooleanMatrix(model,sDofArray,bc,gbc,hz,v);
[FIrr,FIrc,Kcc,Kccg,dr,fcg]=substructureFETI_DP.assembleAllParameters(v,hz,Kcc,Krc,Krr,Bc,Bbr,gfr,gfbc);
[Kbrbr]=substructureFETI_DP.getBoundryReminderMatrix(Krr,sinDofId,v,hz);
[lP,A]=substructureFETI_DP.getLumpedPreconditioner(Bbr,Kbrbr,sinDofId,srDofId,ur2,v,hz);
[lPglobal]=substructureFETI_DP.assembleLumpedPreconditioner(lP,v,hz);
tic
[lmd]=FETI_DPSolver.PCG(FIrr,FIrc,Kccg,dr,fcg,1);  %Preconditioner funktioniert in diesem Bsp nicht--> wird zu 1 gesetzt
[uc]=FETI_DPSolver.solveCornerDofs(Kccg,fcg,FIrc,lmd);
[urem]=FETI_DPSolver.solveReminderDofs(Krr,gfr,Krc,Bc,uc,Bbr,lmd,v,hz);
[ufinal,ufinalred]=FETI_DPSolver.getResultVector(model,uc,urem,ur,ur2,bcdofId,v,hz);
t=toc
SimpleAssembler.assignResultsToDofs(model, ufinalred);
%%



v = Visualization(model);
%v.setScaling(10000000);
f1=figure(1);
%v.plotUndeformed
v.plotDeformed
xlabel('x')
ylabel('y')