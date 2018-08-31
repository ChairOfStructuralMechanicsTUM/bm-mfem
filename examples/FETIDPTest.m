clear;

node01=Node(1,0,0);
node02=Node(2,0,-5);
node03=Node(3,0,-10);
node04=Node(4,0,-15);
node05=Node(5,0,-20);
%
node06=Node(6,5,0);
node07=Node(7,5,-5);
node08=Node(8,5,-10);
node09=Node(9,5,-15);
node10=Node(10,5,-20);
%
node11=Node(11,10,0);
node12=Node(12,10,-5);
node13=Node(13,10,-10);
node14=Node(14,10,-15);
node15=Node(15,10,-20);
%
node16=Node(16,15,0);
node17=Node(17,15,-5);
node18=Node(18,15,-10);
node19=Node(19,15,-15);
node20=Node(20,15,-20);
%
node21=Node(21,20,0);
node22=Node(22,20,-5);
node23=Node(23,20,-10);
node24=Node(24,20,-15);
node25=Node(25,20,-20);
%
node26=Node(26,25,0);
node27=Node(27,25,-5);
node28=Node(28,25,-10);
node29=Node(29,25,-15);
node30=Node(30,25,-20);
%
node31=Node(31,30,0);
node32=Node(32,30,-5);
node33=Node(33,30,-10);
node34=Node(34,30,-15);
node35=Node(35,30,-20);

nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08 node09 ...
    node10 node11 node12 node13 node14 node15 node16 node17 node18 node19 node20 ...
    node21 node22 node23 node24 node25 node26 node27 node28 node29 node30 ...
    node31 node32 node33 node34 node35];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});


ele01 = BarElement2d2n(1,[node01 node02]);
ele02 = BarElement2d2n(2,[node02 node03]);
ele03 = BarElement2d2n(3,[node03 node04]);
ele04 = BarElement2d2n(4,[node04 node05]);
ele05 = BarElement2d2n(5,[node01 node06]);
ele06 = BarElement2d2n(6,[node02 node06]);
ele07 = BarElement2d2n(7,[node02 node07]);
ele08 = BarElement2d2n(8,[node03 node07]);
ele09 = BarElement2d2n(9,[node03 node08]);
ele10 = BarElement2d2n(10,[node04 node08]);
ele11 = BarElement2d2n(11,[node04 node09]);
ele12 = BarElement2d2n(12,[node05 node09]);
ele13 = BarElement2d2n(13,[node05 node10]);
ele14 = BarElement2d2n(14,[node06 node07]);
ele15 = BarElement2d2n(15,[node07 node08]);
ele16 = BarElement2d2n(16,[node08 node09]);
ele17 = BarElement2d2n(17,[node09 node10]);
ele18 = BarElement2d2n(18,[node06 node11]);
ele19 = BarElement2d2n(19,[node07 node11]);
ele20 = BarElement2d2n(20,[node07 node12]);
ele21 = BarElement2d2n(21,[node08 node12]);
ele22 = BarElement2d2n(22,[node08 node13]);
ele23 = BarElement2d2n(23,[node09 node13]);
ele24 = BarElement2d2n(24,[node09 node14]);
ele25 = BarElement2d2n(25,[node10 node14]);
ele26 = BarElement2d2n(26,[node10 node15]);
ele27 = BarElement2d2n(27,[node11 node12]);
ele28 = BarElement2d2n(28,[node12 node13]);
ele29 = BarElement2d2n(29,[node13 node14]);
ele30 = BarElement2d2n(30,[node14 node15]);
ele31 = BarElement2d2n(31,[node11 node16]);
ele32 = BarElement2d2n(32,[node12 node16]);
ele33 = BarElement2d2n(33,[node12 node17]);
ele34 = BarElement2d2n(34,[node13 node17]);
ele35 = BarElement2d2n(35,[node13 node18]);
ele36 = BarElement2d2n(36,[node14 node18]);
ele37 = BarElement2d2n(37,[node14 node19]);
ele38 = BarElement2d2n(38,[node15 node19]);
ele39 = BarElement2d2n(39,[node15 node20]);
ele40 = BarElement2d2n(40,[node16 node17]);
ele41 = BarElement2d2n(41,[node17 node18]);
ele42 = BarElement2d2n(42,[node18 node19]);
ele43 = BarElement2d2n(43,[node19 node20]);
ele44 = BarElement2d2n(44,[node16 node21]);
ele45 = BarElement2d2n(45,[node17 node21]);
ele46 = BarElement2d2n(46,[node17 node22]);
ele47 = BarElement2d2n(47,[node18 node22]);
ele48 = BarElement2d2n(48,[node18 node23]);
ele49 = BarElement2d2n(49,[node19 node23]);
ele50 = BarElement2d2n(50,[node19 node24]);
ele51 = BarElement2d2n(51,[node20 node24]);
ele52 = BarElement2d2n(52,[node20 node25]);
ele53 = BarElement2d2n(53,[node21 node22]);
ele54 = BarElement2d2n(54,[node22 node23]);
ele55 = BarElement2d2n(55,[node23 node24]);
ele56 = BarElement2d2n(56,[node24 node25]);
ele57 = BarElement2d2n(57,[node21 node26]);
ele58 = BarElement2d2n(58,[node22 node26]);
ele59 = BarElement2d2n(59,[node22 node27]);
ele60 = BarElement2d2n(60,[node23 node27]);
ele61 = BarElement2d2n(61,[node23 node28]);
ele62 = BarElement2d2n(62,[node24 node28]);
ele63 = BarElement2d2n(63,[node24 node29]);
ele64 = BarElement2d2n(64,[node25 node29]);
ele65 = BarElement2d2n(65,[node25 node30]);
ele66 = BarElement2d2n(66,[node26 node27]);
ele67 = BarElement2d2n(67,[node27 node28]);
ele68 = BarElement2d2n(68,[node28 node29]);
ele69 = BarElement2d2n(69,[node29 node30]);
ele70 = BarElement2d2n(70,[node26 node31]);
ele71 = BarElement2d2n(71,[node27 node31]);
ele72 = BarElement2d2n(72,[node27 node32]);
ele73 = BarElement2d2n(73,[node28 node32]);
ele74 = BarElement2d2n(74,[node28 node33]);
ele75 = BarElement2d2n(75,[node29 node33]);
ele76 = BarElement2d2n(76,[node29 node34]);
ele77 = BarElement2d2n(77,[node30 node34]);
ele78 = BarElement2d2n(78,[node30 node35]);
ele79 = BarElement2d2n(79,[node31 node32]);
ele80 = BarElement2d2n(80,[node32 node33]);
ele81 = BarElement2d2n(81,[node33 node34]);
ele82 = BarElement2d2n(82,[node34 node35]);


elementArray = [ele01 ele02 ele03 ele04 ele05 ele06 ele07 ele08 ele09 ...
    ele10 ele11 ele12 ele13 ele14 ele15 ele16 ele17 ele18 ele19 ...
    ele20 ele21 ele22 ele23 ele24 ele25 ele26 ele27 ele28 ele29 ele30 ...
    ele31 ele32 ele33 ele34 ele35 ele36 ele37 ele38 ele39 ele40 ...
    ele41 ele42 ele43 ele44 ele45 ele46 ele47 ele48 ele49 ele50 ...
    ele51 ele52 ele53 ele54 ele55 ele56 ele57 ele58 ele59 ele60 ...
    ele61 ele62 ele63 ele64 ele65 ele66 ele67 ele68 ele69 ele70 ...
    ele71 ele72 ele73 ele74 ele75 ele76 ele77 ele78 ele79 ele80 ...
    ele81 ele82
    ];


elementArray.setPropertyValue('CROSS_SECTION',1);
elementArray.setPropertyValue('YOUNGS_MODULUS',1000);

elementIds = elementArray.getId;

node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');

node02.fixDof('DISPLACEMENT_X');
node02.fixDof('DISPLACEMENT_Y');

node04.fixDof('DISPLACEMENT_X');
node04.fixDof('DISPLACEMENT_Y');

node05.fixDof('DISPLACEMENT_X');
node05.fixDof('DISPLACEMENT_Y');

node35.fixDof('DISPLACEMENT_X');
node35.fixDof('DISPLACEMENT_Y');

addPointLoad([node21 node26 node31],10,[0 -1]);


model = FemModel(nodeArray, elementArray);
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
[forceVector,reducedForceVector] = assembling.applyExternalForces(model);

%% testcase substructuring
dim=[5,7];
Ns=6;
v=2;
hz=3;
nodematrixtest=substructureFETI_DP.setupNodeMatrix(model,dim);
[K,bc,br,in,gbc,gbr,gin]=substructureFETI_DP.substructureNodeMatrix(nodematrixtest,Ns,v,hz,dim);
[DoubleNodes]=substructureFETI_DP.getDoubleNodes(gbr);
[sNodeIdArray]=substructureFETI_DP.getSubstructureNodeIdArray(K,v,hz);
[sNodeArray] = substructureFETI_DP.getSubstructureNodeArray(model,sNodeIdArray,v,hz);
[sElementArray,id,nodes] = substructureFETI_DP.getSubstructureElementArray(model,sNodeIdArray,v,hz);
[sDofArray]= substructureFETI_DP.getSubstrucureDofArray(sNodeArray,v,hz);
[gstiffnessMatrix, greducedStiffnessMatrix] = substructureFETI_DP.assembleSubstructureStiffnessMatrix(model,sElementArray,sDofArray,v,hz);
[Ksort,Krr,Kcc,Krc,Kcr,suDofId,srDofId,suDofIdLoc,srDofIdLoc]=substructureFETI_DP.splitMatrix(model,gstiffnessMatrix,sDofArray,v,hz,in,br,bc);
[sForceVector]=substructureFETI_DP.getSubstructureForceVector(model,assembling,suDofId,v,hz);
[gfr,gfbc]=substructureFETI_DP.sortSubstructureForceVector(sForceVector,srDofId,v,hz);
[Bbr,ur,ur2,sinDofId,sbrDofId]=substructureFETI_DP.getInterfaceBooleanMatrix(model,in,sDofArray,srDofId,v,hz);
[Bc,bcgl,bcdof]=substructureFETI_DP.getCornerBooleanMatrix(model,bc,gbc,hz,v);
[FIrr,FIrc,Kcc,Kccg,dr,fcg]=substructureFETI_DP.assembleAllParameters(v,hz,Kcc,Krc,Krr,Bc,Bbr,gfr,gfbc);
[Kbrbr]=substructureFETI_DP.getBoundryReminderMatrix(Krr,sinDofId,v,hz);
[lP]=substructureFETI_DP.getLumpedPreconditioner(Bbr,Kbrbr,sinDofId,srDofId,ur2,v,hz);
[lPglobal]=substructureFETI_DP.assembleLumpedPreconditioner(lP,v,hz);
[lmd]=FETI_DPSolver.PCG(FIrr,FIrc,Kccg,dr,fcg,lPglobal);
[uc]=FETI_DPSolver.solveCornerDofs(Kccg,fcg,FIrc,lmd);
[urem]=FETI_DPSolver.solveReminderDofs(Krr,gfr,Krc,Bc,uc,Bbr,lmd,v,hz);
[ufinal,ufinalred]=FETI_DPSolver.getResultVector(model,uc,urem,ur,ur2,bcdof,v,hz);
SimpleAssembler.assignResultsToDofs(model, ufinalred);
%%



v = Visualization(model);
%v.setScaling(10000000);
f1=figure(1);
v.plotUndeformed
v.plotDeformed







