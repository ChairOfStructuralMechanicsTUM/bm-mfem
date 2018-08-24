clear;

node01=Node(1,0,0);
node02=Node(2,5,0);
node03=Node(3,10,0);
node04=Node(4,15,0);
node05=Node(5,20,0);
%
node06=Node(6,0,5);
node07=Node(7,5,5);
node08=Node(8,10,5);
node09=Node(9,15,5);
node10=Node(10,20,5);
%
node11=Node(11,0,10);
node12=Node(12,5,10);
node13=Node(13,10,10);
node14=Node(14,15,10);
node15=Node(15,20,10);
%
node16=Node(16,0,15);
node17=Node(17,5,15);
node18=Node(18,10,15);
node19=Node(19,15,15);
node20=Node(20,20,15);
%
node21=Node(21,0,20);
node22=Node(22,5,20);
node23=Node(23,10,20);
node24=Node(24,15,20);
node25=Node(25,20,20);
%
node26=Node(26,0,25);
node27=Node(27,5,25);
node28=Node(28,10,25);
node29=Node(29,15,25);
node30=Node(30,20,25);
%
node31=Node(31,0,30);
node32=Node(32,5,30);
node33=Node(33,10,30);
node34=Node(34,15,30);
node35=Node(35,20,30);

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
%


ele32 = BarElement2d2n(32,[node16 node21]);
ele33 = BarElement2d2n(33,[node17 node21]);
ele34 = BarElement2d2n(34,[node17 node22]);
ele35 = BarElement2d2n(35,[node18 node22]);
ele36 = BarElement2d2n(36,[node18 node23]);
ele37 = BarElement2d2n(37,[node19 node23]);
ele38 = BarElement2d2n(38,[node19 node24]);
ele39 = BarElement2d2n(39,[node20 node24]);
ele40 = BarElement2d2n(40,[node20 node25]);
ele41 = BarElement2d2n(41,[node10 node14]);

