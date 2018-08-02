clear;
%3x4 nodes
node01 = Node(1,0,0);
node02 = Node(2,0.033,0);
node03 = Node(3,0.066,0);
node04 = Node(4,0.1,0);
node05 = Node(5,0,0.01);
node06 = Node(6,0.033,0.01);
node07 = Node(7,0.066,0.01);
node08 = Node(8,0.1,0.01);
node09 = Node(9,0,0.02);
node10 = Node(10,0.033,0.02);
node11 = Node(11,0.066,0.02);
node12 = Node(12,0.1,0.02);

nodeArray = [node01 node02 node03 node04 node05 ...
    node06 node07 node08 node09 node10 node11 node12];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});


ele01 = QuadrilateralElement2d4n(1,[node01 node02 node06 node05]);
ele02 = QuadrilateralElement2d4n(2,[node02 node03 node07 node06]);
ele03 = QuadrilateralElement2d4n(3,[node03 node04 node08 node07]);
ele04 = QuadrilateralElement2d4n(4,[node05 node06 node10 node09]);
ele05 = QuadrilateralElement2d4n(5,[node06 node07 node11 node10]);
ele06 = QuadrilateralElement2d4n(6,[node07 node08 node12 node11]);




elementArray = [ele01 ele02 ele03 ele04 ele05 ele06];

elementArray.setPropertyValue('YOUNGS_MODULUS',78000000);
elementArray.setPropertyValue('POISSON_RATIO',0.3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);
elementArray.setPropertyValue('DENSITY',70000);


model = FemModel(nodeArray, elementArray);

solver = BlochInverse1D(model);
assembling = SimpleAssembler(model);

stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);


initialize(solver)


% [Ksorted,Msorted] = sortKandM(solver,stiffnessMatrix,massMatrix);
% 
% [Kred,Mred] = reducedStiffnesAndMass (solver,Ksorted,Msorted);  
% 
% 
% omega = cell(10,1);
% f=zeros(1,10);
% f_2=zeros(1,10);
% % f_3=zeros(1,10);
% % f_4=zeros(1,10);
% 
% for i = 1:10
%     omega{i,1} = solver.calcOmega(Kred{i,1},Mred{i,1});
%     f(i) = omega{i,1}(1,1)/(2*pi);
%     f_2(i) = omega{i,1}(2,1)/(2*pi);
% %     f_3(i) = omega{i,1}(3,1)/(2*pi);
% %     f_4(i) = omega{i,1}(4,1)/(2*pi);
%    
% end
% 
% 
% [kx,miu] = propConst(solver,10);
% Bandnummer=1;
% 
% figure(1);
% plot(kx,f)
% title('Dispersion curves')
% xlabel('Phase k')
% ylabel('frequenzy f')
% xlim([0 pi])
% legend(['bandnumber=' num2str(Bandnummer)],'Location','EastOutside')
% hold on
% plot(kx,f_2)


[Ksorted,Msorted] = sortKandM(solver,stiffnessMatrix,massMatrix);

numberOfPhases = 10;

[Kred,Mred] = reducedStiffnesAndMass (solver,Ksorted,Msorted,numberOfPhases);  

omega = cell(numberOfPhases,1);

nob = 3;
[kx,miu] = propConst(solver,numberOfPhases); %already used in reducedStiffnesAndMass(..)

    
for j = 1:nob
%     f(j,1) = zeros(1,numberOfPhases); 
    for i = 1:numberOfPhases
        omega{i,1} = solver.calcOmega(Kred{i,1},Mred{i,1},nob);
        f(j,i) = omega{i,1}(j,1)/(2*pi);

    end
%     if j == 1
%         figure(1)
%         plot(kx,f(j,:))
%         title('Dispersion curves')
%         xlabel('Phase k')
%         ylabel('frequenzy f')
%         xlim([0 pi])
%         legend('1stBand','Location','EastOutside')
%         hold on
%     elseif j == 2
%         plot(kx,f(j,:))
%         legend({'1stBand','2ndBand'},'Location','EastOutside')        
%     else
%       figure(j-1)
        figure(1)
        plot(kx,f(j,:),'r')
        title('Dispersion curves')
        xlabel('Phase k')
        ylabel('frequenzy f')
        xlim([0 pi])
        legend(['bandnumber: ' num2str(j)],'Location','EastOutside')
        hold on
    

end

v = Visualization(model);
v.setScaling(1);
v.plotUndeformed