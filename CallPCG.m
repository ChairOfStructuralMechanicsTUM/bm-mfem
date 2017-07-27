% Solve Substructure:

function [u1,u2]=CallPCG(Substructure01,Substructure02,interfaceNodes)
StiffnessMatrixSub1=SimpleAssembler(Substructure01);
StiffnessMatrixSub2=SimpleAssembler(Substructure02);

K1=StiffnessMatrixSub1.reducedStiffnessMatrix;
K2=StiffnessMatrixSub2.reducedStiffnessMatrix;
f1=StiffnessMatrixSub1.reducedForceVector';
f2=StiffnessMatrixSub2.reducedForceVector';

    

[R2,KS]=computeRidgedBodyModes(Substructure02);


B1=computeB(Substructure01,interfaceNodes);% Substructure.IntrefaceNodes
interfaceNodes=[1 2];
B2=computeB(Substructure02,interfaceNodes);  

%Gleichungssystem Aufstellen:
F1=B1*K1^-1*B1'+B2*KS*B2';
Gi2=B2*R2;
d=B2*KS*f2-B1*K1^-1*f1 ;  
e=R2'*f2;
Pi=[B1 B2]*[K1, zeros(size(K1,1),size(K2,2));...
    zeros(size(K2,2),size(K1,1)) K2]*[B1';B2'];

[lambda,alpha]= PCG(F1,Pi,Gi2,d,e);
 
%[lambda,alpha]=CG2(F1,Gi2,d,e);

u1=K1\(f1+B1'*lambda);
u2=KS*(f2-B2'*lambda)+R2*alpha;  % KS statt "\" daa KS schon berechnet!!
GG=abs(B1*u1)-abs(B2*u2);
if GG'*GG >10e-10
   disp('Gleichgewicht: B1*u1=B2*u2 nicht eingehalten')
end
 
%    U1=K1^-1*(f1+B1'*x(1:4));
%    U2=KS*(f2-B2'*x(1:4))+R2*x(5:7);
%    % GG=B1*U1-B2*U2;

end