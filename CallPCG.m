% Solve Substructure:

function [u]=CallPCG(Substructure01,Substructure02)
StiffnessMatrixSub1=SimpleAssembler(Substructure01);
StiffnessMatrixSub2=SimpleAssembler(Substructure02);

K1=StiffnessMatrixSub1.reducedStiffnessMatrix;
K2=StiffnessMatrixSub2.reducedStiffnessMatrix;
f1=StiffnessMatrixSub1.reducedForceVector';
f2=StiffnessMatrixSub2.reducedForceVector';

    

[R2,KS]=computeRidgedBodyModes(Substructure02);

[B1]=computeB(Substructure01);% Substructure.IntrefaceNodes
[B2]=computeB(Substructure02);  

%Gleichungssystem Aufstellen:
F1=B1*K1^-1*B1'+B2*KS*B2';
Gi2=B2*R2;
d=B2*KS*f2-B1*K1^-1*f1; 
e=R2'*f2;

%Pi=[B1 B2]*[K1, zeros(size(K1,1),size(K2,2));...
% zeros(size(K2,2),size(K1,1)) K2]*[B1';B2'];

    % matlab optimazition function
%[lambda,alpha]=minimize(F1,Gi2,d,e)
%alpha=(Gi2*((Gi2'*Gi2)^-1))'*(d-F1*x);
a=1;
switch a
    case 1
[lambda,alpha]=CG(F1,Gi2,d,e);
    case 2
[lambda,alpha]= CG2(F1,Gi2,d,e); %CG - reorthogonalization
    case 3
end


u1=K1\(f1+B1'*lambda);
u2=KS*(f2-B2'*lambda)+R2*alpha;  % KS statt "\" da KS schon berechnet!!
% löschen der doppelten interfaceDofs in u2 :

GG=abs(B1*u1)-abs(B2*u2);
if norm(GG) >10e-4
   disp('Gleichgewicht: B1*u1=B2*u2 nicht eingehalten')
end
[k,~]=find(B2'==1);
u2(k)=[];
u=[u1;u2];

end