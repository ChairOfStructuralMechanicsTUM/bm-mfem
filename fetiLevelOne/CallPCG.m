% Solve Substructure:

function [u]=CallPCG(Substructure01,Substructure02)
StiffnessMatrixSubstrucuture01=SimpleAssembler(Substructure01);
StiffnessMatrixSubstructure02=SimpleAssembler(Substructure02);

K1=StiffnessMatrixSubstrucuture01.reducedStiffnessMatrix;
K2=StiffnessMatrixSubstructure02.reducedStiffnessMatrix;
f1=StiffnessMatrixSubstrucuture01.reducedForceVector';
f2=StiffnessMatrixSubstructure02.reducedForceVector';

[R2,KS]=computeRidgedBodyModes(K2);
[B1]=computeB(Substructure01);
[B2]=computeB(Substructure02);

% interface system:
F1=B1*K1^-1*B1'+B2*KS*B2';
Gi2=B2*R2;
d=B2*KS*f2-B1*K1^-1*f1;
e=R2'*f2;

%preconditioner
Pi=[B1 B2]*[K1, zeros(size(K1,1),size(K2,2));...
    zeros(size(K2,2),size(K1,1)) K2]*[B1';B2'];

a=2;
switch a
    case 1
        [lambda,alpha]=CG(F1,Gi2,d,e); % CG
    case 2
        [lambda,alpha]=PCG(F1,Pi,Gi2,d,e); % Precond. CG
end

u1=K1\(f1+B1'*lambda);
u2=KS*(f2-B2'*lambda)+R2*alpha;

[U1]=assignUToDof(Substructure01,u1);
[U2]=assignUToDof(Substructure02,u2);

% delte interfaceDofs on u2
[k,~]=find(abs(B2')==1);
U2(k,:)=[];
U=[U1;U2];
% sorting DofValues by Node.Ids
[~,I]=sort(U(:,1));
u=U(I,2);

end