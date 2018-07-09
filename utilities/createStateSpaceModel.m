function [ sys ] = createStateSpaceModel( femModel, observationDof )
%CREATESTATESPACEMODEL Creates a state space representation of the femModel
%   The output is obtained from at observationDof
[~,K] = SimpleAssembler.assembleGlobalStiffnessMatrix(femModel);
D = SimpleAssembler.assembleGlobalDampingMatrix(femModel);
M = SimpleAssembler.assembleGlobalMassMatrix(femModel);
[~,f] = SimpleAssembler.applyExternalForces(femModel);

[freeDofs, fixedDofs] = femModel.getDofConstraints;
if ~ isempty(fixedDofs)
    fixedDofIds = fixedDofs.getId();
    D = applyMatrixBoundaryConditions(D, fixedDofIds);
    M = applyMatrixBoundaryConditions(M, fixedDofIds);
end
N = length(freeDofs);

C = zeros(1,N+length(fixedDofs));
C(observationDof.getId) = 1;
C = applyVectorBoundaryConditions(C, fixedDofIds);

E = [eye(N) zeros(N,N); zeros(N,N) M]; 
A = [zeros(N,N) eye(N); -K -D] ; 
B = [zeros(N,size(f,1)); f.']; 
C = [C, zeros(size(C,1), N)];

sys = dss(A,B,C,0,E);
end

