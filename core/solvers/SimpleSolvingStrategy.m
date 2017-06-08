classdef SimpleSolvingStrategy < Solver
    %SIMPLESOLVINGSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
%         nodalDisplacments
    end
    


    %     methods
        
        %     constructor
%         function nodalDisplacements = SimpleSolvingStrategy(femModel)
%             if (nargin > 0)
%                 nodalDisplacements.res = SimpleSolvingStrategy.solve(femModel);
%             end
%         end
%     end
    
    
    
methods (Static)
    function x = solve(femModel)
        Kred = SimpleAssembler(femModel).reducedStiffnessMatrix;
        f = SimpleAssembler(femModel).reducedForceVector;
        
        
        %Chol = chol(Kred)
        [R, pseudoInv] = SimpleSolvingStrategy.psInvRBM(Kred);

        Kred2 = Kred*pseudoInv*Kred;
        
        %set quasi zero elements < 10^-12 to zero
        for ii=1:length(Kred2)
            for jj = 1:length(Kred2)
                if abs(Kred2(ii,jj)) < 10^-12
                    Kred2(ii,jj) = 0;
                end
            end
        end
        
        %comparison of stiffness matrices
        Kred2;
        Kred;
        
        %check whether body modes form the null space
        Kred*R;

        x = Kred \ f';
        
        SimpleAssembler.assignResultsToDofs(femModel, x);
    end
    
    function [ R, pseudoInv ] = psInvRBM( mat )
        %UNTITLED6 Summary of this function goes here
        %function coordinating the cholesky decomposition, pseudo inverse
        %generation and also the generation of the rigid body modes

        %do cholesky decomposition
        cholesky = SimpleSolvingStrategy.choleskyDecomp(mat);

        %eliminate zero-pivot rows and columns and find body modes
        [kppFactors,kpr, cnt] = SimpleSolvingStrategy.findBodyModes(cholesky);

        %set full rank submatrix of Kred, named kpp
        kpp = kppFactors'*kppFactors;

        %find the pseudoInverse from kpp
        pseudoInv = SimpleSolvingStrategy.createPseudoInv(kpp, cnt);

        %generate matrix of rigid body modes  from the Factors of kpp and the kpr
        %part of the stiffness matrix found during the elimination of the pivots
        R = kppFactors\kpr;
        R = [R; eye(cnt)];
    end
    
    function [ kppFactors, kpr, cnt ] = findBodyModes( mat )
        %UNTITLED2 Summary of this function goes here
        % Functions that eliminates the redundant rows and columns of the stiffness
        % matrix and finds the rigid body modes
        n = length(mat);
        ll = 1;
        ii = 1;
        jj = 1;

        %factos of Kpp
        kppFactors = mat;
        cnt = 0;

        while ii <= n
            %find the zero-pivot elements
            if  kppFactors(ii,ii) == 1
                kpr(1:ii-1,ll) = -mat(1:ii-1,jj);
                %kpr(jj, ll) = 1;
                kppFactors = removerows(kppFactors, 'ind', ii);
                kppFactors = (removerows(kppFactors', 'ind', ii))';
                cnt = cnt+1;
                ll = ll+1;
                n = n-1;
            else
                ii = ii+1;
            end
            jj = jj+1;
        end
    end
    
    function [ pseudoInv ] = createPseudoInv( kpp, cnt )
        %UNTITLED2 Summary of this function goes here
        %   from kpp the pseudoInvers is created by taking kpp^-1 and adding zero
        %   matrices for as many times as there were zero pivots in the stiffness
        %   matix

        pseudoInv = inv(kpp);
        for ii = 1:cnt
            pseudoInv = [pseudoInv, zeros(length(pseudoInv),1); zeros(1, length(pseudoInv)+1)];
        end        
    end
    
    function [ cholesky ] = choleskyDecomp( mat )
        %UNTITLED Summary of this function goes here
        %   Function performing the cholesky decomposition for positive definite
        %   matrices. For positive semidefinite ones an approcimation of the 
        %   cholesky decomposition is done.

        n = length(mat);
        cholesky = zeros(n,n);

        for ii = 1:n
           %create diagonal elemens of cholesky matrix 
           cholesky(ii,ii) = sqrt(mat(ii,ii) - cholesky(:,ii)'*cholesky(:,ii));

           %in case of a zero pivot element set it to 1
           if cholesky(ii,ii) < 10^-6
               cholesky(ii,ii) = 1;   
           else
               %for non singular parts do normal cholesky decomposition
               for jj=(ii+1):n
                   cholesky(ii,jj) = (mat(ii,jj) - cholesky(:,ii)'*cholesky(:,jj))/cholesky(ii,ii);
               end
           end
        end
    end
 
    function nodalForces = getNodalForces(femModel)
        
        dofs = femModel.getDofArray;
        nDofs = length(dofs);
        nodalForces = SimpleAssembler(femModel).forceVector;
        fixedDofs = [];
        
        for itDof = 1:nDofs
            if (dofs(itDof).isFixed)
                fixedDofs = [fixedDofs itDof];                          % array of fixed dofs and their location
            else
                nodalForces(itDof) = dofs(itDof).getValue;
            end
        end
        
        nodalForces(fixedDofs) = SimpleAssembler(femModel).stiffnessMatrix(fixedDofs, :) * getValue(femModel.dofArray)';
        
    end
    
end
end

