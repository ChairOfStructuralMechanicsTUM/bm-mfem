classdef FetiSolver < SimpleSolvingStrategy
    %   This class can solve system with the Feti method
    
    properties
    end
    
    methods (Static)
        
        function solve = solveFeti(K_01, K_02, f_01, f_02)
            
            rank_01 = rank(K_01);
            rank_02 = rank(K_02);
            
            if rank_01 == length(K_01) && rank_02 == length(K_02)
                %non-singular
                disp('non-singular');

            elseif rank_01 == length(K_01) && rank_02 ~= length(K_02)
                %K_02 singular
                disp('K_02 singular');

                [R_02, K02_plus] = FetiSolver.psInvRBM(K_02);
                
            elseif rank_01 ~= length(K_01) && rank_02 == length(K_02)
                %K_01 singular
                disp('K_01 singular');

                [R_01, K01_plus] = FetiSolver.psInvRBM(K_01);
                
            else
                %all singular
                disp('all singular');

                [R_01, K01_plus] = FetiSolver.psInvRBM(K_01);
                [R_02, K02_plus] = FetiSolver.psInvRBM(K_02);
                
            end
        end
        
        function [R, pseudoInv] = psInvRBM(K)
            %function coordinating the cholesky decomposition, pseudo inverse
            %generation and also the generation of the rigid body modes
            
            %do cholesky decomposition to get cholesky factors and Kpr
            %needed to compute body modes
            [KppFactors, Kpr, clm] = FetiSolver.choleskyDecomp(K);

            %matrix of rigid body modes is created
            R = FetiSolver.createBodyModesMatrix(KppFactors, Kpr, clm);

            %create pseudo inverse
            pseudoInv = FetiSolver.createPseudoInv(KppFactors, clm);

            %tests
%             originalK = K
%             testK = K*pseudoInv*K
%             testR = K*R
        end
        
        function [KppFactors, Kpr, clm] = choleskyDecomp(K)
            %   Function performing the cholesky decomposition. The matrix
            %   Kpr needed to compute the rigid body modes is calculated as
            %   well.
            
            n = length(K);
            KppFactors = zeros(n,n);
            %indice for KppFactors
            ii = 1;
            %indice for stiffness matrix
            mm = 1;
            %counts the times a row and column is deleted
            cnt = 0;
            Kpr = [];
            %safes the indices of rows and columns that are deleted
            clm = [1];

            while ii <= n
               %diagonal elements of cholesky decomposition 
               KppFactors(ii,ii) = sqrt(K(mm,mm) - KppFactors(1:(ii-1),ii)'*KppFactors(1:(ii-1),ii));

               %for zero pivot delete row, column and safe their indice
                if KppFactors(ii,ii) < 10^-6

                    cnt = cnt+1;
                    Kpr(1:(ii-1),cnt) = -KppFactors(1:(ii-1),ii);
                    %delete rows and columns
                    KppFactors = removerows(KppFactors, 'ind', ii);
                    KppFactors = (removerows(KppFactors', 'ind', ii))';
                   
                    %remember which rows, columns are deleted
                    clm = [clm, mm];

                    %decrease ii and n after deleting a row
                    n = n-1;
                    ii = ii-1;
                else
                   
                    %for non singular parts do normal cholesky decomposition
                    jj = ii+1;
                    nn = mm+1;
                    while jj <= n
                        KppFactors(ii,jj) = (K(mm,nn) - KppFactors(1:(ii-1),ii)'*KppFactors(1:(ii-1),jj))/KppFactors(ii,ii);
                        jj = jj+1;
                        nn = nn+1;
                    end
                end
                ii = ii+1;
                mm = mm+1;
            end
        end

        function R = createBodyModesMatrix(KppFactors, Kpr, clm)
            %function creating the matrix of rigid body modes

            %backward substitution on R
            R = KppFactors\Kpr;
            R = [R; zeros(3,3)];
            
            %insert identity matrix under/into body modes. They need to be
            %at the position of the deleted zero-pivot
            for ii = 2:length(clm)
                if clm(ii) < size(R,1)
                    %safe row in which one for zero-pivot is inserted
                    temp = R(clm(ii),:);
                    %set row of zero-pivot to zero and move entries one row
                    %down. this implies only identity one is in that row.
                    R(clm(ii),:) = zeros(1, size(R,2));
                    R(clm(ii)+1,:) = temp;
                    %set zero-pivot location to 1
                    R(clm(ii),ii-1) = 1;
                %in last row just add the one instead of zero-pivot
                else
                     R(clm(ii),ii-1) = 1;
                end
            end
        end
                      
        function pseudoInv = createPseudoInv(KppFactors, clm)
            %function that creates the pseudo inverse
        
            %get full rank submatrix of stiffness matrix
            Kpp = KppFactors'*KppFactors;
            KppInv = inv(Kpp);
            pseudoInv = [];

            for ii = 2:length(clm)      
                %if deleted rows/columns have non redundant
                %rows/columns in between
                if clm(ii)-clm(ii-1) > 1            
                    diff = clm(ii)-clm(ii-1);
                    %first iteration create basic pseudoInv
                    if ii == 2
                        basePart = [KppInv(1:(clm(ii)-1), clm(ii-1):(clm(ii)-1))];
                        pseudoInv = [pseudoInv, basePart];
                    %all further iterations
                    else
                        if clm(ii) <= size(KppInv,2)
                            cols = clm(ii);
                        else
                            cols = size(KppInv,2);
                        end
                        %basePart is column of a non redundant column
                        %between redundant ones. basePart2 same for rows
                        basePart = [KppInv(1:(clm(ii-1)-1), clm(ii-1):cols); zeros(diff-1, diff-1)];
                        basePart2 = [KppInv(clm(ii-1), 1:clm(ii-1)-1), zeros(diff-1, diff-1), KppInv(clm(ii-1),clm(ii-1))];
                        pseudoInv = [pseudoInv, basePart; basePart2];
                    end
                    
                    %dimensions of pseudoInv 
                    high = size(pseudoInv,1);
                    len = size(pseudoInv,2); 
                    %insert zeros after a block of non-redundant rows,
                    %columns to display next redundant equation.
                    pseudoInv = [pseudoInv, zeros(high,1); zeros(1,len+1)];

                %another zero line/column inserted if a deleted row/column
                %directly follows another
                else
                    high = size(pseudoInv,1);
                    len = size(pseudoInv,2);
                    pseudoInv = [pseudoInv, zeros(high,1); zeros(1,len+1)]; 
                end

            end             
        end
        
        function booleanMatrix = createBooleanMatrix(femModel)
            %function creating the boolean matrix needed for the interface
            %problem 
            
            %all degrees of freedom in substructure
            totDof = femModel.getDofArray;
            int = length(totDof);
            
            %get number of interface nodal unknowns
            dofs = [];
            nodes = femModel.getNodesIntf;
            for jj = 1:length(femModel.getNodesIntf)
                dofs = [dofs, nodes(jj).getDofArray];
            end
            %number of total interface nodal unknowns
            intf = length(dofs);
            
            %count how many dof are fixed
            for ii = 1:length(totDof)
                if totDof(ii).isFixed
                    int = int-1;
                end
            end         
            for jj = 1:length(dofs)
                if dofs(jj).isFixed
                    intf = intf-1;
                end
            end
            %subtract interface unknowns from total unknowns
            int = int - intf;
            booleanMatrix = [zeros(intf, int), eye(intf, intf)];
        end

        
        function [lambda, alpha] = projectedConjugateGradient(K_01, f_01, B_01, K_plus, f_02, B_02)
            %function to perform the projected conjugate gradient method in
            %order to solve the singular system of equations
            
            lambda(1) = R_02*(R_02'*R_02)^-1*R_02'*f_02;
            r(1) = [B_02*K_plus*f_02-B_01*inv(K_01)*f_01];
            F_01 = B_01*inv(K_01)*B_01' + B_02*K_plus*B_02';
            kk = 1;
            
            while r(kk+1) > 10^-3
                if kk == 1
                    beta(1) = 0;
                    s(1) = r;
                else
                    beta(kk) = (r(kk)'*r(kk))/(r(kk-1)'*r(kk));
                    s(kk) = r(kk)+beta(kk)*s(kk-1);
                end
                
                gamma(kk) = (r(kk)'*r(kk))/(s(kk)'*F_01*s(kk));
                lambda(kk+1) = lambda(kk)+gamma(kk);
                r(kk+1) = r(kk)-gamma(kk)*F_01*s(kk);
               
                kk = kk+1;
            end
            
            alpha = inv((B_02*R_02)'*(B_02*R_02))*(F_01*lambda-B_02*K_plus*f_02+B_01*inv(K_01)*f_01);
        
        
        end
        
    end
    
end

