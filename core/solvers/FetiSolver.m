classdef FetiSolver < SimpleSolvingStrategy
    %   This class can solve system with the Feti method
    
    properties
    end
    
    methods (Static)
        
        function [u01, u02] = solveFeti(K01, K02, f01, f02, femModel01, femModel02)
            
            rank01 = rank(K01);
            rank02 = rank(K02);
            
            if rank01 == length(K01) && rank02 == length(K02)
                %non-singular
                disp('non-singular');
                
                u01 = K01\f01';
                u02 = K02\f02';
                
            elseif rank01 == length(K01) && rank02 ~= length(K02)
                %K_02 singular
                disp('K02 singular');
                

                [R02, K02Plus] = FetiSolver.psInvRBM(K02);
                B01 = FetiSolver.createBooleanMatrix(femModel01);
                B02 = FetiSolver.createBooleanMatrix(femModel02);
                
                [lambda, alpha] = FetiSolver.projectedConjugateGradient(K01, f01, B01, K02, K02Plus, f02, B02, R02);
                
                u01 = K01\(f01'+B01'*lambda);
                u02 = K02Plus*(f02'-B02'*lambda)+R02*alpha;
                
                B01*u01;
                B02*u02;
                
                
            elseif rank01 ~= length(K01) && rank02 == length(K02)
                %K_01 singular
                disp('K01 singular');

                [R01, K01Plus] = FetiSolver.psInvRBM(K01);
                B01 = FetiSolver.createBooleanMatrix(femModel01);
                B02 = FetiSolver.createBooleanMatrix(femModel02);
                
                [lambda, alpha] = FetiSolver.projectedConjugateGradient(K02, f02, B02, K01, K01Plus, f01, B01, R01);
                
            else
                %all singular
                disp('all singular');

                [R01, K01Plus] = FetiSolver.psInvRBM(K01);
                [R02, K02Plus] = FetiSolver.psInvRBM(K02);
                B01 = FetiSolver.createBooleanMatrix(femModel01);
                B02 = FetiSolver.createBooleanMatrix(femModel02);
                
                
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

            %test for rigid body modes. 10 decimal places
            nullspace = round(K*R, 10);
            if nullspace == 0
                disp('Rigid body modes are fine.');
            else
                disp('Rigid body modes are not correct.');
            end
            
            %test for pseudo inverse. 5 decimal places
            a = round(K,5);
            b = round(K*pseudoInv*K,5);
            if ismember(b,a) == 1
                disp('Pseudo inverse is fine.');
            else 
                disp('Pseudo inverse is not correct.');
            end                    
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
            
            %PROBLEM: maybe the assembly of the R matrix is not completly
            %right. needs to be further tested

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
            
            %PROBLEM: maybe the assembly of the matrix pseudoInv is not
            %excatly right. needs to be tested further.
        
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

        
        function [lambda, alpha] = projectedConjugateGradient(K01, f01, B01, K02, Kplus, f02, B02, R02)
            %function to perform the projected conjugate gradient method in
            %order to solve the singular system of equations. non-singular
            %system is always called 01, singular system 02.
            
            
            %G02 is the restriction of the body modes on the interface
            G02 = B02*R02;
            F01 = B01*inv(K01)*B01' + B02*Kplus*B02';
            %change this preconditioner!!!
            F01Inv = B01*K01*B01' + B02*K02*B02';
            
            %check definitness of F01: is positive definite
            [~, q] = chol(F01);
            rank(F01);
            size(F01,1);
            
            %projection operator 
            Proj = [eye(size(G02,1))-G02*inv(G02'*G02)*G02'];

            
            %From "A Method of FETI..."
            %arbitrary vector
            %lambda(:,1) = G02*inv(G02'*G02)*R02'*f02';
            
            %option 1
            %r(:,1) = B02*Kplus*f02'-B01*inv(K01)*f01';
            
            
            
            %option 2
            %r(:,1) = F01*lambda(:,1)+B01*inv(K01)*f01'-B02*Kplus*f02';
            
            
            %From "Theory and Implementation..."
%             lambda(:,1) = G02*inv(G02'*G02)*(f02*R02)'; %correct
%             r(:,1) = B01*inv(K01)*f01'+B02*Kplus*f02'-F01*lambda(:,1); %correct            
%             w(:,1) = Proj*r(:,1);
%             %y(:,1) = P*inv(F01Pre)*w(:,1);
%             p(:,1) = w(:,1); %correct



                %%%BEST ONE
            %From The second generation
            lambda(:,1) = G02*inv(G02'*G02)*(f02*R02)';
            %w is the residuum
            w(:,1) = Proj'*((B01*inv(K01)*f01'+B02*Kplus*f02')-F01*lambda(:,1));
            
            %constraint check
            G02'*lambda(:,1);
            r02f02= R02'*f02';
            
                
            kk = 1;
            %while r(kk+1) > 10^-3
            
            while kk < 4
                

                %%%BEST ONE
                y(:,kk) = Proj*F01Inv*w(:,kk);
                
                
                if kk > 1
                    temp = 0;
                    for mm = 1:kk-1
                        temp = temp + ((y(:,kk)'*F01*p(:,mm))/(p(:,mm)'*F01*p(:,mm)))*p(:,mm);
                    end
                    p(:,kk) = y(:,kk)-temp;   
                else
                    p(:,kk) = y(:,kk);
                end
                
                
                n(:,kk) = (p(:,kk)'*w(:,kk))/(p(:,kk)'*F01*p(:,kk));
                lambda(:,kk+1) = lambda(:,kk) + n(:,kk)*p(:,kk);
                w(:,kk+1) = w(:,kk) - n(:,kk)*Proj'*F01*p(:,kk);

                
                
%                 %constraint check
                g02 = G02'*lambda(:,kk);
                test = G02'*p(:,kk);
                
                kk = kk+1;
                
            end
            

            
            lambda;
            w;
            
            lambda = lambda(:,kk);
            
            
            %lambda is improvised here
            alpha = inv(G02'*G02)*G02'*(F01*lambda-B02*Kplus*f02'+B01*inv(K01)*f01');  
        
        end
        
    end
    
end




                  %%%From "A Method of FETI..."
%                 %disp('check');
%                 
%                 %From "A Method of FETI..."
%                 if kk == 2
%                     beta(kk) = 0;
%                     s(:,kk) = r(:,kk-1);   
%                 else
%                     beta(kk) = (r(:,kk-1)'*r(:,kk-1))/(r(:,kk-2)'*r(:,kk-2));
%                     s(:,kk) = r(:,kk-1)+beta(:,kk)*s(:,kk-1);
%                     %s(:,kk) = P*s(:,kk);
%                 end
%                 
%                 s(:,kk) = P*s(:,kk);
%                 
%                 gamma(kk) = (r(:,kk-1)'*r(:,kk-1))/(s(:,kk)'*F01*s(:,kk));
%                 
%                 lambda(:,kk) = lambda(:,kk-1)+gamma(kk)*s(:,kk);
%                 r(:,kk) = r(:,kk-1)-gamma(:,kk)*F01*s(:,kk);
%                
%                 %constraint check
%                 kk;
%                 g02 = G02'*lambda(:,kk); 
%                 G02'*s(:,kk);
%                 
%                 
%                 
%                 kk = kk+1;





%                 %%%From "Theory and Implementation..."
% 
%                 
%                 gamma(kk-1) = dot(w(:,kk-1),w(:,kk-1))/dot(p(:,kk-1), F01*p(:,kk-1));
%                 
%                 lambda(:,kk) = lambda(:,kk-1)+gamma(kk-1)*p(:,kk-1);
%                 
%                 r(:,kk) = r(:,kk-1)-gamma(kk-1)*F01*p(:,kk-1);
%                 
%                 w(:,kk) = Proj*r(:,kk);
%                 
%                 %y(:,kk) = P*inv(F01Pre)*w(:,kk);
%                 
%                 beta(kk) = dot(w(:,kk), w(:,kk))/dot(w(:,kk-1), w(:,kk-1));               
%                 p(:,kk) = w(:,kk)+beta(kk)*p(:,kk-1);
%                 
%                 %term(kk) = (y(:,kk)'*F01*p(:,kk-1))/(p(:,kk-1)'*F01*p(:,kk-1));
%                 %sum(term,2);
%                 %p(:,kk) = y(:,kk)- sum(term,2)*p(:,kk-1);
%                 
%                 %constraint check
% %                 kk
%                 g02 = G02'*lambda(:,kk);
%                 G02'*p(:,kk);
%                 
%                 kk = kk+1;


