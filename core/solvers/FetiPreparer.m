classdef FetiPreparer < SimpleSolvingStrategy
    %class that prepares problem for the PCPG Algorithm of Feti Level 1 and
    %2 Methods. E.g. find pseudo Inverse, Body Modes,... and later sends
    %the generated matrices to the Feti1/Feti2Solver to solve.
    
    properties
    end
    
    methods (Static)
        function [u] = solveFeti(K,f,substructures)
            
            %see whether the subdomains are singular or not
            singulars = [];
            for ii = 1:length(K)
                if length(K{1,ii}) ~= rank(K{1,ii})
                    method = 'singular';
                    singulars = [singulars ii];
                end
            end
            
            intfNodes = FetiPreparer.findInterfaceNode(substructures);
            
            if strcmp(method, 'singular')    
                %Rigid Body Modes
                R = cell(1,length(K));
                %pseudo inverse
                Kinv = cell(1,length(K));
                %Boolean Matrix
                B = cell(1,length(K));
                %test variable for correct pseudo Inverse/Body Modes
                test = 0;
                for it = 1:length(K)
                    if ismember(it, singulars)
                        %find pseudo Inverse/body modes for singulars
                        [r, kinv, test] = FetiPreparer.psInvRBM(K{1,it},test);
                        R{1,it} = r;
                        Kinv{1,it} = kinv;
                    else
                        %pseudo Inverse/body modes set to zero for
                        %non-singular subdomain
                        R{1,it} = 0;
                        Kinv{1,it} = 0;
                    end
                    %create Boolean Matrix
                    B{1,it} = FetiPreparer.createBooleanMatrix(substructures(it), intfNodes, it);
                end
                
                %f{1,1}
                
                %test for correct Pseudo Invers and Rigid Body Modes
                if test ~= 2*length(singulars)
                    disp('ERROR: Pseudo Invers and/or Rigid Body Modes are not correct!');
                end
                
                %all information given to the pCG is in cell array format
                %with one row and number of columns corresponding to the
                %number of subdomains
                %[lambda, alpha] = Feti1Solver.pCG(K,f,B,R,Kinv);
                [lambda, alpha] = Feti2Solver.pCG(K,f,B,R,Kinv);
                           
                %u01 = inv(K01)*(f01'+B01'*lambda);
                %u02 = K02Plus*(f02'-B02'*lambda)+R02*alpha;
                

             else
                lenU = 0;
                %define size for u
                for jj = 1:length(f)
                    lenU = lenU+length(f{jj});
                end
                
                u = zeros(lenU, 1);
                ll = 1;
                %assign displacements to u by doing matrix left division
                for kk = 1:length(K)
                    u(ll:ll+length(f{kk})-1,1) = K{kk}\f{kk}';
                    ll = ll+length(f{kk});
                end
            end
        end
        
        function intfNodes = findInterfaceNode(substructures)
            
            %find interface nodes of the substructures: intfNodes is a cell
            %array which contains all interface nodes of a substructure in
            %a row. The columns indicate with which other substructure the
            %interface exists.
            intfNodes = {};
            %iterate over all substructures
            for itS = 1:length(substructures)
                nodes = substructures(itS).getAllNodes;
                %iterate over all substructures not chosen now 
                for itOn = 1:length(substructures)
                    copies = [];
                    if itOn ~= itS
                        subsNodes = substructures(itOn).getAllNodes;
                    else
                        continue;
                    end
                    %iterate over all nodes of substructure of interest
                    for itN = 1:length(nodes)
                        %find interface nodes
                        nodeIntf = findIntfNode(nodes(itN),subsNodes);
                        if nodeIntf.getId ~= 0
                            copies = [copies nodeIntf];
                        end
                    end
                    intfNodes{itS, itOn} = [copies];
                end
            end
        end
        
        function [R, pseudoInv, test] = psInvRBM(K, test)
            %function coordinating the cholesky decomposition, pseudo inverse
            %generation and also the generation of the rigid body modes
            
            %do cholesky decomposition to get cholesky factors and Kpr
            %needed to compute body modes
            
            [KppFactors, Kpr, clm] = Feti1Solver.choleskyDecomp(K);

            %matrix of rigid body modes is created
            R = Feti1Solver.createBodyModesMatrix(KppFactors, Kpr, clm);

            %create pseudo inverse
            pseudoInv = Feti1Solver.createPseudoInv(KppFactors, clm);

            %test for rigid body modes. 10 decimal places
            nullspace = round(K*R, 10);
            if nullspace == 0
                test = test+1;
            end

            %test for pseudo inverse. 5 decimal places
            a = round(K,5);
            b = round(K*pseudoInv*K,5);
            if ismember(b,a) == 1
                test = test+1;
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
            len = length(clm)-1;
            R = KppFactors\Kpr;
            R = [R; zeros(len,len)];
            
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
        
        %Function creating the Boolean Matrix, named 'B' here, for each 
        %substructure. As we are dealing with subdomains, which have 
        %interfaces with multiple other subdomains we first have to 
        %evaluate the Boolean Matrices of each interface configuration 
        %between the different subdomainsm these Boolean Matrices are 
        %called 'b' here. Once these have been created one can calculate 
        %the total Boolean Matrix for the substructure by adding up all the 
        %individual Boolean Matrices.

        function B = createBooleanMatrix(substructure, intfNodes, it)

            %all nodes of substructure
            allNodes = substructure.getAllNodes;
            
            %find all interface nodes for all interface configurations
            nodesIntf = [];
            for numSubs = 1:length(intfNodes)
                nodesIntf = unique([nodesIntf intfNodes{it,numSubs}]);
            end
            %sort interface nodes by Id
            [~,idx] = sort([nodesIntf.getId]);
            nodesIntf = nodesIntf(idx);
            
            %find all dofs of substructure which are free and put them in 
            %array
            dofsFree = [];
            ii = 0;
            for numNode = 1:length(allNodes)
                dofs = allNodes(numNode).getDofArray;
                %only valid for 2D
                for dof = 1:2
                    if ~dofs(dof).isFixed
                        ii = ii+1;
                        dofsFree(1,ii) = ii;
                    end
                end
            end
            %total number of free dofs
            numDofs = numel(dofsFree);
            
            %find all interface dofs which are free
            dofsIntf = [];
            ii = 0;
            for numNode = 1:length(nodesIntf)
                dofs = nodesIntf(numNode).getDofArray;
                %only valid for 2D
                for dof = 1:2
                    if ~dofs(dof).isFixed
                        ii = ii+1;
                        dofsIntf(1,ii) = ii;
                    end
                end
            end
            %total number of interface dofs
            numIntfDofs = numel(dofsIntf);
 
            %find interface dofs for each interface configuration by
            %comparing them with all free dofs. This results in an array
            %giving the horizontal positions for the Boolean Matrix entries
            %corresponding to the interface dofs. 
            intfIndices2 = cell(1,length(intfNodes));
            for numSubs = 1:length(intfNodes)
                jj = 0;
                indices2 = [];
                if numSubs ~= it
                    for numNode = 1:length(allNodes)
                        dofs = allNodes(numNode).getDofArray;
                        %only valid for 2D
                        for dof = 1:2
                            if ~dofs(dof).isFixed
                                jj = jj+1;
                                %check whether the current node is an
                                %interface node for this configuration
                                if ismember(allNodes(numNode),intfNodes{it,numSubs})
                                    indices2 = [indices2 jj];
                                end
                            end
                        end
                    end
                end
                intfIndices2{1,numSubs} = indices2;
            end
            
            %find interface dofs for each interface configuration by
            %comparing them with all free interface dofs. This results in
            %an array giving the vertical positions for the Boolean Matrix
            %entries corresponding to the interface dofs.
            intfIndices1 = cell(1,length(intfNodes));
            for numSubs = 1:length(intfNodes)
                jj = 0;
                indices1 = [];
                if numSubs ~= it
                    for numNode = 1:length(nodesIntf)
                        dofs = nodesIntf(numNode).getDofArray;
                        %only valid for 2D
                        for dof = 1:2
                            if ~dofs(dof).isFixed
                                jj = jj+1;
                                %check whether the current node is an
                                %interface node for this configuration
                                if ismember(nodesIntf(numNode),intfNodes{it,numSubs})
                                    indices1 = [indices1 jj];
                                end
                            end
                        end
                    end
                end
                intfIndices1{1,numSubs} = indices1;
            end
             
            %number of interface nodes in the total system before this 
            %substructure. meaning all the interface nodes of the
            %subdomains preceding this one.
            dofsBefore = 0;
            for numNodes = 1:it-1
                nodes = [];
                for numrows = 1:length(intfNodes)
                    nodes = unique([nodes intfNodes{numNodes,numrows}]);
                end
                for ii = 1:length(nodes)
                    dofs = nodes(ii).getDofArray;
                    for jj = 1:length(dofs)
                        if ~dofs(jj).isFixed
                            dofsBefore = dofsBefore + 1;
                        end
                    end
                end
            end

            %number of interface nodes in the total system after this 
            %substructure. meaning all the interface nodes of the
            %subdomains following this one.
            dofsAfter = 0;
            for numNodes = it+1:length(intfNodes)
                nodes = [];
                for numrows = 1:length(intfNodes)
                    nodes = unique([nodes  intfNodes{numNodes,numrows}]);
                end
                for ii = 1:length(nodes)
                    dofs = nodes(ii).getDofArray;
                    for jj = 1:length(dofs)
                        if ~dofs(jj).isFixed
                            dofsAfter = dofsAfter + 1;
                        end
                    end
                end
            end

            %knowing the positions of all dofs and the interface dofs one
            %can assemble the c matrix for each interface. to do this one
            %sets all the entries of c to 1, which correspond to an
            %interface dof.
            c = cell(1,length(intfNodes));
            for numSubs = 1:length(intfNodes)
                d = [];
                if numSubs ~= it
                    %only valid for dimension 2, this is to exclude the
                    %crosspoints, they are evaluated later on
                    if numel(intfIndices1{1,numSubs}) ~= 2
                        d = zeros(numIntfDofs, numDofs);
                        intfDofs1 = intfIndices1{1,numSubs};
                        intfDofs2 = intfIndices2{1,numSubs};
                            for numIntfDof = 1:length(intfDofs2)
                                d(intfDofs1(numIntfDof), intfDofs2(numIntfDof)) = 1;
                            end
                    else
                        %include case where there is a crosspoint!!!
                    end
                end
                c{1,numSubs} = d;
            end

            %now a matrix b, which is the Boolean Matrix for each interface 
            %configuration of the whole subdomain, is created and saved for
            %each configuration in a cell array.
            b = cell(1,length(intfNodes));
            for numSubs = 1:length(intfNodes)
                if ~isempty(c{1,numSubs})
                    b{1,numSubs} = [zeros(dofsBefore,numDofs); c{1,numSubs}; zeros(dofsAfter,numDofs)];     
                else
                    b{1,numSubs} = [];
                end
            end
            
            %create B, by adding the various Boolean Matrices for the
            %different configurations of interfaces for the substructure
            B = 0;
            for numb = 1:length(b)
                if ~isempty(b{1,numb})
                    %B = b{1,numb};
                    B = B + b{1,numb};
                end
            end
        end 
    end 
end