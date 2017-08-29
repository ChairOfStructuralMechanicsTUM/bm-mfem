classdef Feti1Solver < FetiPreparer
    %   This class can solve system with the FETI level 1 method
    
    properties
    end
    
    methods (Static)
        
        function [lambda, alpha] = pCG(K,f,B,R,Kinv)
            
            %%%Set matrices
            %ns = number of subdomains
            ns = length(K);
            %nf = number of floating subdomains
            nf = 0;
            for numSubs = 1:ns
                if Kinv{1,numSubs} == 0
                else
                    nf = nf+1;
                end
            end

            %F-Matrix: sum(B*K^-1*B') from 1 to ns
            F = 0;
            for numSubs = 1:ns
                if Kinv{1,numSubs} == 0
                    F = F + B{1,numSubs}*inv(K{1,numSubs})*B{1,numSubs}';
                else
                    F = F + B{1,numSubs}*Kinv{1,numSubs}*B{1,numSubs}';
                end
            end
            F;
            %definitness check
            rank(F);
            size(F,1);
          
            
            %d-Matrix: sum(B*K^-1*f') from 1 to ns
            d = 0;
            for numSubs = 1:ns
                if Kinv{1,numSubs} == 0
                    d = d + B{1,numSubs}*inv(K{1,numSubs})*f{1,numSubs}';
                else
                    d = d + B{1,numSubs}*Kinv{1,numSubs}*f{1,numSubs}';
                end
            end
            
            %G-Matrix: [B*R ... B*R] from 1 to nf
            G = [];
            for numSubs = 1:ns
                if Kinv{1,numSubs} == 0
                else
                    G = [G B{1,numSubs}*R{1,numSubs}];
                end
            end
            
            %e-Matrix: [f*R ... f*R]' from 1 to nf
            preE = [];
            for numSubs = 1:ns
                if Kinv{1,numSubs} == 0
                else
                    preE = [preE (f{1,numSubs}*R{1,numSubs})];
                end
            end
            e = preE';
            
            %P-Projektor: P = I-G*(G'*G)^-1*G'
            I = eye(size(G*G'));
            P = I-G*inv(G'*G)*G';
            
            %Lumped Preconditioner: try (sum(B*K*B'))
            pre = 0;
            for numSubs = 1:ns
                pre = pre + B{1,numSubs}*K{1,numSubs}*B{1,numSubs}';
            end
            lumpedF = pre;
            %Preconditioner Check
            %F*lumpedF;
            
            %%%Initializing
            %lambda
            lambda(:,1) = G*inv(G'*G)*e;
            %residuum-w
            w(:,1) = P'*(d-F*lambda(:,1));
            
            %%%Tests
%             test1 = G'*lambda(:,1);
%             for ii = 1:ns
%                 test2 = R{1,ns}'*f{1,ns}';
%             end
            
            
            %%%Iterate
            var = 0;
            kk = 1;
            while var < 1
            %for kk = 1:4
                y(:,kk) = P*lumpedF*w(:,kk);
                
                if kk > 1
                    temp = 0;
                    for mm = 1:kk-1
                        temp = temp + ((y(:,kk)'*F*p(:,mm))/(p(:,mm)'*F*p(:,mm)))*p(:,mm);
                    end
                    p(:,kk) = y(:,kk)-temp;   
                else
                    p(:,kk) = y(:,kk);
                end
                
                eta(:,kk) = (p(:,kk)'*w(:,kk))/(p(:,kk)'*F*p(:,kk));
                
                lambda(:,kk+1) = lambda(:,kk)+eta(:,kk)*p(:,kk);
                
                w(:,kk+1) = w(:,kk)-eta(:,kk)*P'*F*p(:,kk);
                
                for numEntries = 1:size(w,1)
                    if abs(w(numEntries,kk)) > 10^-13
                        var = 0;
                        kk = kk+1;
                        break;
                    else
                        var = 1;
                    end
                end 
            end
            lambda = lambda(:,kk);
            
%             %Tests
%             test3 = R{1,2}'*(f{1,2}'-B{1,2}'*lambda);

            %calculate alpha
            alphA = inv(G'*G)*G'*(F*lambda-d);
            row = 1;
            %this should be valid for all cases. n is used to assign the
            %correct values of alpha to the corresponding subdomain
            for numSubs = 1:ns
                n = size(R{1,numSubs},2);
                if Kinv{1,numSubs} == 0
                    alpha{1,numSubs} = 0;
                else
                    alpha{1,numSubs} = alphA(row:row+n-1,1);
                    row = row+n;
                end
            end
        end
    end
end
