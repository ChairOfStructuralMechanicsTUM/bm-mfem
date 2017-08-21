classdef Feti1Solver < FetiPreparer
    %   This class can solve system with the FETI level 1 method
    
    properties
    end
    
    methods (Static)
        
        function [lambda, alpha] = pCG(K,f,B,R,Kinv)
            
            B{1,1} = [0 0 0 0 0 0 0 1 0 0 0
                      0 0 0 0 0 0 0 0 1 0 0
                      0 0 0 0 0 0 0 0 0 1 0
                      0 0 0 0 0 0 0 0 0 0 1];
                  
            B{1,2} = [0 0 0 0 0 0 0 0 0 0 -1 0 0 0
                      0 0 0 0 0 0 0 0 0 0 0 -1 0 0
                      0 0 0 0 0 0 0 0 0 0 0 0 -1 0
                      0 0 0 0 0 0 0 0 0 0 0 0 0 -1];
       
            
            
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
            
            %%%Initializing
            %lambda
            lambda(:,1) = G*inv(G'*G)*e;
            %residuum-w
            w(:,1) = P'*(d-F*lambda(:,1));
            
            %%%Tests
            test1 = G'*lambda(:,1);
            for ii = 1:ns
                test2 = R{1,ns}'*f{1,ns}';
            end
            
            
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
                    if abs(w(numEntries,kk)) > 10^-3
                        var = 0;
                        kk = kk+1;
                        break;
                    else
                        var = 1;
                    end
                end        
            end
            lambda = lambda(:,kk);
            w = w(:,kk);
            kk;
            
            %Tests
            test3 = R{1,2}'*(f{1,2}'-B{1,2}'*lambda);
            
            %is alpha general quantity or does one have to calculate it for
            %every subdomain individually?
%             for numSubs = 1:ns
%                 if Kinv{1,numSubs} == 0
%                     
%                 else
%                     g = B{1,numSubs}*R{1,numSubs};
%                     fi = B{1,numSubs}*Kinv{1,numSubs}*B{1,numSubs}';
%                     di = B{1,numSubs}*Kinv{1,numSubs}*f{1,numSubs}';
%                     alpha{1,numSubs} = (inv(g'*g)*g')*(fi*lambda-di);
%                 end
%             end
                g = B{1,2}*R{1,2};
                %fi = B{1,1}*inv(K{1,1})*B{1,1}'+B{1,2}*Kinv{1,2}*B{1,2}';
                di = B{1,1}*inv(K{1,1})*f{1,1}'+B{1,2}*Kinv{1,2}*f{1,2}';
                alpha{1,2} = (inv(g'*g)*g')*(F*lambda-di);
        end
    end
end
