classdef Feti2Solver < FetiPreparer
    %UNTITLED Summary of this class goes here
    %   This class can solve system with the FETI level 2 method
    
    properties
    end
    
    methods (Static)
        function [lambda, alpha] = pCG(K,f,B,R,Kinv,C,W)
            
            %change B matrix to try...
            B{1,1} = [0 1 zeros(1,13)
                      0 0 1 zeros(1,12)
                      zeros(1,7) 1 zeros(1,7)
                      zeros(1,8) 1 zeros(1,6)
                      zeros(1,9) 1 zeros(1,5)
                      zeros(1,10) 1 zeros(1,4)
                      zeros(1,11) 1 zeros(1,3)
                      zeros(1,12) 1 zeros(1,2)
                      zeros(1,13) 2 zeros(1,1)
                      zeros(1,14) 2
                      zeros(8,15)];

            B{1,2} = [zeros(1,12) -1 zeros(1,5)
                      zeros(1,13) -1 zeros(1,4)
                      zeros(1,14) -1 zeros(1,3)
                      zeros(1,15) -1 zeros(1,2)
                      zeros(4,18)
                      zeros(1,16) -2 0
                      zeros(1,17) -2
                      zeros(1,8) -1 zeros(1,9)
                      zeros(1,9) -1 zeros(1,8)
                      zeros(1,10) -1 zeros(1,7)
                      zeros(1,11) -1 zeros(1,6)
                      zeros(4,18)];
                  
            B{1,3} = [zeros(4,18)
                      zeros(1,12) -1 zeros(1,5)
                      zeros(1,13) -1 zeros(1,4)
                      zeros(1,14) -1 zeros(1,3)
                      zeros(1,15) -1 zeros(1,2)
                      zeros(1,16) -2 0
                      zeros(1,17) -2
                      zeros(4,18)
                      zeros(1,4) -1 zeros(1,13)
                      zeros(1,5) -1 zeros(1,12)
                      zeros(1,10) -1 zeros(1,7)
                      zeros(1,11) -1 zeros(1,6)];
                  
            B{1,4} = [zeros(8,18)
                      zeros(1,16) 2 0
                      zeros(1,17) 2
                      zeros(1,8) 1 zeros(1,9)
                      zeros(1,9) 1 zeros(1,8)
                      zeros(1,10) 1 zeros(1,7)
                      zeros(1,11) 1 zeros(1,6)
                      zeros(1,12) 1 zeros(1,5)
                      zeros(1,13) 1 zeros(1,4)
                      zeros(1,14) 1 zeros(1,3)
                      zeros(1,15) 1 zeros(1,2)];
            
%                   B{1,1}
%                   B{1,2}
%                   B{1,3}
%                   B{1,4}

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
            
            %C-Matrix
            C = [0;0;0;0;0;0;0;0;1;1;0;0;0;0;0;0;0;0];
            
            %Lumped Preconditioner: try W*(sum(B*K*B'))*W
            pre = 0;
            for numSubs = 1:ns
                pre = pre + B{1,numSubs}*K{1,numSubs}*B{1,numSubs}';
            end
            lumpedF = W*pre*W;
            %lumpedF = pre;
            
            %%%Initializing
            nomi = C'*P'*(d-F*G*inv(G'*G)*e);
            denomi = C'*P'*F*P*C;
            %gamma: only works for this case with only one crosspoint
            %otherwise a more complicated equation has to be solved for
            %gamme...
            gamma(1) = nomi/denomi;
            %lambda
            lambda(:,1) = G*inv(G'*G)*e+P*C*gamma(1);
            %residuum-w
            w(:,1) = P'*(d-F*lambda(:,1));
            
            %Preconditioner Check
            %F*lumpedF
            
            %%%Iterate
            var = 0;
            kk = 1;
            while var < 1
            %for kk = 1:100
                y(:,kk) = P*lumpedF*w(:,kk);
                
                if kk > 1
                    temp = 0;
                    for mm = 1:kk-1
                        temp = temp + ((y(:,kk)'*F*p2(:,mm))/(p2(:,mm)'*F*p2(:,mm)))*p2(:,mm);
                    end
                    p(:,kk) = y(:,kk)-temp;   
                else
                    p(:,kk) = y(:,kk);
                end
                
                p(:,kk);
                
                nomi = -C'*P'*F*p(:,kk);
                denomi = C'*P'*F*P*C;
                gamma(kk) = nomi/denomi;
                
                p2(:,kk) = p(:,kk)+P*C*gamma(kk);
                
                eta(:,kk) = (p2(:,kk)'*w(:,kk))/(p2(:,kk)'*F*p2(:,kk));
                
                lambda(:,kk+1) = lambda(:,kk)+eta(:,kk)*p2(:,kk);
                
                w(:,kk+1) = w(:,kk)-eta(:,kk)*P'*F*p2(:,kk);
                
                for numEntries = 1:size(w,1)
                    if abs(w(numEntries,kk)) > 10^-10
                        var = 0;
                        kk = kk+1;
                        break;
                    else
                        var = 1;
                    end
                end
            end
            lambda = lambda(:,kk);
           
            %test
%             R{1,2}'*(f{1,2}'-B{1,2}'*lambda)
%             R{1,3}'*(f{1,3}'-B{1,3}'*lambda)
%             R{1,4}'*(f{1,4}'-B{1,4}'*lambda)
            
            %is alpha general quantity or does one have to calculate it for
            %every subdomain individually?
%             for numSubs = 1:ns
%                 if Kinv{1,numSubs} == 0
%                 else
%                     g = B{1,numSubs}*R{1,numSubs};
%                     %fi = B{1,numSubs}*Kinv{1,numSubs}*B{1,numSubs}';
%                     di = B{1,numSubs}*Kinv{1,numSubs}*f{1,numSubs}';
%                     alpha{1,numSubs} = (inv(g'*g)*g')*(F*lambda-d);
%                 end
%             end

%             alpha{1,1} = 0;
% 
%             di = 2*B{1,2}*Kinv{1,2}*f{1,2}'+B{1,1}*inv(K{1,1})*f{1,1}'+B{1,4}*Kinv{1,4}*f{1,4}';
%             g = B{1,2}*R{1,2};
%             %fi = 2*B{1,2}*Kinv{1,2}*B{1,2}'+B{1,1}*inv(K{1,1})*B{1,1}'+B{1,4}*Kinv{1,4}*B{1,4}';
%             alpha{1,2} = (inv(g'*g)*g')*(F*lambda-di);
%             
%             di = 2*B{1,3}*Kinv{1,3}*f{1,3}'+B{1,1}*inv(K{1,1})*f{1,1}'+B{1,4}*Kinv{1,4}*f{1,4}';
%             g = B{1,3}*R{1,3};
%             %fi = 2*B{1,3}*Kinv{1,3}*B{1,3}'+B{1,1}*inv(K{1,1})*B{1,1}'+B{1,4}*Kinv{1,4}*B{1,4}';
%             alpha{1,3} = (inv(g'*g)*g')*(F*lambda-di);
%             
%             di = B{1,2}*Kinv{1,2}*f{1,2}'+B{1,3}*Kinv{1,3}*f{1,3}'+2*B{1,4}*Kinv{1,4}*f{1,4}';
%             g = B{1,4}*R{1,4};
%             %fi = B{1,2}*Kinv{1,2}*B{1,2}'+B{1,3}*Kinv{1,3}*B{1,3}'+2*B{1,4}*Kinv{1,4}*B{1,4}';
%             alpha{1,4} = (inv(g'*g)*g')*(F*lambda-di);
% 
%             
%             
%             alpha{1,2}
%             alpha{1,3}
%             alpha{1,4}
            
            alphaI = inv(G'*G)*G'*(F*lambda-d);
            alpha{1,1} = 0;
            alpha{1,2} = alphaI(1:3,1);
            alpha{1,3} = alphaI(4:6,1);
            alpha{1,4} = alphaI(7:9,1);
            
%             alpha{1,2}
%             alpha{1,3}
%             alpha{1,4}
            
%             kk
%             w(:,kk-1)
%             w(:,kk)
%             w(:,kk+1)
            
        end
    end
end

