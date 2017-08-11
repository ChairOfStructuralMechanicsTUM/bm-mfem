classdef Feti2Solver < FetiPreparer
    %UNTITLED Summary of this class goes here
    %   This class can solve system with the FETI level 2 method
    
    properties
    end
    
    methods (Static)
        function [lambda, alpha] = pCG(K,f,B,R,Kinv)
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
                        
            

            

            %Lumped preconditioner
            %F01Inv = B01*K01*B01' + B02*K02*B02';
            %Dirichlet preconditioner
            %Q01 = [zeros(8,11);
            %       zeros(3, 8), K01(9:11, 9:11) - K01(1:8,9:11)'*inv(K01(1:8,1:8))*K01(1:8,9:11)];
            %Q02 = [zeros(11,14);
            %       zeros(3, 11), K02(12:14, 12:14) - K02(1:11,12:14)'*inv(K02(1:11,1:11))*K02(1:11,12:14)];
               
            %F01Inv = B01*Q01*B01' + B02*Q02*B02';
            
            
            %SECOND LEVEL
            %From The second generation
            %C is matrix localizing the jump equations at the corner of the
            %problem. rectangular matrix with a number of columns equal to
            %the number of corners in the problem mulitplied by the number
            %of Lagrange multipliers defined at these corners. number of
            %rows equal to interface dofs (4x4)
            
%             C = [0 0 0 0;
%                  0 1 0 0;
%                  0 0 0 0;
%                  0 0 0 1];
%                 
%               gamma = (C'*P'*F01*P*C)
% 
%             %correct
%             lambda(:,1) = G02*inv(G02'*G02)*(f02*R02)'+P*C*gamma;
%             %w is the residuum
              %correct
%             w(:,1) = P'*((B02*Kplus*B02'+B01*inv(K01)*B01')-F01*lambda(:,1));


%             while sum(abs(w(:,kk))) > num*10^-3
%                 
%                 y(:,kk) = Proj*F01Inv*w(:,kk);
%   
%                 if kk > 1
%                     temp = 0;
%                     for mm = 1:kk-1
%                         temp = temp + ((y(:,kk)'*F01*p2(:,mm))/(p2(:,mm)'*F01*p2(:,mm)))*p2(:,mm);
%                     end
%                     p(:,kk) = y(:,kk)-temp;   
%                 else
%                     p(:,kk) = y(:,kk);
%                 end
%                 
%                 eqn = (C'*Proj'*F01*Proj*C)*beta == -C'*Proj'*F01*p(:,kk);
%                 [gamma1 gamma2 gamma3 gamma4] = solve(eqn,beta);
%                 gamma = [gamma1; gamma2; gamma3; gamma4];
%                 
%                 p2(:,kk) = p(:,kk)+Proj*C*gamma;
%                 
%                 n(:,kk) = (p2(:,kk)'*w(:,kk))/(p(:,kk)'*F01*p(:,kk));
%                 
%                 lambda(:,kk+1) = lambda(:,kk) + n(:,kk)*p2(:,kk);
%                 
%                 w(:,kk+1) = w(:,kk) - n(:,kk)*Proj'*F01*p2(:,kk);
%                 
%                 
%                 %constraint check
%                 g02 = G02'*lambda(:,kk);
%                 test = G02'*p(:,kk);
%                 
%                 kk = kk+1;
%             end

%             lambda = lambda(:,kk);            
%             alpha = (inv(G02'*G02)*G02')*(F01*lambda-B02*Kplus*f02'+B01*inv(K01)*f01');
%             %is this possible, according to "a unified framework..." p.265
%             %eq. 30
%             alpha = (inv(G02'*G02)*G02')*(B02*Kplus*f02'+B01*inv(K01)*f01'-F01*lambda);
            
            
            
        end
    end
end

