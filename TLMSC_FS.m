function [Q] = TLMSC_FS(alpha,gamma,lambda,beta,X,numCls)
% %%  1/2||X-PG||_F^2+Wvj*Tr((G(v)-S(v))(G(v)-S(v))')+alpha*||C||_*-lambda*Tr(Q'S(i)'S(i)Q)
% rng('default');
m = length(X);
n = size(X{1},2);
for v=1:m
   X{v} = X{v}./(repmat(sqrt(sum(X{v}.^2,1)),size(X{v},1),1)+10e-10);
end

M = numCls;
Q=zeros(n,numCls); 

for in=1:m
   G{in}=zeros(M,n);
   S{in}=full(G{in});
   P{in}=zeros(M,n);
   Z{in}=S{in};
   K{in}=zeros(M,n);
end


sX = [M, n, m];
cross_b = beta; %beta
self_b = gamma; %gamma
epson = 1e-6; 
W = cross_b*ones(m) - diag(cross_b*ones(1,m)) + diag(self_b*ones(1,m));

mu=0.001;
max_mu=10e10;
pho_mu=2;

iter=1;
Isconverg = 0;

while(Isconverg == 0)
    fprintf('----processing iter %d--------\n', iter);
    
    for pj=1:m
        %% update P
        tw = X{pj}*G{pj}';
        tw(isnan(tw)) = 0;
        tw(isinf(tw)) = 0;
        [Up,Sp,Vp] = svd(tw,'econ');
        P{pj} = Up*Vp';
         
         
        %% update G
        bwwc = 0;
        for j = 1:m
            bwwc = bwwc + W(pj,j)*S{j};
        end
        bwws = 0;
        for j = 1:m
            bwws = bwws + W(pj,j)*G{j};
        end  
        ts = P{pj}'*P{pj};
        ts1 = P{pj}'*X{pj}+bwwc-bwws;
        ts(isnan(ts))=0;
        ts(isinf(ts))=0;
        G{pj} = pinv(ts)*ts1;
        %% update S
        sum_G = 0;
        for j=1:m
             sum_G = sum_G+S{j};
        end
         bwwc2 = 0;
         bwws2 = 0;
         for j = 1:m
             bwws2 = bwws2 + W(pj,j)*G{j};
         end
         for j = 1:m
             bwwc2 = bwwc2 + W(pj,j)*S{j};
         end
         tc2 = bwws2-bwwc2+mu*Z{pj}-K{pj}+2*lambda*sum_G*Q*Q';
         S{pj} = tc2/(mu+0.0001);
         S{pj} = max(0,S{pj});
         S{pj} = min(G{pj},S{pj});
    end
    
    %% update Z
    S_tensor = cat(3, S{:,:});
    K_tensor = cat(3, K{:,:});
    Z_tensor = cat(3, Z{:,:});
    z = Z_tensor(:);
    k = K_tensor(:);
    s = S_tensor(:);

    [z, objV] = wshrinkObj(s + 1/mu*k,alpha/mu,sX,0,3); 
    Z_tensor = reshape(z, sX);    
    
    %% update K,mu
    k = k + mu*(s - z);  
    mu = min(mu*pho_mu, max_mu);
   %% update H

    
    %% Update Q
     L = 0;
     for j = 1:m
         L = L+S{j};
     end

     [U,sigma,V] = svd(L,'econ');
     Q = V;
    
    %% check the convergence
    X_PG_total = 0;
    S_Z_total = 0;
    for i=1:m
        if (norm(X{i}-P{i}*G{i},'fro')>epson)
            norm_PG = norm(X{i}-P{i}*G{i},'fro');
             X_PG_total=X_PG_total+norm_PG;          
            fprintf('X_PG %7.10f    \n', norm_PG);
            Isconverg = 0;
        end
        
        Z{i} = Z_tensor(:,:,i);
        K_tensor = reshape(k, sX);
        K{i} = K_tensor(:,:,i);
        
        nor = norm(S{i}-Z{i},'fro');
        if (norm(S{i}-Z{i},'fro')>epson)
            norm_S_Z = norm(S{i}-Z{i},'fro');
            S_Z_total=S_Z_total+norm_S_Z;
            fprintf('norm_S_Z %7.10f    \n', norm_S_Z);
            Isconverg = 0;
        end

        
    end
    err = max(X_PG_total,S_Z_total);
    fprintf('——————err_total %7.10f——————    \n', err);
    if (iter>20)
        Isconverg  = 1;
    end
    iter = iter + 1;
    
end
end


