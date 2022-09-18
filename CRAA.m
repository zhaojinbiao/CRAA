function [C,max_iter] = CRAA(X,opts,M,lambda_1,K)


V = opts.V;
N = opts.N;
SD = opts.SD;

P = zeros(SD,K);
H = zeros(K,N); Er = zeros(K,N);
Z = zeros(N,N);J = zeros(N,N);C = zeros(N,N);
Eh = zeros(SD,N);Ec = zeros(N,N);
Y1 = zeros(SD,N);Y2 = zeros(K,N);Y3 = zeros(N,N);Y4 = zeros(N,N);

H = rand(K,N);

IsConverge = 0;
mu = 1e-4;
pho = 1.31;
max_mu = 1e6;
max_iter = 100;


iter = 1;

H = H./repmat(sqrt(sum(H.^2,1)),size(H,1),1);

while (IsConverge == 0&&iter<max_iter+1)
    % update P
    P=updateP(Y1,mu,M,H,Eh);
    % update H
    A = mu*((P'*P)+eye(K));
    B = mu*(Z*Z'-Z-Z')+eye(N)*1e-6;
    K_p = (P'*M+Er-Er*Z'-P'*Eh)*mu+(P'*Y1-Y2+Y2*Z');
    
    A(isinf(A))=0;
    A(isnan(A))=0;
    B(isinf(B))=0;
    B(isnan(B))=0;
    K_p(isinf(K_p))=0;
    K_p(isnan(K_p))=0;
    H = sylvester(A,B,K_p);
    
    % update J
    
    J = softth((C-Y4/mu)+eye(N)*1e-8,lambda_1/mu);
    
    % update Z
    Z = (H'*H+eye(N))\((H'*Y2-Y3)/mu+H'*H-H'*Er+C+Ec);
    C = (Y3-Y4)/(2*mu)+(J+Z-Ec)/2;
    
    % update E
    G = [M-P*H+Y1/mu; H-H*Z+Y2/mu; Z-C+Y3/mu];
    [E] = solve_l1l2(G,1/mu);
    Eh = E(1:SD,:);
    Er = E(1+SD:SD+K,:);
    Ec = E(1+SD+K:SD+K+N,:);
    
    % update multipliers
    Y1 = Y1+ mu*(M-P*H-Eh);
    Y2 = Y2+ mu*(H-H*Z-Er);
    Y3 = Y3+ mu*(Z-C-Ec);
    Y4 = Y4+ mu*(J-C);
    
    mu = min(pho*mu, max_mu);
    thrsh = 0.00001;
    % check the convergence conditions
    thrsh = 0.0001;
    norm(M-P*H-Eh,inf);
    norm(H-H*Z-Er,inf);
    norm(J-Z,inf);
    if(norm(M-P*H-Eh,inf)<thrsh && norm(H-H*Z-Er,inf)<thrsh && norm(J-C,inf)<thrsh&&norm(Z-C-Ec,inf)<thrsh)
        IsConverge = 1;
    end
    iter = iter + 1;
end
end