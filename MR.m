function VecFld = MR(X, Y, conf,ind)
%% Initialization
gamma = conf.gamma;
beta = conf.beta;
lambda1 = conf.lambda1;
lambda2 = conf.lambda2;
theta = conf.theta;
a = conf.a;
MaxIter = conf.MaxIter;
ecr = conf.ecr;
M = conf.M;
r = conf.r;

%% Initialization gamma, C = 0,P = I_N*N
[N, D]=size(X); 
L = length(ind);
Xc = X(ind,:);
Yc = Y(ind,:);
tmp_X = unique(X, 'rows'); idx = randperm(size(tmp_X,1)); 
idx = idx(1:min(M,size(tmp_X,1)));ctrl_pts=tmp_X(idx,:);
if(size(tmp_X,1) < M || size(tmp_X,1) == M)
    ctrl_pts = X;
else
    ctrl_pts=tmp_X(idx,:);
end
K=con_K(ctrl_pts,ctrl_pts,beta);% M*M
U = con_K(X, ctrl_pts, beta);%N*M
Uc = U(ind,:);%L*M
C = zeros(M, D);%M*2
P = eye(L);%L*L
iter=1;  tecr=1; E=1;

sigma2=sum(sum((Yc-Xc).^2))/(L*D);

%%
%compute the graph laplacian matrix A 
X2= sum(X.^2,2);   %N*1
distance = repmat(X2,1,N)+repmat(X2',N,1)-2*X*X';
index = find(distance(:) < r);
W = zeros(N*N,1);
W(index) = exp(-distance(index)/r);
W = reshape(W,N,N);
Dia = sum(W, 2);
A = diag(Dia) - W; %N*N


%%
% EM iteration
%figure;
PS = [];
EE = [];
while ( (iter<MaxIter) & (abs(tecr) > ecr) & (sigma2 > 1e-8) )
    %% E-step.
    % Update P
    E_old = E;
    %Tx = X + U*C;
    Tx = X + U*C; %N*M
    %Tx = U*C;
    Tc = Tx(ind,:);
    
    [P1, E] = get_P(Yc, Tc, sigma2 ,gamma, a);  
    PS = [PS, P1];
   % P1 = max(P1, minP);
    %DrawArrow2(X,Y,P1);
    P = diag(P1);
    %update E 
    %E = E + lambda1/2*trace(C'*K*C)+lambda2/2*trace(Tx'*A*Tx);
    E1 = lambda1/2*trace(C'*K*C);
    E2 = lambda2/2*trace(Tx'*A*Tx);
    E = E + E1+E2;
    tecr=(E-E_old)/E;
    EE = [EE;E1,E2];
    
    %% M-step.
    % update C by solving linear system
    % update sigma^2 and gamma by tr(v_TPv)/D.tr(P)  and tr(P)/N

    %C = (U'*P*U + eta*sigma2*U'*Q*U)\(U'*P*Y-U'*(P + eta*sigma2*Q)*X);%M*2
    % C = (U'*(P + eta*sigma2*Q)*U) \ (U'*P*Y-U'*(P + eta*sigma2*Q)*X); 
    C = (Uc'*P*Uc + lambda1*sigma2*K + lambda2*U'*A*U)\(Uc'*P*(Yc-Xc));
    Vc = Yc - Xc - Uc*C;
    %Vc = Yc - Uc*C;
    sigma2 = trace(Vc'*P*Vc)/(D*trace(P));
    numcorr = length(find(P > theta));
    gamma=numcorr/N;
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
    iter=iter+1;
end
%%
VecFld.X = X;
VecFld.Y = Y;
VecFld.beta = beta;
VecFld.TX=  X +  U*C;
VecFld.C=C;
VecFld.ctrl_pts = ctrl_pts;
VecFld.P = diag(P);
VecFld.PS = PS;
VecFld.EE = EE;
VecFld.VFCIndex = find(VecFld.P > theta);

% disp('Removing outliers succesfully completed.');


%%%%%%%%%%%%%%%%%%%%%%%%
function K=con_K(x,y,beta)
% CON_K constructs the kernel K, 
%   where K(i, j) = k(x, y) = exp(-beta*||x-y||^2).

[n, d]=size(x); [m, d]=size(y);

K=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);
K=squeeze(sum(K.^2,2));
K=-beta * K;
K=exp(K);

%%%%%%%%%%%%%%%%%%%%%%%%
function [P, E]=get_P(Y, Tx, sigma2 ,gamma, a)
% GET_P estimates the posterior probability and part of the energy.

D = size(Y, 2);
temp1 = exp(-sum((Y-Tx).^2,2)/(2*sigma2));
temp2 = (2*pi*sigma2)^(D/2)*(1-gamma)/(gamma*a);
P=temp1./(temp1+temp2);
E=P'*sum((Y-Tx).^2,2)/(2*sigma2)+sum(P)*log(sigma2)*D/2 - log(gamma)*sum(P) - log(1-gamma)*sum(1-P);
