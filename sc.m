    function [newX,newY,newV,ind,ifxmore,out] = sc(X,Y,Z)
    mean_dist_global=[]; % use [] to estimate scale from the data
    eps_dum=0.15;
    nbins_theta=12;
    nbins_r=5;
    r=1; % annealing rate
    beta_init=1;
    nsamp1=size(X,1);
    nsamp2=size(Y,1);
    ndum1=0;
    ndum2=0;
    r_inner=1/8;
    r_outer=2;
    k = 1;
    out_vec_1=zeros(1,nsamp1);
    out_vec_2=zeros(1,nsamp2);
    if nsamp2>nsamp1
        ndum1=ndum1+(nsamp2-nsamp1);
    end
    Xk  = X;
    [BH1,mean_dist_1]=sc_compute(Xk',zeros(1,nsamp1),mean_dist_global,nbins_theta,nbins_r,r_inner,r_outer,out_vec_1);
    [BH2,mean_dist_2]=sc_compute(Y',zeros(1,nsamp2),mean_dist_1,nbins_theta,nbins_r,r_inner,r_outer,out_vec_2);
    % compute regularization parameter
    %beta_k=(mean_dist_1^2)*beta_init*r^(k-1);
    % compute pairwise cost between all shape contexts
    costmat=hist_cost_2(BH1,BH2);



    % pad the cost matrix with costs for dummies
    nptsd=nsamp1+ndum1;
    costmat2=eps_dum*ones(nptsd,nptsd);
    costmat2(1:nsamp1,1:nsamp2)=costmat;
    cvec=hungarian(costmat2);

    [a,cvec2]=sort(cvec);
    %out_vec_1=cvec2(1:nsamp1)>nsamp2;%
    %out_vec_2=cvec(1:nsamp2)>nsamp1;

    % format versions of Xk and Y that can be plotted with outliers'
    % correspondences missing
    X2=NaN*ones(nptsd,2);
    X2(1:nsamp1,:)=Xk;
    X2=X2(cvec,:);
    X2b=NaN*ones(nptsd,2);
    X2b(1:nsamp1,:)=X;


    X2b=X2b(cvec,:);
    Y2=NaN*ones(nptsd,2);
    Y2(1:nsamp2,:)=Y;
    
   

    ifxmore = 0;
    if(nsamp1 < nsamp2)
        ind_err = find(isnan(X2b(:,1)));
        all = 1:nptsd;
        ind = setdiff(all, ind_err);
        newY = [Y2(ind,:); Y2(ind_err,:)];
        newX = [X2b(ind,:); X2b(ind_err,:)];
        V=NaN*ones(nptsd,2);
        V(1:min(nsamp1,nsamp2),:)=Z;
        V=V(cvec,:);
        newV = [V(ind,:)];
        ifxmore = 1;
         out = cvec([ind,ind_err']);
    else
        ind_err = find(isnan(Y2(:,1)));
        all = 1:nptsd;
        ind = setdiff(all, ind_err);
        newY = [Y2(ind,:); Y2(ind_err,:)];
         newX = [X2b(ind,:); X2b(ind_err,:)];
        newV = Z;
        out = cvec([ind,ind_err']);
    end
    ind = 1:length(ind);