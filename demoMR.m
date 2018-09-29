%%This is a demo for removing outliers. In this demo, the SIFT matches  is known in advance.
clear;
close all;
warning off all;
normalize = 1;
conf.lambda1 = 3;
conf.lambda2 = 0.05;
conf.r = 0.08;
conf.a = 10;
conf.M = 15;
%conf.beta = 0.1

show_result = 1;

load save_chinese_occlusion_6_97.mat;
addpath('sc');

X = y1;
Y = y2a;
tic;

conf = LLT_init(conf);
[X2,Y2,V,ind,ifxmore,out] = sc(X,Y,y2);
%[X2,Y2,V,ind,ifxmore,out] = idsc(X,Y,y2);

if ifxmore
    X2 = X2(ind,:);
end
lenX = size(X,1);
lenY = size(Y,1);
if lenX < lenY
    Xtmp=NaN*ones(lenY,2);
    Xtmp(1:lenX,:)=X;
    X=Xtmp(out,:);
else
    X=X(out,:);
end
%X = X(out,:);
normal.xm=0; normal.ym=0;
normal.xscale=1; normal.yscale=1;

if normalize
    [nX, nY, normal]=norm_ind(X2,Y2,ind);
end

if ~exist('conf'), conf = []; end

VecFld=MR(nX, nY, conf,ind);
Xtemp = [];
means = [];
for iii = 1:30
    if normalize,VecFld.TX=(VecFld.TX)*normal.yscale+repmat(normal.ym,size(VecFld.TX,1),1);end
    Xk = VecFld.TX;
    
   [X2,Y2,V,ind,~,out] = sc(Xk,Y2,V);
   %X2 = Xk;
    if ifxmore
        X2 = X2(ind,:);
    end
    %X = X(out,:);
    YY = Xk(ind,:);
    mea = sum(sqrt(sum((YY-V).^2,2)))/size(YY,1);
    means = [means;mea];
    Xtemp = [Xtemp;Xk];

    if normalize
        [nX, nY, normal]=norm_ind(X2,Y2,ind);
    end
    VecFld=MR(nX, nY, conf,ind);
end

if normalize,VecFld.TX=(VecFld.TX)*normal.yscale+repmat(normal.ym,size(VecFld.TX,1),1);end
Xk = VecFld.TX;
toc

if show_result
    figure
    plot(X(:,1),X(:,2),'b+',Y(:,1),Y(:,2),'ro');
    title( 'original data');
    figure
    plot(Xk(:,1),Xk(:,2),'b+',Y(:,1),Y(:,2),'ro');
    title( 'registration result');
end