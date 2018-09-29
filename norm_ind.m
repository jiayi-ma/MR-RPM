function  [X, Y, normal] =norm_ind(x,y,ind)
% NORM2 nomalizes the data to have zero means and unit covariance

x = double(x);
y = double(y);

n=size(x,1);
m=size(y,1);

normal.xm=mean(x(ind,:));
normal.ym=mean(y(ind,:));
l = length(ind);

x=x-repmat(normal.xm,n,1);
y=y-repmat(normal.ym,m,1);

normal.xscale=sqrt(sum(sum(x(ind,:).^2,2))/l);
normal.yscale=sqrt(sum(sum(y(ind,:).^2,2))/l);

X=x/normal.xscale;
Y=y/normal.yscale;