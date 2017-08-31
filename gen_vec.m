function  X=gen_vec(n,p,Nsamples,CPLX);
%
% function  X=gen_vec(n,p,Nsamples,CPLX);
%
% generate Nsamples n-dimensional vectors uniformly distributed 
% in the unit radius l_p norm-ball.
% X is a (n,Nsamples) matrix.
%
% CPLX=1 : complex vectors
% CPLX=0 : real vectors
%
% Reference:
% G. Calafiore, F. Dabbene, R. Tempo. 
% ``Uniform Sample Generation in $\ell_p$ Balls for Probabilistic Robustness Analysis.'' 
% Proceedings of CDC 1998, Tampa, Florida, Dec. 1998.

% $Id: gen_vec.m 6 2006-10-06 11:57:36Z atremba $

if ~exist('CPLX'), CPLX=1; end;
if ~exist('Nsamples'), Nsamples=1; end;
X=zeros(n,Nsamples);

% complex case

if CPLX==1
   switch p
   case 2
      y=gen_vec(2*n,2,Nsamples,0);
      X=y(1:n,:)+sqrt(-1)*y(n+1:2*n,:);
   case inf
      for k=1:n
         X(k,:)=gen_vec(1,2,Nsamples,1);
      end;         
   otherwise
      for k=1:Nsamples
         j=sqrt(-1);
         x=exp(j*(rand(n,1)*2*pi));
         x=x.*gamrnd(2/p,1,n,1).^(1/p);
         w=rand(1)^(1/(2*n));
         X(:,k)=w*x/norm(x,p);
      end;
   end;
   
else
   
   % real case
   switch p
   case 2
      X=randn(n,Nsamples);
      X=X.*( ones(n,1)*(rand(1,Nsamples).^(1/n)./sqrt(sum(X.^2)) ) );      
   case inf
      X=2*rand(n,Nsamples)-1;
   otherwise
      for k=1:Nsamples
         x=gamrnd(1/p,1,n,1).^(1/p);
         s=rand(n,1)>0.5;       % generates random signs
         x=x.*(-ones(n,1)).^s;
         w=rand(1)^(1/(n));
         X(:,k)=w*x/norm(x,p);
      end;
   end;
end;
