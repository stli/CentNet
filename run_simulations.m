% script to simulate Laplace neural networks
% Author: Steffen Limmer (steffen.limmer@tu-berlin.de)
% Last update: 31.08.2017
clear

%% set problem parameters
Ni = 1e2;
N = 9;
M = 3;

%% load input examples
x = abs( gen_vec(N,1,Ni,0) );

%% load/generate matrices and measurements
A = randn(M,N);
try % symbolic svd to increase numerical stability
    [Us,Ss,Vs] = svd(sym(A),'econ');
catch
    [Us,Ss,Vs] = svd(A,'econ');
end
% use householder reflectors to sparsify Vs
% faster but usually less numerically stable
%[Q,R] = qr(sym(Vs.'));
%Vs = (Vs*Q);

% compute measurements
Vs(abs(double(Vs))<1e-16) = 0;
t = double(Vs.'*x);
Vs = double(Vs);

tic  
for i=1:Ni
   [ x_smp(:,i), vol_smp(i), ~ ] = simp_cme( t(:,i), Vs, x(:,i) );

   [ x_dnn(:,i), vol_dnn(i) ] = eval_centnet( t(:,i), Vs );
end
figure
semilogy( sort(abs(vol_dnn - vol_smp)) )
figure
semilogy( sort(sum( (x - x_dnn).^2, 1) ) )
hold on
semilogy( sort(sum( (x - x_smp).^2, 1) ) )

toc
%save results