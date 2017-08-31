function [ cent, vol ] = eval_centnet( t, Vs )
% function to evaluate centroid of intersection polytope
% Author: Steffen Limmer (steffen.limmer@tu-berlin.de)
% Last update: 31.08.2017
% algorithm uses sum of M-times inverse laplace transform of the function 
% F(la(1),...,la(M)) = exp(-a.'*la)/prod(B*la)
% inputs:
% y: affine constraint vector, size M x 1
% V: basis of an M-dim subspace, size M x N
% outputs:
% cent: centroid of the intersection polytope
% vol: volume of the intersection polytope

[N,M] = size(Vs);
if norm(Vs.'*Vs - eye(M),'fro') >= 1e-12
    warning('input Vs should be orthogonal')
end

% compute the volume of terms
S = [eye(N),zeros(N,1)];
% allocate
cent1 = zeros(N,1);
cent2 = zeros(N,N);

for k = 1:N+1
    a1{k} = Vs.'*S(:,k);
    %B1{k} = (S(:,setdiff(1:N+1,k)) - S(:,k)*ones(1,N)).'*Vs;
    B1{k} = (S(:,[1:k-1,k+1:N+1]) - S(:,k)*ones(1,N)).'*Vs;
    % cent1 is similar to volnet
    cent1(k,1) = eval_laplace((t),(a1{k}),(B1{k}),0); 
    for n=[1:k-1,k+1:N+1]
        a2{k,n} = Vs.'*S(:,n);
        %B2{k,n} = [(S(:,setdiff(1:N+1,n)) - S(:,n)*ones(1,N)), S(:,k) - S(:,n)].'*Vs;
        B2{k,n} = [S(:,k) - S(:,n), S(:,k) - S(:,n), (S(:,[1:min(n,k)-1,min(n,k)+1:max(n,k)-1,max(n,k)+1:N+1]) - S(:,n)*ones(1,N-1))].'*Vs;
        cent2(k,n) = eval_laplace((t),(a2{k,n}),(B2{k,n}),1);
    end
end

vol = double( sum(cent1) );
cent = double((cent1(1:end-1) + sum(cent2(1:end-1,:),2) - sum(cent2(:,1:end-1).',2))/sum(cent1));

end