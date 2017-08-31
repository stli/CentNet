function [ vol ] = eval_volnet( t, Vs )
% evaluate vol( simp \sup Vx=t ), i.e. the volume network 
% algorithm uses sum of M-times inverse laplace transform of the function 
% F(la(1),...,la(M)) = exp(-a.'*la)/prod(B*la)
% inputs:
% t: affine constraint vector, size M x 1
% Vs: (sparse) basis of an M-dim subspace, size N x M

% val: volume of the intersection (simp \sup Vs.'*x=t) 
% (w.r.t to M-dim Lebesbue measure) 
% vol is M-times inverse Laplace transform of
% \int_{\simp} exp(-la(:).'*Vs.'*x) with

[N,M] = size(Vs);
if norm(Vs.'*Vs - eye(M),'fro') >= 1e-12
    warning('input Vs should be orthogonal')
end

S = [eye(N),zeros(N,1)];
for n = 1:N+1
    a{n} = Vs.'*S(:,n);
    %B{n} = (S(:,setdiff(1:N+1,n)) - S(:,n)*ones(1,N)).'*Vs;
    B{n} = (S(:,[1:n-1,n+1:N+1]) - S(:,n)*ones(1,N)).'*Vs;
    voln(n) = eval_laplace(t,a{n},B{n},0);
end

% real(sum(complex())) increases numerical stability
vol = sum( voln ); 
end