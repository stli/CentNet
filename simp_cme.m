function [ x_cm, volpoly, nsimp ] = simp_cme( t, Vs, x )
% function to evaluate centroid of intersection polytope via simplicial
% decompotision
% Author: Steffen Limmer (steffen.limmer@tu-berlin.de)
% Last update: 31.08.2017

[N,M] = size(Vs);
% check orthogonality
if norm(Vs.'*Vs - eye(M),'fro') >= 1e-12
    warning('input Vs should be orthogonal')
end

dim_ker = N-M;

% H-polytope: {x | Hx <= b}
% estimator given by x_cm(i) = vol(U)^{-1}*int_{x \in U} x(i) dx
H = [-eye(N);ones(1,N)];
Heq = [Vs.'];
b = [zeros(N,1);1];
beq = [t];


V = bsxfun(@max,lcon2vert(H,b,Vs.',t,1e-12,1).',0);

if isempty(V) || size(V,2)<dim_ker+1 
    % use genie knowledge of interior point
    V = bsxfun(@max,qlcon2vert(x,H,b,Heq,beq,1e-24,1).',0);
end

% normalize
num_verts = size(V,2);

% compute N-linear forms of type (e_i.*x) according to Remark (9) in "how
% to integrate a polynomial over a simplex", by baldoni et. al.
L = eye(N);
if num_verts == 1
    volpoly = 0;
    x_cm = V(:,1);
    nsimp = 0;   
elseif num_verts == 2
    % line segment
    volpoly = norm(V(:,1)-V(:,2));
    x_cm = (V(:,1) + V(:,2))/2;
    nsimp = 1;
elseif num_verts == dim_ker + 1 
    % proper simplex, i.e., simplex with (dim_ker+1)-affinely
    % independent vertices in affine subspace of dimension (dim_ker)   
    d = num_verts - 1;
    x_cm = factorial(d)/factorial(1+d)*sum(L*V,2);
    nsimp = 1;
    linV = V(:,2:end) - repmat(V(:,1),1,size(V,2)-1);
    volpoly = prod(svds(linV))/factorial(size(linV,2));
elseif num_verts > dim_ker + 1 
    % perform simplicial decomposition to obtain disjoint set of size(T,1)
    % proper simplices
    T = [];
    i = 0;
    while isempty(T) & (i<1e3) % loop over different subspace projections
    idxT = randperm(N);
    %try
    T = qhullmx(V(idxT(1:dim_ker),:), 'd ', 'Qt Qbb Qc Qx Qz');
    i = i+1;
    end

    nsimp = size(T,1);
    for i = 1:size(T,1)
        % obtain volumes and individual conditional mean estimates for each
        % simplex
        Vi{i} = V(:,T(i,:));
        linVi{i} = Vi{i}(:,2:end) - repmat(Vi{i}(:,1),1,size(Vi{i},2)-1);
        volVi(1,i) = prod(svds(linVi{i}))/factorial(size(linVi{i},2));
        x_cmi(:,i) = 1/(size(linVi{i},2)+1)*sum(L*Vi{i},2);
    end
    % final cm estimate is convex combination with weights determined by 
    % volumes of individual simplices
    x_cm = sum( x_cmi.*(ones(N,1)*volVi), 2)/sum(volVi);
    volpoly = sum(volVi);
end
    
end

