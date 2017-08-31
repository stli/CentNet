function [ val ] = eval_laplace( t, a, B, dpole )
% function to evaluate M-times inverse laplace transform of the function
% F(la(1),...,la(M)) = exp(-a.'*la)/prod(B*la)
% Author: Steffen Limmer (steffen.limmer@tu-berlin.de)
% Last update: 31.08.2017
% inputs:
% t: function variable, size M x 1
% a: numerator exponential weight vector, size M x 1, !!a(1)>=0!!
% B: denominator product weight matrix, size N x M
% dpole: flag for double pole, matrix B supposed to have equal 1st & 2nd
% row
% outputs:
% val: inverse laplace transform evaluated for the vector t;
% (global s: only for evaluation purpose)

%if a(1)<0 % optional for evaluation purpose
%    disp('Warning: negative value detected')
%end

if (t(1)-a(1))<=0
    val = 0;
    return;
end

% collect relevant nominator terms, 
a = a(:);
% truncate denominator terms
nnzs = B(:,1)~=0;
if dpole==1 & any( B(1,:)~=B(2,:) )
    disp('error: 1st and 2nd row of B should be equal')
end

Bt = B(nnzs,:);
[N,M] = size(Bt);
Ct = Bt./repmat(Bt(:,1),1,M);

if N>=2 && dpole==1 & all( B(1,:)==B(2,:) ) & any( Bt(1,:)~=Bt(2,:) );
    % no dependence of la(1) in double pole la(1): B*la=0
    dpole=0;
end 

    
if dpole == 0
 
    if length(t)>1 % recursive evaluation
        % 3rd dimension contains child terms
        aout = reshape( (ones(N,1)*a(2:end).' + (Bt(:,2:end)*(t(1)-a(1)))./repmat( Bt(:,1),1,M-1) ).', M-1, 1, N);
       
        for n = 1:N
            Bout{n} = [Ct([1:n-1,n+1:N],2:end) - repmat(Ct(n,2:end),N-1,1) ;B(~nnzs,2:end)];
            %Bout(:,:,n) = [Bt([1:n-1,n+1:N],2:end)./repmat( Bt([1:n-1,n+1:N],1),1,M-1) - repmat(1/Bt(n,1)*Bt(n,2:end),N-1,1) ;B(~nnzs,2:end)];
            tmp(n) = eval_laplace( t(2:end), aout(:,:,n), Bout{n},0 );
        end
        val = heaviside(t(1)-a(1))*real(sum(sort(complex(tmp))))*1/prod(B(nnzs,1),1);
    else
        val = heaviside(t(1)-a(1))*(t(1)-a(1))^(N-1)*1/(factorial(N-1)*prod(B(nnzs,1)));
    end
elseif dpole == 1

    if length(t)>1 % recursive evaluation
        aout = reshape( (ones(N-1,1)*a(2:end).' + (Bt(2:end,2:end)*(t(1)-a(1)))./repmat( Bt(2:end,1),1,M-1) ).', M-1, 1, N-1);
        %B1 = [C(setdiff(1:N-1,1),2:end) - ones(N-2,1)*C(1,2:end);BU(~nnzs,2:end)];
        Bout{1} = [Ct(3:end,2:end) - repmat(Ct(1,2:end),N-2,1);B(~nnzs,2:end)];
        tmp(1) = (t(1)-a(1))*eval_laplace( t(2:end), aout(:,:,1), Bout{1}, 0);
        for n = 2:N-1
            Bout{n} = [[1;1]*(Ct(n+1,2:end)-Ct(1,2:end)); Ct([3:n,n+2:end],2:end) - repmat(Ct(n+1,2:end),N-3,1);B(~nnzs,2:end)];        
            tmp(n) = eval_laplace( t(2:end), aout(:,:,n), Bout{n},1 );
            tmp(n+N-2) = -eval_laplace( t(2:end), aout(:,:,1), Bout{n},1 );
        end
        val = heaviside(t(1)-a(1))*real(sum(sort(complex(tmp))))/prod(B(nnzs,1),1);
        % testing log calculus to improve stability
        %val = heaviside(t(1)-a(1))*sign(real(sum(sort(complex(tmp)))))*prod(sign(B(:,1)))*exp( log( abs(real(sum(sort(complex(tmp))))) ) - sum(log(abs(B(:,1)))));
    else
        val = heaviside(t(1)-a(1))*(t(1)-a(1))^(N-1)/(factorial(N-1)*prod(B(nnzs,1)));
        % testing log calculus to improve stability
        %val = heaviside(t(1)-a(1))*prod(sign(B(:,1)))*exp( sum( [(N-1)*log(t(1)-a(1)), - log(1:N-1), - log(abs(B(:,1))).']) );
    end
else
    disp('error: dpole should be 0 or 1')
end
            
end