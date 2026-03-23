%% model generation

close all
clear
clc

debug = true ;

% number of lags (both for f and z)
n = 40 ;

% input
x_sym = sym('x',[n,1]) ;

% jacobian
J_sym = sym('J',[1,n]) ;

% hessian (symmetric)
H_sym = sym('H',[n,n]) ;
H_sym = tril(H_sym,-1) + tril(H_sym,0).' ;

% output
y_sym = J_sym*x_sym + x_sym.'*H_sym*x_sym/2 ;

% vectorization of the coefficients
a_sym = [symvar(J_sym), symvar(H_sym)] ;

% design coefficients
[p_sym,~] = coeffs(y_sym,a_sym) ;

% vector of regressors
matlabFunction(p_sym,'file','fun_p','vars',{x_sym.'})

% jacobian vector (full)
matlabFunction(J_sym,'file','fun_J','vars',{a_sym})

% hessian matrix (full)
matlabFunction(H_sym,'file','fun_H','vars',{a_sym})

%% data analysis

% load data4GA ;
% t1 = data4GA(1,:) ; dt1 = (t1(2)-t1(1)) ;
% z1 = data4GA(2,:) ;
% f1 = data4GA(3,:) ;

% 10000s, 5m
data = importfile10000s5m() ;

t1 = data(:,1) ; dt1 = (t1(2)-t1(1)) ;
z1 = data(:,2) ;
f1 = data(:,3) + data(:,4) ; % first and second order loading

t1 = t1.'; z1 = z1.'; f1 = f1.';
stdZ1 = std(z1) ;
stdF1 = std(f1) ;

z1 = z1/stdZ1 ;
f1 = f1/stdF1 ;

dt2 = 1.00 ;

[Xt,Xz,Xf,Yf,t3,z3,f3] = fun_dataPacking_ga(t1, z1, f1, dt2, n/2, numSegments=1, overlapping=0.0) ;

%%

if debug

    close all

    figure
    hold all
    plot(t1,z1,'--k')
    for i = 1:1:size(t3,1)
        plot(t3(i,:),z3(i,:),'-.')
    end

    figure
    hold all
    plot(t3(1,:),z3(1,:),'--k')
    for j = 1:20:size(Xz{1},1)
        plot(Xt{1}(j,:),Xz{1}(j,:),'-.')
    end

end

%% pre-processing

% data setup
x = [Xf{1}, Xz{1}] ;
y = Yf{1} ;

% evaluate the regressors
p = fun_p(x) ;

%% LS regression

% ls solution (non-sparse so all elements of p are non-zero)
a_ls = pinv(p)*y ;

% validation plot
figure
plot(y,p*a_ls,'x')
% axes('yscale','log')
% plot(sort(abs(p)),'-o')

%% LASSO regression

% this should do automatically some sort of cross-validation

rng default % For reproducibility
[B,FitInfo] = lasso(p,y,CV=4,NumLambda=200,Intercept=false);

idxLambda1SE = FitInfo.Index1SE;
a_lasso = B(:,idxLambda1SE);
a_lasso_0 = FitInfo.Intercept(idxLambda1SE); % this must be zero as intercept is set to false

lassoPlot(B,FitInfo,'PlotType','CV');
legend('show') % Show legend

figure
hold all
plot(y,p*a_ls,'x')
plot(y,p*a_lasso+a_lasso_0,'x')

% indicies of nonzero coefficients to use in the symbolic reconstruction of
% the J and H so that one can construct both matrices using only nonzero
% coefficients.
nonzero_indices = a_lasso~=0 ;

% sparsity
sum(nonzero_indices)/numel(nonzero_indices)

% here I generate a regressor that contains only the term retained by Lasso
matlabFunction(p_sym(nonzero_indices),'file','fun_p_lasso','vars',{x_sym.'})

%% estimation of Jacobian and Hessian for multiple segments (boostrap)

[Xt, Xz, Xf, Yf, t3, z3, f3] = fun_dataPacking_ga(t1, z1, f1, dt2, n/2, numSegments=1, overlapping=0.4) ;

% a NARX model is estimated for each segment considering the model
% structure produced by the LASSO regression on the entire dataset.
for i = 1:1:numel(Xf)

    % data setup
    x = [Xf{i},Xz{i}] ;
    y = Yf{i} ;

    % design matrix
    p1 = fun_p_lasso(x) ;

    % ls solution (non-sparse so all elements of p are non-zero)
    a_ls1 = pinv(p1)*y ;

    % validation plot
    figure
    hold all

    for j = 1:1:numel(Xf)

        % data setup
        x = [Xf{j},Xz{j}] ;
        y = Yf{j} ;

        % regressors on another segment
        p1 = fun_p_lasso(x) ;

        % NRSME from Worden
        errors(i,j) = 100 * sum((y-p1*a_ls1).^2)/(numel(y)*std(y)^2) ;

        % add the last plot
        plot(y,p1*a_ls1,'x')

    end

    pause
    delete(gcf)

    % reconstruction of Jacobian and Hessian for probing
    a_ls1_ext = zeros(size(a_ls)) ;
    a_ls1_ext(nonzero_indices,1) = a_ls1 ;

    % reconstruction of sparse Jacobian and Hessian
    J(i,:) = fun_J(a_ls1_ext.') ;
    H(:,:,i) = fun_H(a_ls1_ext.') ;

end
%  save polyNARXdata J H dt2 stdF1 stdZ1
%%

% for i = 1:1:size(H,1)
%     for j = 1:1:size(H,2)
%         H_cov(i,j) = std(squeeze(H(i,j,:)))/abs(mean(squeeze(H(i,j,:)))) ;
%     end
% end
% 
% figure
% stem(std(a_ls1_vec)./abs(mean(a_ls1_vec)))

%% generalization to arbitrary order NARX

close all
clear
clc

debug = true ;

% number of lags (both for f and z)
n = 4 ;

% input
x_sym = sym('x',[n,1]) ;

% test (matlab shall be smart enough to understand that the terms are identical)
unique([x_sym(1)*x_sym(2)*x_sym(3),x_sym(1)*x_sym(3)*x_sym(2)])

% list of monomials of linear expansion
h1 = x_sym.' ;

% list of monomials included in a 2rd order expansion
h2 = [] ;
for i = 1:1:numel(x_sym)
    for j = 1:1:numel(x_sym)
        h2 = unique([h2,x_sym(i)*x_sym(j)]) ;
    end
end

% list of monomials included in a 3rd order expansion
h3 = [] ;
for i = 1:1:numel(x_sym)
    for j = 1:1:numel(x_sym)
        for k = 1:1:numel(x_sym)
            h3 = unique([h3,x_sym(i)*x_sym(j)*x_sym(k)]) ;
        end
    end
end

a1 = sym('a1',size(h1)) ;
a2 = sym('a2',size(h2)) ;
a3 = sym('a3',size(h3)) ;

h = [h1,h2,h3] ;
a = [a1,a2,a3] ;

% polynarx
y = sum(a.*h) ;

%% more optimized implementation

close all
clear
clc

x = sym('x',[1,3]) ;

p1 = x ;
p2 = [] ;
for i = 1:numel(x)
    for j = 1:numel(x)
        p2 = unique([p2,x(i)*x(j)]) ;
    end
end
p2_mod = [] ;
for i = 1:numel(x)
    for j = i:numel(x)
        p2_mod = [p2_mod,x(i)*x(j)] ;
    end
end

p3 = [] ;
for i = 1:numel(x)
    for j = 1:numel(x)
        for k = 1:numel(x)
            p3 = unique([p3,x(i)*x(j)*x(k)]) ;
        end
    end
end

p3_mod = [] ;
for i = 1:numel(x)
    for j = i:numel(x)
        for k = j:numel(x)
            p3_mod = [p3_mod,x(i)*x(j)*x(k)] ;
        end
    end
end

p4 = [] ;
for i = 1:numel(x)
    for j = 1:numel(x)
        for k = 1:numel(x)
            for l = 1:numel(x)
                p4 = unique([p4,x(i)*x(j)*x(k)*x(l)]) ;
            end
        end
    end
end

p4_mod = [] ;
for i = 1:numel(x)
    for j = i:numel(x)
        for k = j:numel(x)
            for l = k:numel(x)
                p4_mod = [p4_mod,x(i)*x(j)*x(k)*x(l)] ;
            end
        end
    end
end

c1 = sym('c1',size(p1)) ;
c2 = sym('c2',size(p2)) ;
c3 = sym('c3',size(p3)) ;

% polynarx
y = sum(c1.*p1) + sum(c2.*p2) + sum(c3.*p3) ;

% monominal coefficients
c = [c1,c2,c3] ;

% vector of regressors
p = jacobian(y,c) ;

% here we can check the tensor expression of the taylor expansion
diff(diff(y,x(1)),x(2),x(3))
diff(diff(y,x(1)),x(3),x(2))


for i = 1:1:numel(x)
    h1(i) = subs(diff(y,x(i)),x,zeros(size(x))) ;
end

for i = 1:1:numel(x)
    for j = 1:1:numel(x)
        h2(i,j) = subs(diff(diff(y,x(i)),x(j)),x,zeros(size(x))) ;
    end
end

for i = 1:1:numel(x)
    for j = 1:1:numel(x)
        for k = 1:1:numel(x)
            h3(i,j,k) = subs(diff(diff(diff(y,x(i)),x(j)),x(k)),x,zeros(size(x))) ;
        end
    end
end

% here I perform the taylor series expansion to check that our expression
% is correct

y_taylor_expansion = 0 ;
% linear
for i = 1:1:numel(h1)
    y_taylor_expansion = y_taylor_expansion + h1(i) * x(i) ;
end
% quadratic
for i = 1:1:size(h2,1)
    for j = 1:1:size(h2,2)
        y_taylor_expansion = y_taylor_expansion + h2(i,j) * x(i) * x(j) / 2 ;
    end
end
% cubic
for i = 1:1:size(h3,1)
    for j = 1:1:size(h3,2)
        for k = 1:1:size(h3,3)
            y_taylor_expansion = y_taylor_expansion + h3(i,j,k) * x(i) * x(j) * x(k) / 6 ;
        end
    end
end

% this should be one as polynax coincides with its taylor series expansion
test = simplify(y/y_taylor_expansion)

%%

close all
clear
clc

% this is the key to simplify the expressions of the volterra series
% expansion:
% https://en.wikipedia.org/wiki/Multinomial_theorem

syms H2(w1,w2)
assume(H2(w1,w2)==H2(w2,w1)) % this is true only for w1 and w2 and not for any arbitrary w1, w2
simplify(H2(w1,w2)+H2(w2,w1))

w = sym('w',[3,1]) ;

assume(H2(w(1),w(2))==H2(w(2),w(1))) % this is true only w(1), w(2)
simplify(H2(w(1),w(2))+H2(w(2),w(1)))


H2_exp1 = 0 ;

for i = 1:numel(w)
    for j = 1:numel(w)
        H2_exp1 = H2_exp1 + H2(w(i),w(j)) ;
    end
end

% this is still not correct
H2_exp2 = 0 ;
for i = 1:numel(w)
    for j = i:numel(w)
        n = 2 ;
        for k1 = 0:1:n
            for k2 = 0:1:n
                if k1+k2 == n
                    H2_exp2 = H2_exp2 + multinomial(n,[k1,k2])*H2(w(i),w(j)) ;
                end
            end
        end
    end
end


function coeff = multinomial(n,k)
    % This function computes the multinomial coefficient
    % n is the total number of items
    % k is a vector of the counts in each group

    if sum(k) ~= n
        error('Sum of elements in k must equal n');
    end

    coeff = factorial(n);
    for i = 1:length(k)
        coeff = coeff / factorial(k(i));
    end
end


