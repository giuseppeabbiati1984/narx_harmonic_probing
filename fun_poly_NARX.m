function NARX = fun_poly_NARX(Xy, Xu, Yy, ord, varargin)
%This function trains a polynomial NARX with a LASSO regression. 
% This function calls up symbolically generated functions from "script_symbolic_poly_gen.m"  

if numel(varargin) > 0
	Lamb_val = varargin{1} ;
else
	Lamb_val = 0.01 ;
end


% lag samples
NARX.lagNumY = size(Xy, 2) ;  % response signal
NARX.lagNumU = size(Xu, 2) ;  % loading signal

%% Train Poly-NARX

% data setup
x = [Xy, Xu] ; %vector of input and output lags
y = Yy ; %vector of outputs

%Create a string with the name of the function based on nf and nz
string_fun_p = ['fun_p_nf' num2str(NARX.lagNumY), 'nz' , num2str(NARX.lagNumU), '_ord', num2str(ord), '(x)'] ;

% evaluate the regressors
p = eval(string_fun_p) ; %coeff. from x combined in front of each term in J and H
% p = fun_p(x) ;


% ls regression solution (non-sparse so all elements of p are non-zero)
a_ls = pinv(p)*y ;

%% LASSO regression

% this should do automatically some sort of cross-validation

rng(6) ; % For reproducibility
if numel(varargin) > 0
    [B, FitInfo] = lasso(p, y, CV=5, Lambda=Lamb_val, Intercept = false) ;						% Specify the penalty value
elseif numel(varargin) == 0
   [B, FitInfo] = lasso(p, y, CV=5, NumLambda=200, Intercept = false, RelTol=1e-4) ;				% Use a range of penalty values
    % [B, FitInfo] = lasso(p, y, CV=5, NumLambda=1000, Intercept = false, RelTol=1e-4, DFmax = 10) ;		% DFmax = number of nonzero terms
    % [B, FitInfo] = lasso(p, y, CV=5, LambdaRatio=0.1, NumLambda=400, Intercept = false, RelTol=1e-4) ;	% LambdaRatio =smallest/largest lambda value in the sequence
    % [B, FitInfo] = lasso(p, y, CV=5, Lambda=linspace(0.001, 1, 1000), Intercept = false, RelTol=1e-4) ;
end

idxLambda1SE = FitInfo.Index1SE;
a_lasso = B(:, idxLambda1SE);
a_lasso_0 = FitInfo.Intercept(idxLambda1SE); % this must be zero as intercept is set to false

nonzero_indices = a_lasso~=0 ;

% sparsity
% sum(nonzero_indices)/numel(nonzero_indices) 
sum(nonzero_indices) 
%% estimation of Jacobian and Hessian for multiple segments (boostrap)
% A NARX model is estimated for each segment considering the model
% structure produced by the LASSO regression on the entire dataset.
% data setup

% x = [Xf, Xz] ; % usde in p_full
% y = Yf ;

% design matrix
p_full = eval(string_fun_p);
p1 = p_full(:, nonzero_indices);

% p1 = fun_p_lasso(x) ;
% p1 = p_full(:, nonzero_indices);


% ls solution (non-sparse so all elements of p are non-zero)
% a_ls1 = @(x,y) pinv(p1(x))*y ;
a_ls1 = pinv(p1)*y ;


% reconstruction of Jacobian and Hessian for probing
a_ls1_ext = zeros(size(a_ls)) ;
a_ls1_ext(nonzero_indices, 1) = a_ls1 ;

NARX.nonzero_indices = nonzero_indices;
NARX.FitInfo = FitInfo;

%Create a function name strings for calling for J and H based on nf and nz
string_fun_J = strcat('fun_J_nf', num2str(NARX.lagNumY), 'nz' , num2str(NARX.lagNumU), '_ord', num2str(ord), "(a_ls1_ext.')") ;
NARX.J(:) =  eval(string_fun_J)  ; % reconstruction of sparse Jacobian

if ord == "2" || ord == "3" || ord =="2O" || ord =="23O"
string_fun_H = strcat('fun_H_nf', num2str(NARX.lagNumY), 'nz' , num2str(NARX.lagNumU), '_ord', num2str(ord), "(a_ls1_ext.')") ;
NARX.H(:, :) = eval(string_fun_H) ; % reconstruction of sparse Hessian
end

if ord == "3" || ord =="23O"
string_fun_T = strcat('fun_T_nf', num2str(NARX.lagNumY), 'nz' , num2str(NARX.lagNumU), '_ord', num2str(ord), "(a_ls1_ext.')") ;
NARX.T(:, :, :) = eval(string_fun_T) ; % reconstruction of sparse Tensor
end


end