% (D. Stamenov, Aug 30, 2023) This script produces the symbolic functions
% needed for the poly NARX and saves them to a file.
clear;  clc; %close all;

%########### Develop Symbolic Functions ################
Ord = "3" ;
lagNumF = 3 ;
lagNumZ = 4 ;
n = lagNumZ + lagNumF ; % total number of lags

% Input
x_sym = sym('x', [n, 1]) ;

% Compute only up to the order selected
switch Ord
    case "1" %#####################1st ORDER##############################
        % jacobian
        J_sym = sym('J' ,[1,n]) ;

        y_sym = J_sym*x_sym ;
        a_sym = [symvar(J_sym)] ;

    case "2" %#####################2nd ORDER##############################
        % jacobian
        J_sym = sym('J',[1,n]) ;
        % hessian (symmetric)
        H_sym = sym('H', [n, n]) ;
        H_sym = tril(H_sym,-1) + tril(H_sym,0).' ; % enforce symmetric

        y_sym = J_sym*x_sym + x_sym.'*H_sym*x_sym/2 ;

        % vectorization of the coefficients
        a_sym = [symvar(J_sym), symvar(H_sym)] ;

    case "3" %#####################3rd ORDER##############################
        % jacobian
        J_sym = sym('J', [1,n]) ;

        % hessian (symmetric)
        H_sym = sym('H', [n, n]) ;
        H_sym = tril(H_sym, -1) + tril(H_sym, 0).' ;									% Enforce symmetric

        % Tensor Derivative (symmetric)
        T_sym = sym('T', [n, n, n]) ;
        % T_sym = tril(T_sym, -1) + tril(T_sym, 0).' ;									% Enforce symmetric

        for i =1:size(T_sym, 1)
	  S(:, i) = reshape(T_sym(:, :, i), n, []) * x_sym ;
        end
        D3 = 1/6 * x_sym.' * S * x_sym ;
        y_sym = J_sym*x_sym + 1/2 * x_sym.'*H_sym*x_sym + D3 ;

        % vectorization of the coefficients
        a_sym = [symvar(J_sym), symvar(H_sym), symvar(T_sym)] ;

    case "2O" %#####################ONLY 2nd ORDER###########################
        J_sym = 0*sym('J', [1,n]) ;

        H_sym = sym('H', [n, n]) ;
        H_sym = tril(H_sym,-1) + tril(H_sym,0).' ;										% enforce symmetric

        y_sym =  x_sym.'*H_sym*x_sym/2 ;

        % vectorization of the coefficients
        a_sym = [symvar(H_sym)] ;

    case "23O" %#####################ONLY 2nd and 3th ORDER###########################

       % jacobian
        J_sym = 0*sym('J', [1,n]) ;

        % hessian (symmetric)
        H_sym = sym('H', [n, n]) ;
        H_sym = tril(H_sym, -1) + tril(H_sym, 0).' ;									% Enforce symmetric

        % Tensor Derivative (symmetric)
        T_sym = sym('T', [n, n, n]) ;
        % T_sym = tril(T_sym, -1) + tril(T_sym, 0).' ;									% Enforce symmetric

        for i =1:size(T_sym, 1)
	  S(:, i) = reshape(T_sym(:, :, i), n, []) * x_sym ;
        end
        D3 = 1/6 * x_sym.' * S * x_sym ;
        y_sym = J_sym*x_sym + 1/2 * x_sym.'*H_sym*x_sym + D3 ;

        % vectorization of the coefficients
        a_sym = [symvar(J_sym), symvar(H_sym), symvar(T_sym)] ;

    otherwise
        errordlg('Select a viable order (1, 2, 2O, 23O, or 3)')
        return
end


% design coefficients
[p_sym,~] = coeffs(y_sym, a_sym) ;

% vector of regressors
matlabFunction(p_sym, 'file',  ['fun_p_nf' num2str(lagNumF), 'nz' , num2str(lagNumZ), '_ord', num2str(Ord)], 'vars', {x_sym.'})

% jacobian vector (full)
matlabFunction(J_sym, 'file', ['fun_J_nf' num2str(lagNumF), 'nz' , num2str(lagNumZ), '_ord', num2str(Ord)], 'vars', {a_sym})


if Ord == "2" || Ord == "3" || Ord =="2O" || Ord =="23O"
    % hessian matrix (full)
    matlabFunction(H_sym, 'file', ['fun_H_nf' num2str(lagNumF), 'nz' , num2str(lagNumZ), '_ord', num2str(Ord)], 'vars', {a_sym})
end

if Ord == "3" || Ord =="23O"
    % Tensor(full)
    matlabFunction(T_sym, 'file', ['fun_T_nf' num2str(lagNumF), 'nz' , num2str(lagNumZ), '_ord', num2str(Ord)], 'vars', {a_sym})
end




