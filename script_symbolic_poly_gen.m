% (D. Stamenov, Aug 30, 2023) This script produces the symbolic functions
% needed for the poly NARX and saves them to a file.
clear;  clc; %close all;

%########### Develop Symbolic Functions ################
Ord = "3" ;
ry = 3 ;
ru = 3 ; % the additional lag is for the present term of the input.
r = ru + ry + 1; % total number of lags

% Input
x_sym = sym('x', [r, 1]) ;

% Compute only up to the order selected
switch Ord
    case "1" %#####################1st ORDER##############################
        % jacobian
        J_sym = sym('J' ,[1,r]) ;

        y_sym = J_sym*x_sym ;
        a_sym = [symvar(J_sym)] ;

    case "2" %#####################2nd ORDER##############################
        % jacobian
        J_sym = sym('J',[1,r]) ;
        % hessian (symmetric)
        H_sym = sym('H', [r, r]) ;
        H_sym = tril(H_sym,-1) + tril(H_sym,0).' ; % enforce symmetric

        y_sym = J_sym*x_sym + x_sym.'*H_sym*x_sym/2 ;

        % vectorization of the coefficients
        a_sym = [symvar(J_sym), symvar(H_sym)] ;

    case "3" %#####################3rd ORDER##############################
        % jacobian
        J_sym = sym('J', [1,r]) ;

        % hessian (symmetric)
        H_sym = sym('H', [r, r]) ;
        H_sym = tril(H_sym, -1) + tril(H_sym, 0).' ;									% Enforce symmetric

        % Tensor Derivative (symmetric)
        T_sym = sym('T', [r, r, r]) ;
        % T_sym = tril(T_sym, -1) + tril(T_sym, 0).' ;									% Enforce symmetric

        for i =1:size(T_sym, 1)
	  S(:, i) = reshape(T_sym(:, :, i), r, []) * x_sym ;
        end
        D3 = 1/6 * x_sym.' * S * x_sym ;
        y_sym = J_sym*x_sym + 1/2 * x_sym.'*H_sym*x_sym + D3 ;

        % vectorization of the coefficients
        a_sym = [symvar(J_sym), symvar(H_sym), symvar(T_sym)] ;

    case "2O" %#####################ONLY 2nd ORDER###########################
        J_sym = 0*sym('J', [1,r]) ;

        H_sym = sym('H', [r, r]) ;
        H_sym = tril(H_sym,-1) + tril(H_sym,0).' ;										% enforce symmetric

        y_sym =  x_sym.'*H_sym*x_sym/2 ;

        % vectorization of the coefficients
        a_sym = [symvar(H_sym)] ;

    case "23O" %#####################ONLY 2nd and 3th ORDER#####################

       % jacobian
        J_sym = 0*sym('J', [1,r]) ;

        % hessian (symmetric)
        H_sym = sym('H', [r, r]) ;
        H_sym = tril(H_sym, -1) + tril(H_sym, 0).' ;									% Enforce symmetric

        % Tensor Derivative (symmetric)
        T_sym = sym('T', [r, r, r]) ;
        % T_sym = tril(T_sym, -1) + tril(T_sym, 0).' ;									% Enforce symmetric

        for i =1:size(T_sym, 1)
	  S(:, i) = reshape(T_sym(:, :, i), r, []) * x_sym ;
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
matlabFunction(p_sym, 'file',  ['fun_p_ry' num2str(ry), 'ru' , num2str(ru), '_ord', num2str(Ord)], 'vars', {x_sym.'})

% jacobian vector (full)
matlabFunction(J_sym, 'file', ['fun_J_ry' num2str(ry), 'ru' , num2str(ru), '_ord', num2str(Ord)], 'vars', {a_sym})


if Ord == "2" || Ord == "3" || Ord =="2O" || Ord =="23O"
    % hessian matrix (full)
    matlabFunction(H_sym, 'file', ['fun_H_ry' num2str(ry), 'ru' , num2str(ru), '_ord', num2str(Ord)], 'vars', {a_sym})
end

if Ord == "3" || Ord =="23O"
    % Tensor(full)
    matlabFunction(T_sym, 'file', ['fun_T_ry' num2str(ry), 'ru' , num2str(ru), '_ord', num2str(Ord)], 'vars', {a_sym})
end




