function [H1_theor, H2_theor, H3_theor, w_dbl] = fun_analytical_duffing(odePars, w)

% This function evaluates the transfer functions for the Duffing oscillator and saves them in arrays.
% The results from this script are used only for comparison.
% % Theoretical Equations (Worden and Tomlison, Sec 8.3, pg 408)

% odePars: structure holding the parameters of the duffing oscillator
% w: single-sided frequency axis



H1_eqn = @(w) 1 / ( -odePars.m * w.^2 + 1i * odePars.c * w + odePars.k1 ) ; %duffing
H2_eqn = @(w1, w2) -odePars.k2 * H1_eqn(w1) * H1_eqn(w2) * H1_eqn(w1+w2);
H3_eqn = @(w1, w2, w3) -1/6*H1_eqn(w1+w2+w3)*...
                    (...
                    4*odePars.k2*(H1_eqn(w1)*H2_eqn(w2,w3)  + H1_eqn(w2)*H2_eqn(w3,w1) + H1_eqn(w3)*H2_eqn(w1,w2)) + ...
                    6*odePars.k3*H1_eqn(w1)*H1_eqn(w2)*H1_eqn(w3) ...
                    ) ;

% TF exact computation:

w_dbl = [-fliplr(w(2:end)), w] ; % Construct the negative part of the freq. axis (avoid doubling of 0)

% Declare arrays in memory
H1_theor = zeros(size(w_dbl)); 
H2_theor = zeros(numel(w_dbl), numel(w_dbl)); 
H3_theor = zeros(numel(w_dbl), numel(w_dbl), numel(w_dbl));

for mm = 1:1:numel(w_dbl)
    H1_theor(mm) = H1_eqn(w_dbl(mm)) ; %LTF
    for nn = 1:1:numel(w_dbl)
        H2_theor(mm, nn) = H2_eqn(w_dbl(mm), w_dbl(nn)) ; %QTF
        for kk = 1:1:numel(w_dbl)
            H3_theor(mm, nn, kk) = H3_eqn(w_dbl(mm), w_dbl(nn), w_dbl(kk)) ; %CTF
        end 
    end
    disp(strcat(num2str(mm), " of ", num2str(numel(w_dbl))) )
end


end 