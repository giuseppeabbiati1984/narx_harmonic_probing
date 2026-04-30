function [w_check] = compute_freq_check(n, w)
%This function produces a list of frequencies that cannot be probed with the numerical probing algorithm.

%INPUT:
% w: freqeuncies invovled in the probing
% n: order of the probing

%OUTPUT:
% w_check: frequencies that cannot be probed

if n==1
    errordlg('All frequencies are identifiable for n=1. Choose n>1.')
    return
end

iq = 0:1:n ; 

iq_cell = repmat({iq}, 1, n); % repeat the iq n times
iq_eq = combinations(iq_cell{:}); % generate combinations of the iq vector n times
iq_eq = iq_eq.Variables ; % convert table to array


iq_eq = iq_eq(sum(iq_eq, 2)==n, :) ; % retain only the rows that sum to n
iq_eq(sum(iq_eq==ones(1,n), 2)==n,:) = [] ;
w_check = [iq_eq*w.'] ; 

end
