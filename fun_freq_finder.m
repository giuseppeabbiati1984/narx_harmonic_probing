function [dt, t] = fun_freq_finder(w, dw)
% This function determines the length of the signal needed to avoid leakage in the DFT computation.

if numel(w) > 3 % check if the supplied frequncies are fewer than 3
    errordlg('Error in number of harmonics sent. Should be 1, 2 or 3')
    return
end

dt = (2*pi)/(2*numel(w)*(sum(w))); % Upper bound computation m*sum(\Omega_p). Line 7. Algorithm 1. (pg. 6)
r =  single((2*pi/dw)/dt); % Length of the time-vector (must be an integer). Line 8+9. Algorithm 1. (pg. 6)
t = 0:dt:(r-1)*dt ; % Time vector, last element removed for periodicity

if ~mod(r, 1) == 0 % Checks if the r is NOT an integer (if not enter if-statement)
    errordlg("Length is not an integer")
    return
end


if isempty(t) % Check for empty vector
    errordlg("Time vector is empty!")
    return
end


end