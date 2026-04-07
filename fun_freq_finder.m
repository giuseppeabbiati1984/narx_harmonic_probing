function [dt, t] = fun_freq_finder(w, dw)
% This function determines the length of the signal needed to avoid leakage in the DFT computation.

if numel(w) > 3
    errordlg('Error in number of harmonics sent. Should be 1, 2 or 3')
    return
end

    dt = (2*pi)/(2*numel(w)*(sum(w))); % Upper bound computation m*sum(\Omega_p)
    t_len =  single((2*pi/dw)/dt); % Length of the time-vector (must be an integer)

if ~mod(t_len,1) == 0													% Checks if the t_len is NOT an integer (if not enter if-statement)
    errordlg("length is not an integer")
    return
end

%time vecotr
t = 0:dt:(t_len-1)*dt ;													% Remove last element to make it periodic

if isempty(t)
	disp("time vector is empty!!!!")
end


end