function [dt, t] = fun_freq_finder_V2(w, dw)

% List of frequencies produced by a quadratic interaction with
% bi-chromatice wave

if numel(w) == 2
    % rational frequencies
    w_comb = [w(1), 2*w(1), 3*w(1), 4*w(1), ...
        w(2), 2*w(2), 3*w(2), 4*w(2), ...
        w(1) + w(2), 2*w(1) + 2*w(2), ...
        2*w(1) + w(2), 3*w(1) + w(2), ...
        w(1) + 2*w(2), w(1) + 3*w(2), ...
        3*w(1) + 2*w(2), 2*w(1) + 3*w(2)] ;
elseif numel(w) == 3

w1 = w(1); w2 = w(2); w3 = w(3);

    w_comb = [w1, 2*w1, w2, 3*w1, 2*w2, w3, 3*w2, 2*w3, 3*w3,...
        (w1+w2), (w1+2*w2), (w1+w3), (2*w1+w2), (w1+3*w2), ...
        (w1+2*w3), (2*w1+2*w2), (2*w1+w3), (w2+3*w1), (w2+w3), ...
        (w1+3*w3), (2*w1+3*w2), (2*w1+2*w3), (w2+2*w3), (3*w1+2*w2),...
        (3*w1+w3), (2*w2+w3), (2*w1+3*w3), (w2+3*w3), (3*w1+3*w2), (3*w1+2*w3),...
        (2*w2+2*w3), (w3+3*w2), (3*w1+3*w3), (2*w2+3*w3), (3*w2+2*w3), (3*w2+3*w3),....
        (w1+w2+w3), (w1+w2+2*w3), (w1+2*w2+w3), (2*w1+w2+w3), (w1+w2+3*w3), ...
        (w1+2*w2+2*w3), (w1+3*w2+w3), (2*w1+w2+2*w3), (2*w1+2*w2+w3), (3*w1+w2+w3),...
        (w1+2*w2+3*w3), (w1+3*w2+2*w3), (2*w1+w2+3*w3), (2*w1+2*w2+2*w3), (2*w1+3*w2+w3), ...
        (3*w1+w2+2*w3), (3*w1+2*w2+w3), (w1+3*w2+3*w3), (2*w1+2*w2+3*w3), (2*w1+3*w2+2*w3), ...
        (3*w1+w2+3*w3), (3*w1+2*w2+2*w3), (3*w1+3*w2+w3), (2*w1+3*w2+3*w3), (3*w1+2*w2+3*w3),...
        (3*w1+3*w2+2*w3), (3*w1+3*w2+3*w3)] ;

else
        errordlg('Error in number of harmonics sent. Should be 2 or 3')
        return
end


dt = (2*pi)/(2*numel(w)*(sum(w)))  ;											% Upper bound computation m*sum(\Omega_p)
% t_len = single(abs(2*numel(w)*(sum(w)))/dw);									% Length of the time-vector (must be an integer)
t_len =  single((2*pi/dw)/dt) ;


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