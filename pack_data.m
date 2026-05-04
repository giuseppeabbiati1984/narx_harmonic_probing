function [Xy, Xu, Y, t_ss, u_ss, y_ss] = pack_data(tr, ur, yr, dts, numSegments, lenSegment, ru, ry)
% This function packages the data into arrays sutable for LASSO regression.

% INPUT:
% ry: Lags on output signal
% ru: Lags on input signal
% numSegments: number of segments for the data to be partitioned in
% lenSegment: length of each segment
% dts: time-step for the sub-sampling
% tr: time-vector [s]
% ur: input-vector, raw/unsampled
% yr: output-vector, raw/unsampled

% OUTPUT:
% Xy: cell array of lagged output data {}[ lenSeg-ry X ry]
% Xu: cell array of lagged input data {}[ lenSeg-ru X ru+1]
% Y: cell array of present-step output data {}[ lenSeg-ry X 1]
% t_ss: sub-sampled time-vector
% u_ss: sub-sampled input-vector
% y_ss: sub-sampled output-vector


%Compute scaling factors (input and output):
stdU = std(ur) ; % input std. dev. for scaling
stdY = std(yr) ; % output std. dev. for scaling


% Scale the signals by the std. dev.
ur_scl = ur/stdU ; % Scaled input
yr_scl = yr/stdY ; % Scaled output


%ensure column vectors
tr = reshape(tr,[length(tr), 1]);
ur_scl = reshape(ur_scl,[length(ur_scl), 1]);
yr_scl = reshape(yr_scl,[length(yr_scl), 1]);

dt = tr(2)-tr(1);

SS = floor(dts/dt); %integer subsampling rate
% SS = dts/dt; %integer subsampling rate

if max(SS)*lenSegment*numSegments > numel(tr)
    error('Number of segments or length of segment too large')
end


% Construct a causal filter
[b, a] = butter(3, 1/SS, 'low') ;
ff = filter(b, a, yr_scl) ; % filtered force signal (time domain)
zf = filter(b, a, ur_scl) ; % filtered wave elevation signal (time domain)

for ii = 1:numel(SS)

    t_ss = tr(1:SS(ii):end); %time
    u_ss = zf(1:SS(ii):end); %input (wave elevation)
    y_ss = ff(1:SS(ii):end); %output (force)

    %% data packking

    % Xy, Xf and Y cell arrays contains input and output for training a NARX
    % model for each pair of i) time step and ii) segment.

    for jj = 1:1:numSegments

        % sample indices corresponding the the segment jj-th
        index = (1:lenSegment) + (jj-1) * lenSegment ;
        t_temp = t_ss(index);
        f_temp = y_ss(index);
        zeta_temp = u_ss(index);

        Xy{ii,jj} = [] ;
        Xu{ii,jj} = [] ;

        for kk = 1:1:ry
	  %left-most is first .... [t, t-1, t-2, t-3 .., t-ny]
	  Xy{ii,jj} = [Xy{ii,jj}, circshift(f_temp, kk)] ; %+ve k shifts down
        end

        for kk = 0:1:ru
	  Xu{ii,jj} = [Xu{ii,jj},circshift(zeta_temp, kk)] ;
        end

        lagNum = max(ry,ru) ;

        % remove first lagNum points
        Xy{ii,jj} = Xy{ii,jj}(lagNum+1:end,:) ;
        Xu{ii,jj} = Xu{ii,jj}(lagNum+1:end,:) ;
        Y{ii,jj} =  f_temp(lagNum+1:end,:) ;

    end

end

end