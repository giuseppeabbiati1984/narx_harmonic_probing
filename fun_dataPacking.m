function [Xf, Xz, F, tr_ss, zf_ss, ff_ss] = fun_dataPacking(tr, zr, fr, dts, numSegments, lenSegment, lagNumZ, lagNumF)


        %ensure column vectors
        tr = reshape(tr,[length(tr), 1]);
        zr = reshape(zr,[length(zr), 1]);
        fr = reshape(fr,[length(fr), 1]);

        dt = tr(2)-tr(1);

        SS = floor(dts/dt); %integer subsampling rate
        % SS = dts/dt; %integer subsampling rate
        
        if max(SS)*lenSegment*numSegments > numel(tr)
            error('Number of segments or length of segment too large')
        end


        % Construct a causal filter
        [b, a] = butter(3, 1/SS, 'low') ;
        ff = filter(b, a, fr) ; % filtered force signal (time domain)
        zf = filter(b, a, zr) ; % filtered wave elevation signal (time domain)

for ii = 1:numel(SS)

        tr_ss = tr(1:SS(ii):end); %time
        zf_ss = zf(1:SS(ii):end); %input (wave elevation)
        ff_ss = ff(1:SS(ii):end); %output (force)

    %% data packking

    % Xy, Xf and Y cell arrays contains input and output for training a NARX
    % model for each pair of i) time step and ii) segment.
    
    for jj = 1:1:numSegments

        % sample indices corresponding the the segment jj-th
        index = (1:lenSegment) + (jj-1) * lenSegment ;
        t_temp = tr_ss(index);
        f_temp = ff_ss(index);
        zeta_temp = zf_ss(index);

        Xf{ii,jj} = [] ;
        Xz{ii,jj} = [] ;
        
        for kk = 1:1:lagNumF
            %left-most is first .... [t, t-1, t-2, t-3 .., t-ny]
            Xf{ii,jj} = [Xf{ii,jj}, circshift(f_temp, kk)] ; %+ve k shifts down
        end
        
        for kk = 0:1:lagNumZ
            Xz{ii,jj} = [Xz{ii,jj},circshift(zeta_temp, kk)] ;
        end
        
        lagNum = max(lagNumF,lagNumZ) ;
        
        % remove first lagNum points
        Xf{ii,jj} = Xf{ii,jj}(lagNum+1:end,:) ;
        Xz{ii,jj} = Xz{ii,jj}(lagNum+1:end,:) ;
        F{ii,jj} =  f_temp(lagNum+1:end,:) ;
           
    end
    
end

end