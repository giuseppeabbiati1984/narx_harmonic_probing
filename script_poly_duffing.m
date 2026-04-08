% The following code illustrates:
% i) Training of Polynomial-NARX on input/output data of a Duffing oscillator 
% ii) Compute LTF, QTF, and CTF using the numerical harmonic probing algorithm 
% iii) Comparison against the theoretical solutions for the Duffing oscillator

% This implementation is based on:
% Stamenov et al. (2025) "Numerical estimation of generalized frequency response functions from time series data using NARX"
% https://doi.org/10.1016/j.ymssp.2025.113278

close all; clear; clc;
rng(10, 'twister') ;

%% ######## GENERATE NUMERICAL DUFFING DATA  ########################

%------------------------------------- Oscillator input: ----------------------------------------------------------------------------
odePars.m = 1 ; % mass [kg]
odePars.k1 = 1e4 ; % stiffness (linear) [N/m] 
odePars.k2 = 1e7 ; % stiffness (quadratic) [N/m2]
odePars.k3 = 5e9 ; % stiffness (cubic) [N/m3]
odePars.w0 = sqrt(odePars.k1/odePars.m) ; % modal frequency [rad/s]
odePars.f0 = odePars.w0/(2*pi) ; % modal frequency [Hz]
odePars.xi = 0.10 ; % damping ratio [  ]
odePars.c = 2*odePars.xi*sqrt(odePars.k1*odePars.m) ; % damping [Ns/m]

% Time Series Input:
ry = 3 ; % lag samples on output
ru = 3 ; % lag samples on input
ord = "3" ; % order of the poly-NARX
numSeg = 3; % number of segments ( >1 )
lenSeg = 100 ; % number of points in each segment
fsr = 5000 ; % sampling frequency for ODE
dtr = 1/fsr ; % reference time step for ODE [s]
dts = 0.005 ; % subsampling time step(s) for NARX [s]   opt = 5.236e-3
tmax = lenSeg*numSeg*dts;
tr = 0:dtr:tmax ; % high frequency time axis [s]
ts = 0:2*max(dts):tmax ; % sampling time axis [s]
fs = normrnd(0,5,[numel(ts), 1]) ; % input signal [N]
ur = interp1(ts, fs, tr) ; % Nyquist-proof input
clear ts fs


%------------------------------------- Harmonic Probing Input: -----------------------------------------------------------------
diag_offset = 1; % diagonal offset 
dw_res  = 5; % freq-axis discretization
w_min = 30; % start frequency on the axis
w_max = 100; % end frequency on the axis
w1 = w_min:dw_res:floor(w_max/dw_res)*dw_res ; % axis of the LTF, QTF, CTF

% ODE response and force with added noise
SNRf = 100 ; % signal-to-noise ratio for input
SNRy = 100 ; % signal-to-noise ratio for output
odeHandle = @(t, y) compute_ODEduffing(t, y, odePars, tr, ur) ; % state-space definition of the duffing OED
														
[yr, ur] = compute_ODE(odeHandle, tr ,ur, SNRf, SNRy) ; % solving the OED

%Compute scaling factors (input and output)
stdF1 = std(ur) ; % input std. dev. for scaling
stdY1 = std(yr) ; % output std. dev. for scaling

%Compute transfer function scaling factors (nonlinear scaling)
scale0 = stdY1/stdF1^0 ; 
scale1 = stdY1/stdF1^1 ; % LTF Scale
scale2 = stdY1/stdF1^2 ; % QTF Scale
scale3 = stdY1/stdF1^3 ; % CTF Scale

% Scale the signals by the std. dev.
ur_scl = ur/stdF1 ; % Scaled input 
yr_scl = yr/stdY1 ; % Scaled output

% Segment the data into arrays sutiable for training the NARX model
[Xy, Xu, Y, tr_ss, zeta_ss, yr_ss] = pack_data(tr, ur_scl, yr_scl, dts, numSeg, lenSeg, ru, ry) ;

%% Training of Poly-NARX
NARX = cell(size(Xy));

    for i = 1:numel(Xy)
        NARX{i} = compute_poly_NARX(Xy{i}, Xu{i}, Y{i}, ord) ; %Train a NARX model for each of the segments. NARX one-step ahead is loaded based on the input data size
    end
%% Compute/Load theoretical transfer functions
freq_max = 200 ; %max limit of the frequency axis

%Uncomment to compute theoretical values at a particular frequency discretization
% freq_disc = 4; %discretization of the frequency axis [rad/s]
% w_theo = [0:freq_disc:freq_max] ; % Discretization of the frequncy axis
% [H1_course, H2_course, H3_course, w_dbl_course] = fun_analytical_duffing(odePars, w_theo) ;
% save DuffTF_dw=4.mat H1_course H2_course H3_course w_dbl_course

load DuffTF_dw=4.mat % Loads the precomputed exact transfer functions with resolution dw = 2 rad/s

% Interpolate to a finer grid (dw=1 rad/s) for plotting:
[X_int, Y_int, Z_int] = meshgrid(-freq_max:1:freq_max);
H1_exact = interp1(w_dbl_course, H1_course, X_int(1,:,1), "spline");
H2_exact = interp2(w_dbl_course, w_dbl_course, H2_course, X_int(:,:,1), Y_int(:,:,1), "spline");
H3_exact = interp3(w_dbl_course, w_dbl_course, w_dbl_course,H3_course, X_int,  Y_int, Z_int, "spline");
w_dbl_exact = -freq_max:1:freq_max; %frequency axis used for comparison in the plots

%% ############### Probing #########################

% Allocating memory for arrays:
H1_scl = zeros(numSeg, numel(w1)) ; %vector for the SCALED output of the probing
H1_sol = zeros(numSeg, numel(w1)) ;  %vector for the UNSCALED solution of the probing
H2_mat_scl = zeros(numel(w1), numel(w1),numSeg) ;
H2_mat_sol = zeros(numel(w1), numel(w1),numSeg) ;
H3_mat_scl = zeros(numel(w1), numel(w1),numel(w1),numSeg) ;
H3_mat_sol = zeros(numel(w1), numel(w1),numel(w1),numSeg) ;

for i = 1:1:size(Y, 2)
disp(strcat('HP of segment : ', num2str(i), '/', num2str(size(Y, 2)))) ;

    C1 = NARX{i}.C1(:).' ; % First order array of coefficients (Jacobian). [1 x r]
    C2 = NARX{i}.C2 ; % Second order array of coefficients (Hessian). [r x r]

    % Numerical Probing
    % ####### H1 (LTF) ################################
    H1_scl(i, :) = compute_probing_H1(C1, ry, ru, w1, dts) ; % Numerical probing of LTF
    H1_sol(i, :) = scale1 .* H1_scl(i, :) ; % Unscale the LTF


    % ####### H2 (QTF) ################################
    H2_mat_scl(:, :, i) = compute_probing_H2(C1, C2, ry, ru , H1_scl(i, :), w1, dts) ; % Numerical probing of QTF
    H2_mat_sol(:, :, i) = scale2 .* H2_mat_scl(:, :, i); % Unscale the QTF


    % ####### H3 (CTF) ################################
    if strcmp(ord, "3")
        C3 = NARX{i}.C3 ; % Tensor of NARX coefficients. [r x r x r]
        H3_mat_scl(:,:,:,i) = compute_probing_H3(C1, C2, C3, ry, ru , H1_scl(i, :), H2_mat_scl(:, :, i), w1, dts) ; % % Numerical probing of CTF
        H3_mat_sol(:,:,:,i) = scale3 .* H3_mat_scl(:,:,:,i) ; % Unscale the CTF
    end

end

%% Plot LTF
close all ;
PlotFontSize = 27 ;
NumStd = 3 ; % How many std dev to plot
Color1 = [.55 .55 .55]; Color2 = [.25 .95 .95]; Color3 = [.05 .05 .05]; color4 = [.25 .25 .25] ; % Define colors


%Obtain LTF mean and std  values
H1_mean = mean(H1_sol, 1) ;
H1_std_real  = std(real(H1_sol),0,1) ; H1_std_imag  = std(imag(H1_sol),0,1) ; H1_std = H1_std_real + 1i*H1_std_imag ;
y_max = max([imag(H1_exact),real(H1_exact)], [], 'all') ;
y_min  = min([imag(H1_exact),real(H1_exact)], [], 'all') ;

H1_2std = (H1_mean+ 2*[H1_std;  -H1_std]) ;
H1_3std = (H1_mean+ 3*[H1_std;  -H1_std]) ;
y_plot2=[H1_2std(1, :), fliplr(H1_2std(2, :))] ;
y_plot3=[H1_2std(1, :), fliplr(H1_3std(2, :))] ;
x_plot =[w1, fliplr(w1)];


% ########PLOT LTF REAL PART ##################################
    figure('Renderer', 'painters', 'Position', [10, 100, 800, 900]) ;
     subplot(2,1,1);  hold on;  box on; set(gca, 'FontSize',  PlotFontSize) ;

     f3 = fill(x_plot, real(y_plot3), 1,'facecolor', Color3, 'edgecolor', 'none', 'facealpha', 0.2);
     f2 = fill(x_plot, real(y_plot2), 1,'facecolor', Color2, 'edgecolor', 'none', 'facealpha', 0.3);

     for i = 1:size(H1_sol,1)
         h1 = plot(w1, real(H1_sol(i, :)),  'kx',  'Color', color4, 'HandleVisibility','off');
     end

     h2 = plot(w_dbl_exact, real(H1_exact),'r-','linewidth', 1.5) ;  % first column is the diagonal
     hm = plot(w1, real(H1_mean), 'bx','linewidth', 1.5) ;
     ylim(1.2*[y_min, y_max]); xlim([0, w_dbl_exact(end)]);

     title({'Linear'},'Interpreter','latex',"FontSize", PlotFontSize) 
     xlabel('$\omega_1$ [rad/s]','Interpreter','latex')
     ylabel('$\Re(H^{(1)})$','Interpreter','latex')

     legend([h2, hm, f2, f3], { "Theoretical", 'Probing Mean' , "2 std. dev", "3 std. dev"}, 'location', 'northeast', 'fontsize', PlotFontSize-12) ;
     set(gca, 'FontSize',  PlotFontSize) ;

     % ########PLOT LTF IMAGINARY PART ##################################
     subplot(2,1,2); hold on; box on; set(gca,'FontSize',  PlotFontSize)

     f3 = fill(x_plot, imag(y_plot3), 1,'facecolor', Color3, 'edgecolor', 'none', 'facealpha', 0.2);
     f2 = fill(x_plot, imag(y_plot2), 1,'facecolor', Color2, 'edgecolor', 'none', 'facealpha', 0.3) ; % fill plot of the std. dev.

     for i = 1:size(H1_sol,1)
          h1 = plot(w1, imag(H1_sol(i, :)),  'kx',  'Color', color4, 'HandleVisibility','off'); 
     end

     h2 = plot(w_dbl_exact, imag(H1_exact), 'r-', 'linewidth', 1.5) ; % exact solution 
     hm = plot(w1, imag(H1_mean), 'bx','linewidth', 1.5) ; % mean of all samples
     ylim(1.2*[y_min, y_max]); xlim([0, w_dbl_exact(end)]);
     xlabel('$\omega_1$ [rad/s]','Interpreter','latex')
     ylabel('$\Im(H^{(1)})$','Interpreter','latex')
     set(gca, 'FontSize',  PlotFontSize) ;

%% Plot Diagonal of QTF

H2_sol = zeros(size(Y, 2),numel(w1)-diag_offset) ;
H2_exact_diag = zeros(1,numel(w_dbl_exact)) ;


for i = 1:1:size(Y, 2)
	H2_sol(i, :) = diag(H2_mat_sol(:,:,i), diag_offset) ;
end

diag_offset_exact = diag_offset*(dw_res/(w_dbl_exact(2)-w_dbl_exact(1))) ;
H2_exact_diag(1:end-diag_offset_exact)  = diag(H2_exact, diag_offset_exact) ;


y_max = max([imag(H2_exact_diag),real(H2_exact_diag)], [], 'all');
y_min  = min([imag(H2_exact_diag),real(H2_exact_diag)], [], 'all');           
x_plot =[w1(1:end-diag_offset), fliplr(w1(1:end-diag_offset))];

H2_mean = mean(H2_sol, 1) ;
H2_std_real  = std(real(H2_sol),0,1) ; H2_std_imag  = std(imag(H2_sol),0,1) ; 
H2_std = H2_std_real + 1i*H2_std_imag ;


H2_2std = (H2_mean + 2*[H2_std;  -H2_std]) ;
H2_3std = (H2_mean + 3*[H2_std;  -H2_std]) ;
y_plot2=[H2_2std(1, :), fliplr(H2_2std(2, :))] ;
y_plot3=[H2_3std(1, :), fliplr(H2_3std(2, :))] ;

% ########PLOT QTF REAL PART ##################################
    figure('Renderer', 'painters', 'Position', [600, 100, 800, 900]);
     subplot(2,1,1);  hold on;  box on; set(gca, 'FontSize',  PlotFontSize);

     f3 = fill(x_plot, real(y_plot3), 1,'facecolor', Color3, 'edgecolor', 'none', 'facealpha', 0.2);
     f2 = fill(x_plot, real(y_plot2), 1,'facecolor', Color2, 'edgecolor', 'none', 'facealpha', 0.3);

     for j = 1:size(H2_sol,1)
         h1 = plot(w1(1:end-diag_offset), real(H2_sol(j, :)),  'kx', 'Color', color4, 'HandleVisibility','off'); 
     end

     h2 = plot(w_dbl_exact, real(H2_exact_diag), 'r-', 'linewidth', 1.5) ; % exact solution
     hm = plot(w1(1:end-diag_offset), real(H2_mean), 'bx','linewidth', 1.5) ; % mean of all samples

     ylim(1.2*[y_min, y_max]); xlim([0, w_dbl_exact(end)]);
     % title({[strcat(' Freq. $\Delta \omega =$', num2str(dw), ' rad/s, ')]; ...
     % [strcat(' $\omega_n =$', num2str(round(wn)), ' rad/s')]}, 'Interpreter','latex')
     title({["Quadratic"]; ...
          [strcat("$ |\omega_2 - \omega_1| =", num2str(dw_res*diag_offset) ,"~r/s$")]}, 'Interpreter','latex',"FontSize", PlotFontSize) 
     xlabel('$\omega_1$ [rad/s]','Interpreter','latex')
     ylabel('$\Re(H^{(2)})$','Interpreter','latex')
     legend([h2, hm, f2, f3], { "Theoretical", 'Probing Mean', "2 std. dev", "3 std. dev"}, 'location', 'northeast', 'fontsize', PlotFontSize-12) ;
     set(gca, 'FontSize',  PlotFontSize) ;

     % ########PLOT QTF IMAGINARY PART ##################################
     subplot(2,1,2); hold on; box on; set(gca,'FontSize',  PlotFontSize)

     f3 = fill(x_plot, imag(y_plot3), 1,'facecolor', Color3, 'edgecolor', 'none', 'facealpha', 0.2);
     f2 = fill(x_plot, imag(y_plot2), 1,'facecolor', Color2, 'edgecolor', 'none', 'facealpha', 0.3);

     for j = 1:size(H2_sol,1)
        h1 =  plot(w1(1:end-diag_offset), imag(H2_sol(j, :)), 'kx', 'Color', color4, 'HandleVisibility','off'); 
     end

     h2 = plot(w_dbl_exact, imag(H2_exact_diag), 'r-', 'linewidth', 1.5) ; % exact solution 
     hm = plot(w1(1:end-diag_offset), imag(H2_mean), 'bx','linewidth', 1.5) ; % mean of all samples

     ylim(1.2*[y_min, y_max]); xlim([0, w_dbl_exact(end)]);
     xlabel('$\omega_1$ [rad/s]','Interpreter','latex')
     ylabel('$\Im(H^{(2)})$','Interpreter','latex')
     set(gca, 'FontSize',  PlotFontSize) ;
     % legend([h1, h1_qtf, f_qtf], {'Harmonic Probing', "Theoretical", strcat(num2str(NumStd), ' std. dev') }, 'location', 'northeast',"FontSize",PlotFontSize-7) ;
     
%% Plot Diagonal of CTF
ax1 = w1(2):dw_res:w1(end)-3*dw_res*diag_offset ; % comparison axis 1
ax2 = ax1+diag_offset*dw_res ; % comparison axis 2
ax3 = ax2+diag_offset*dw_res  ; % comparison axis 3


H3_sol = zeros(1,numel(ax1)) ;
for i = 1:1:size(Y, 2)
    for j = 1:numel(ax1)
        [~, ind1] = min(abs(ax1(j)-w1)) ;
        [~, ind2] = min(abs(ax2(j)-w1)) ;
        [~, ind3] = min(abs(ax3(j)-w1)) ;

        H3_sol(i, j) = H3_mat_sol(ind1, ind2, ind3, i) ;
    end
end


H3_exact_diag = zeros(1, numel(ax1)) ;
ax1_exact = w_dbl_exact ;
ax2_exact = ax1_exact+diag_offset*dw_res/(w_dbl_exact(2)-w_dbl_exact(1)) ;
ax3_exact = ax2_exact+diag_offset*dw_res/(w_dbl_exact(2)-w_dbl_exact(1))  ;

for j = 1:numel(ax1_exact)
    [~, ind1] = min(abs(ax1_exact(j)-w_dbl_exact)) ;
    [~, ind2] = min(abs(ax2_exact(j)-w_dbl_exact)) ;
    [~, ind3] = min(abs(ax3_exact(j)-w_dbl_exact)) ;

    H3_exact_diag(j) =H3_exact(ind1, ind2, ind3) ;
end


y_max = max([imag(H3_exact_diag), real(H3_exact_diag)], [], 'all');
y_min =  min([imag(H3_exact_diag), real(H3_exact_diag)], [], 'all');   
x_plot =[ax1, fliplr(ax1)];

H3_mean = mean(H3_sol, 1) ;
H3_std_real  = std(real(H3_sol),0,1,"omitmissing") ; H3_std_imag  = std(imag(H3_sol),0,1,"omitmissing") ; 
H3_std = H3_std_real + 1i*H3_std_imag ;

H3_2std = (H3_mean + 2*[H3_std;  -H3_std]) ;
H3_3std = (H3_mean + 3*[H3_std;  -H3_std]) ;
y_plot2=[H3_2std(1, :), fliplr(H3_2std(2, :))] ;
y_plot3=[H3_3std(1, :), fliplr(H3_3std(2, :))] ;

% ########PLOT CTF REAL PART ##################################
    figure('Renderer', 'painters', 'Position', [1200, 100, 800, 900]);
     subplot(2,1,1);  hold on;  box on; set(gca, 'FontSize',  PlotFontSize);

     f3 = fill(x_plot, real(y_plot3), 1,'facecolor', Color3, 'edgecolor', 'none', 'facealpha', 0.2);
     f2 = fill(x_plot, real(y_plot2), 1,'facecolor', Color2, 'edgecolor', 'none', 'facealpha', 0.3);

     for j = 1:size(H3_sol,1)
        h1 = plot(ax1, real(H3_sol(j, :)),  'kx', 'Color', color4, 'HandleVisibility','off'); 
     end

     h2 = plot(w_dbl_exact, real(H3_exact_diag) , 'r-', 'linewidth', 1.5) ; % exact solution 
     hm = plot(ax1, real(H3_mean), 'bx', 'linewidth', 1.5) ; % mean of all samples

     ylim(1.2*[y_min, y_max]); xlim([0, w_dbl_exact(end)]);

     title({['Cubic'] ...
         [strcat("$ |\omega_2 - \omega_1| = $", num2str(diag_offset*dw_res), " rad/s ")]...
         [strcat("$ |\omega_3 - \omega_2| = $", num2str(diag_offset*dw_res), " rad/s ")]}, 'Interpreter','latex',"FontSize", PlotFontSize) 
     xlabel('$\omega_1$ [rad/s]','Interpreter','latex')
     ylabel('$\Re(H^{(3)})$','Interpreter','latex')
    legend([h2, hm, f2, f3], { "Theoretical", 'Probing Mean',  "2 std. dev", "3 std. dev"}, 'location', 'northeast', 'fontsize', PlotFontSize-12) ;
     set(gca, 'FontSize',  PlotFontSize) ;

  % ########PLOT CTF IMAGINARY PART ##################################
     subplot(2,1,2); hold on; box on; set(gca,'FontSize',  PlotFontSize)

     f3 = fill(x_plot, imag(y_plot3), 1,'facecolor', Color3, 'edgecolor', 'none', 'facealpha', 0.2);
     f2 = fill(x_plot, imag(y_plot2), 1,'facecolor', Color2, 'edgecolor', 'none', 'facealpha', 0.3);

     for j = 1:size(H2_sol,1)
         h1 = plot(ax1, imag(H3_sol(j, :)),  'kx', 'Color', color4, 'HandleVisibility','off');
     end

     h2 = plot(w_dbl_exact, imag(H3_exact_diag), 'r-', 'linewidth', 1.5) ; % exact solution 
     hm = plot(ax1, imag(H3_mean), 'bx','linewidth', 1.5) ; % mean of all samples

     ylim(1.2*[y_min, y_max]); xlim([0, w_dbl_exact(end)]);
     xlabel('$\omega_1$ [rad/s]','Interpreter','latex')
     ylabel('$\Im(H^{(3)})$','Interpreter','latex')
     set(gca, 'FontSize',  PlotFontSize) ;