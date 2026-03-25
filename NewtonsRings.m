clear all;clc; format long;

%{
---------------------------------------------------------------------------
This program is designed to import, store, convert, calculate and display
data recorded thoughout an experiment where we measured the distance
between the consecutive dark fringes of an interference pattern for a
convex lens, producing Newton's Rings, and two flat lenses with a thin
wedge in between, producing a thin film interference pattern.
---------------------------------------------------------------------------
%}



% --- Thin Film Interference Analysis ---

data_TW = readtable('data\TW_dataREF.csv', 'NumHeaderLines', 2);               % Imports data from Thin Wedge observations

Len = 5.67*10^(-2);                         % Length of glass slide, measured directly
lambda = 589*10^(-9);                       % Wavelength of sodium lamp source

TW_h = cell(3, 1);                          % Stores data from multiple trials (t)
TW_d = [3, 1]';                             % Stores the distance (d) between the inital position on the slide and the final position

for t = 1:3

    TW_h{t} = table2array(data_TW(t, 2:3));                 
    TW_d(t) = TW_h{t}(2) - TW_h{t}(1);

end


TW_d_avg = 1/3 * (TW_d(1) + TW_d(2) + TW_d(3));
delta_x = TW_d_avg / 40;                                        % Average distance between each dark fringe (should be constant)

h = (lambda * Len) / (2 * delta_x);                               % Thickness of thin wedge



    % Uncertainty anaylsis
sigma_Ver = ((1/2 * 1/50*10^(-3))) * ones(12, 1);               % Uncertainty in physical measurements taken by travelling miscroscope

sigma_TW_d_stat = std(TW_d);
sigma_TW_d_mean = sigma_TW_d_stat / sqrt(3);

sigma_TW_d_instr = sqrt(2) * sigma_Ver(1);
sigma_TW_d_total = sqrt(sigma_TW_d_mean^2 + sigma_TW_d_instr^2);

sigma_delta_x = sigma_TW_d_total / 40;

sigma_L = 0.01 * 10^(-3);
sigma_h = h * sqrt((sigma_L / Len)^2 + (sigma_delta_x / delta_x)^2);




% --- Newton's Rings Analysis ---

data_NR = readtable('data\NR_dataREF.csv', 'NumHeaderLines', 3);       % Imports data from Newton's Rings observations

m = table2array(data_NR(1:12, 1));          % m^th dark fringe

NR_L = cell(2, 1);

NR_L{1} = table2array(data_NR(1:12, 2));            % Newton's Rings: Left measuremet trial 1
NR_L{2} = table2array(data_NR(1:12, 4));


NR_R = cell(2, 1);

NR_R{1} = table2array(data_NR(1:12, 7));            % Newton's Rings: Right measuremet trial 1
NR_R{2} = table2array(data_NR(1:12, 9));


r_L = cell(2, 1);
r_L_avg = zeros(12, 1);

r_R = cell(2, 1);
r_R_avg = zeros(12, 1);

r = zeros(12, 1);                                   % Average r from left and right measurements

for i = 1:12

    r_L{1}(i, 1) = abs(NR_L{1}(i) - NR_L{1}(1));
    r_L{2}(i, 1) = abs(NR_L{2}(i) - NR_L{2}(1));
    r_L_avg(i) = 1/2 * (r_L{1}(i) + r_L{2}(i));
    
    r_R{1}(i, 1) = abs(NR_R{1}(i) - NR_R{1}(1));
    r_R{2}(i, 1) = abs(NR_R{2}(i) - NR_R{2}(1));
    r_R_avg(i) = 1/2 * (r_R{1}(i) + r_R{2}(i));

    r(i) = 1/2 * (r_L_avg(i) + r_R_avg(i));

end


    % Uncertainty analysis
sigma_r_L = sqrt((sigma_Ver).^2 + (sigma_Ver).^2);
sigma_r_R = sigma_r_L;

r_L_mat = [r_L{1}, r_L{2}];   
r_R_mat = [r_R{1}, r_R{2}];

sigma_r_L_stat = std(r_L_mat, 0, 2);   % sample standard deviation across columns
sigma_r_R_stat = std(r_R_mat, 0, 2);

sigma_r_L_mean = sigma_r_L_stat / sqrt(2);
sigma_r_R_mean = sigma_r_R_stat / sqrt(2);

sigma_r_L_total = sqrt(sigma_r_L_mean.^2 + sigma_r_L.^2);
sigma_r_R_total = sqrt(sigma_r_R_mean.^2 + sigma_r_R.^2);

sigma_r = 0.5 * sqrt(sigma_r_L_total.^2 + sigma_r_R_total.^2);


r = r(2: end);
m = m(2: end);
sigma_r = sigma_r(2: end);
sigma_r2 = 2 * r .* sigma_r;


    % Plot of m^th dark fringe against radius squared to give linear relationship
figure(1)
hold on
plot(m, r.^2, 'k*')
errorbar(m , r.^2, sigma_r2, 'ko')

p = polyfit(m, r.^2, 1);
m_fit = polyval(p, m);

plot(m, m_fit, 'b')

xlabel('m', fontsize = 18)
ylabel('r^2 (m^2)', fontsize = 18)
title('Radial Distance Squared vs m^{th} Dark Fringe')

set(gca, fontsize = 16)
grid on
exportgraphics(gcf, 'R2vsM.png', 'ContentType', 'vector')


    % Radius of curvature for convex lens
Rad = p(1) / lambda;        


    %Residual Analysis for Newton's Rings
residuals = r.^2 - m_fit;
chi2 = sum((residuals ./ sigma_r2) .^2);
nu = length(r) - length(p);
chi2_reduced = chi2 / nu;

    figure(2)
    plot(m_fit, residuals, 'ko', 'MarkerFaceColor', 'k')
    yline(0, 'k--')

    xlabel('r^2 (m^2)', fontsize = 18)
    ylabel('Residuals', fontsize = 18)
    title('Residuals vs Radial Distance Squared')

    set(gca, fontsize = 16)
    grid on
    axis padded
    exportgraphics(gcf, 'ResidualsVsR2.png', 'ContentType', 'vector')




% --- Extra plots which weren't included; diagnotics ---
    
    % Plot of radial distance vs m^th dark fringe
figure(3)
hold on
plot(m, r, 'k*')
errorbar(m , r, sigma_r, 'ko')

p_quad = polyfit(m, r, 2);
m_fit_quad = polyval(p_quad, m);

p_lin = polyfit(m, r, 1);
m_fit_lin = polyval(p_lin, m);


quad_fit = plot(m, m_fit_quad, 'b', 'DisplayName', 'Quadratic Regression');         % quadratic fit of data
lin_fit = plot(m, m_fit_lin, 'r', 'DisplayName', 'Linear Regression');              % linear fit of data

xlabel('m', fontsize = 18)
ylabel('r (m)', fontsize = 18)
title('Radial Distance vs m^{th} Dark Fringe')

legend([quad_fit, lin_fit], 'Location', 'east')

set(gca, fontsize = 16)
grid on
exportgraphics(gcf, 'RvsM.png', 'ContentType', 'vector')

residuals_quad = r - m_fit_quad;
chi2_quad = sum((residuals_quad ./ sigma_r) .^2);
nu_quad = length(r) - length(p_quad);
chi2_reduced_quad = chi2_quad / nu_quad;

    figure(4)
    plot(m_fit_quad, residuals_quad, 'ko', 'MarkerFaceColor', 'k')
    yline(0, 'k--')

    xlabel('r (m)', fontsize = 18)
    ylabel('Residuals', fontsize = 18)
    title('Residuals vs Radial Distance for Quadratic Regression')

    set(gca, fontsize = 16)
    grid on
    axis padded
    exportgraphics(gcf, 'RvsMQR.png', 'ContentType', 'vector')


residuals_lin = r - m_fit_lin;
chi2_lin = sum((residuals_lin ./ sigma_r) .^2);
nu_lin = length(r) - length(p_lin);
chi2_reduced_lin = chi2_lin / nu_lin;

    figure(5)
    plot(m_fit_lin, residuals_lin, 'ko', 'MarkerFaceColor', 'k')
    yline(0, 'k--')

    xlabel('r (m)', fontsize = 18)
    ylabel('Residuals', fontsize = 18)
    title('Residuals vs Radial Distance for Linear Regression')

    set(gca, fontsize = 16)
    grid on
    axis padded
    exportgraphics(gcf, 'RvsMLR.png', 'ContentType', 'vector')

    %{
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    %}