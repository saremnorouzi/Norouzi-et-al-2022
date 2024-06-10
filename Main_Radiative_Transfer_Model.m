% Code to calibrate the new proposed radiative transfer model. Please feel
% free to contact me via email: sarem.nrz@gmail.com

clc;
clear;
close all


% -------------------------------------------------------------------------
% Please enter the proper values for the following flags and run the code
% -------------------------------------------------------------------------
flag_opt = 1; % Put 1 if a new optimization is needed. For values other than 1 the code loads the lastest saved results (if exist)
% For second run or if you do not need to repeat the optimization, put 0
linear_model_flag = 0; % When linear_model_flag = 1, the code uses a linea mixing model,i.e., p_a = p_c = 1, otherwise it uses the general nonlinear model

soilname = 'AZ12'; 
Bands = [1650, 2210]; % Please specify the wavelengths at which model predictions are plotted.
sug_Sw = [0.1, 0.3, 0.6, 0.8].'; % Please specify degree of saturations (theta/theta_sat) used in the last plot of predicted vs. measured reflectance


% -------------------------------------------------------------------------
% Importing the measured reflectance-moisture data
% -------------------------------------------------------------------------
% Here we load the .mat file containing hyperspectral measured reflectance
% during evaporation experiment. The first row of the file contains the
% water contents at which reflectance values have been recorded. Each of
% the following rows contains the measured reflectance values corresponding
% to a wavelength. The first column contains the wavelengths.

load([soilname '_Moisture_Reflectance_Data.mat'], 'theta_sat','R_Mois_All') % This file contains the measured hyperspectral reflectance data during evaporation experiment
m = size(R_Mois_All,1); % Number of wavelengths + 1
n = size(R_Mois_All,2); % Number of theta values + 1
theta_meas = R_Mois_All(1,2:n); % Theta vector
Wave = R_Mois_All(2:m,1); % Wavelengths
refl = R_Mois_All(2:m,2:n); % All the measured reflectances
rtran = (1-refl).^2./(2.*refl); % Converting all reflectance values to transformed reflectance [See Eq. (4) of our paper]



% -------------------------------------------------------------------------
% Partitioning the water content values of evaporation experiment into
% capillary and adsorbed water
% -------------------------------------------------------------------------
% Here and prior to optimizing the proposed radiative transfer model for
% finding model parameters (i.e., c_a, c_c, p_a and p_c) we need to find
% the capillary and adsorbed water components at water contents
% corresponding to evaporation experiment. To that end, we load the Lebeau
% and Konrad (2010) parameters that we optimized in previous code. Then we
% interpolate over the capillary and adsorbed curve to find the water
% Components at the desired intervals.

version = 'Lebeau_Konrad'; % Please do not change this
% Loading the retention curve parameters
load(['All_AZ_Soils_' soilname '_' version '.mat'],'x')

% Solution
theta_0 = x(1);
h_median = x(2);
sigma = x(3);

% Generating a very refined curve of SWRC and water components
a = -2; b = 5; % Range of matric potential: 10^a to 10^b
head = -logspace(a,b,1000); % A very refined and logarithmically distributed vector of matric potential
h_dry = -10^5; % Matric potential at oven dryness

theta_c = theta_sat.*(0.5.*erfc( log(head./h_median)./(sqrt(2).*sigma) )); % Capillary component
theta_a = theta_0.*(1 - log(abs(head))./log(abs(h_dry)) ) .* (1 - theta_c./theta_sat); % Adsorptive component
theta_p = theta_c + theta_a; % Total water content


% Finding the capillary and adsorptive components at water contents of
% reflectance experiment by interpolating over the refined curves.
theta_cap = interp1(theta_p,theta_c,theta_meas,'spline'); % Capillary water at measured water contents of evaporation experiment
theta_ads = interp1(theta_p,theta_a,theta_meas,'spline'); % Adsorptive water at measured water contents of evaporation experiment



% Plot of water components at measured intervals
figure('name','Capillary versus Water content')

plot(theta_p,theta_c,'b-','Linewidth',2.0) % Theoretical Capillary Water
hold on
plot(theta_meas, theta_cap, 'bo') % Capillary Water at Measured Intervals
plot(theta_p,theta_a,'r-','Linewidth',2.0) % Theoretical Adsorbed Water
plot(theta_meas, theta_ads, 'ro') % Adsorbed Water at Measured Intervals

xlabel('Volumetric Water Content [m^{3} m^{-3}]','FontSize',14,'Color','k')
ylabel('Capillary and Adsorbed Water [m^{3} m^{-3}]','FontSize',14)
legend({'Capillary Water', 'Capillary Water at Measured Intervals' ...
    , 'Adsorbed Water', 'Adsorbed Water at Measured Intervals'},'fontsize',12,'Location','northwest')
axis([0 0.45 0 0.45])
title(soilname)
pbaspect([1.3 1 1])



% -------------------------------------------------------------------------
% Optimization algorithm for finding optical coefficients
% -------------------------------------------------------------------------
% Here we find the optimum c_a, c_c, p_a and p_c via least squares fit of
% Eq(8) to the measurements and repeat this process at each wavelength to
% obtain the SWIR spectrum of these coefficients.

result = zeros(length(Wave), 6); % Empty matrix to store the results
x_old = [10 10 1 1]; % Initial guess for optical coefficients

% Nonlinear mixing model, i.e., p_a and p_c could take any values

x_lower = [eps eps eps eps]; % Lower bounds corresponding to [c_a, c_c, p_a, p_c]
x_upper = [500 500 2 2]; % Upper bounds corresponding to [c_a, c_c, p_a, p_c]


% linear mixing model, i.e., p_a = p_c = 1
if linear_model_flag == 1
    
    x_lower = [eps eps 1 1]; % Lower bounds corresponding to [c_a, c_c, p_a, p_c]
    x_upper = [500 500 1 1]; % Upper bounds corresponding to [c_a, c_c, p_a, p_c]
    
end


wavenum_start = 1200; % Starting wavelength
i_start = wavenum_start - 349;

wavenum_end = 2500; % End wavelength
i_end = wavenum_end - 349;

SWIR_range = (wavenum_start <= Wave) & (Wave <= wavenum_end); % Bolean vector determining the SWIR range


% Optimization process
if flag_opt == 1
    
    for i = i_start:i_end
        
        rng default
        
        rtran_meas = rtran(i,:); % Measured transformed reflectance
        r_d = rtran(i, end); % Dry transformed reflectacne
        
        fun = @(x)Fit_General_Form(x, r_d, rtran_meas, theta_ads, theta_cap); % Function to be minimized
        
        gs = GlobalSearch('FunctionTolerance',1e-20,'NumTrialPoints',500,'XTolerance',1e-20,'Display','off');
        problem = createOptimProblem('fmincon','x0',x_old,...
            'objective',fun,'lb',x_lower,'ub',x_upper);
        ms = MultiStart(gs);
        x = run(ms,problem,1);
        x_old = x;
        
        % Optimized optical properties
        c_a = x(1);
        c_c = x(2);
        p_a = x(3);
        p_c = x(4);
        
        
        result(i,1) = Wave(i);
        result(i,2) = r_d;
        result(i,3) = c_a;
        result(i,4) = c_c;
        result(i,5) = p_a;
        result(i,6) = p_c;
        
        disp(Wave(i))
        
    end
    
    
    save([soilname '_FRW_Model_Fitted_Spectrum_Coeffs'  '.mat'], 'result')
else
    load([soilname '_FRW_Model_Fitted_Spectrum_Coeffs'  '.mat'], 'result')
end

% -------------------------------------------------------------------------
% Plot of optical coefficients vs. Wavelength
% -------------------------------------------------------------------------
figure('name','c_a and c_c');
plot(result(SWIR_range,1), result(SWIR_range,3), 'r-','Linewidth',1.5)
hold on;
plot(result(SWIR_range,1), result(SWIR_range,4), 'b-','Linewidth',1.5)
xlabel('Wavelength [nm]','FontSize',14,'Color','k')
ylabel('Coefficients of the Radiative Transfer Equation, c [-]','FontSize',14)
legend({'c_a Spectrum (Adsorbed Water)','c_c Spectrum (Capillary Water)'},'fontsize',13,'location','northwest')
axis([1100 2600 0 25])  %
title(soilname)
xticks([1300 1700 2100 2500])
yticks([0 5 10 15 20 25])
legend boxoff


figure('name','p_a and p_c');
plot(result(SWIR_range,1), result(SWIR_range,5), 'r-','Linewidth',1.5);
hold on;
plot(result(SWIR_range,1), result(SWIR_range,6), 'b-','Linewidth',1.5)
xlabel('Wavelength [nm]','FontSize',14,'Color','k')
ylabel('Coefficients of the Radiative Transfer Equation, p [-]','FontSize',14,'Color','k')
legend({'p_a Spectrum (Adsorbed Water)','p_c Spectrum (Capillary Water)'},'fontsize',13,'location','northwest')
axis([1100 2600 0 3])  %
title(soilname)
xticks([1300 1700 2100 2500])
yticks([0 0.5 1.0 1.5 2.0 2.5])
legend boxoff


% -------------------------------------------------------------------------
% Plot of predicted versus measured transformed reflectance at selected
% wavelengths
% -------------------------------------------------------------------------
l = length(Bands);

figure('name','Fitted versus Measured')
for j = 1:l
    
    wavenum = Bands(j);
    index = wavenum - 349;
    
    % Wave(i) = result(i,1);
    r_d = result(index,2);
    c_a = result(index,3);
    c_c = result(index,4);
    p_a = result(index,5);
    p_c = result(index,6);
    
    
    r_pre = r_d + (c_a .* theta_a .^ p_a) + (c_c .* theta_c .^ p_c); % Predicted transformed reflectance [Eq. (8)]
    
    % The final equation as text
    txt = ['\it r \rm = ' num2str(r_d,2) ' + ' num2str(c_a,3) ' \theta_{\it a}'...
        ' + ' num2str(c_c,3) ' \theta_{\it c}']; 
    
    r_pre_ad = real(c_a .* theta_a .^ p_a); % Second Term of Eq. (8)
    r_pre_ca = real(c_c.*theta_c.^p_c); % Third Term of Eq. (8)
    r_pre = real(r_pre);
    
    max_r_pre_band(j) = max(r_pre); %#ok<SAGROW>
    
    % Plot of transformed reflecatance at selected bands
    subplot(2, l, j)
    plot(theta_meas,rtran(index,:),'ro','Markersize',7,'Linewidth',1.5)
    hold on
    plot(theta_p,r_pre,'k-','Linewidth',2.0)
    plot(theta_p,r_pre_ad,'r--','Linewidth',2.0) % Adsorptive reflectance component [second term on the right hand side of Eq. (8)]
    plot(theta_p,r_pre_ca,'b-','Linewidth',2.0) % Capillary reflectance component [third term on the right hand side of Eq. (8)]
    
    xlabel('Volumetric Water Content [m^{3} m^{-3}]','FontSize',12,'Color','k')
    ylabel('Transformed Reflectance [-]','FontSize',12)
    max_vertical = 2.5;
    axis([0 1.2 * max(theta_p) 0 max_vertical])
    title([soilname ' at ' num2str(wavenum) ' nm'])
    pbaspect([1 1 1])
    
    if j == 1
        legend({'Measurements','Proposed model Model','Second Term of Eq. (8)','Third Term of Eq. (8)'},'fontsize',12,'location','northwest')
        legend boxoff
    end
    
    % Fitted linear equations inside plots
    if j == 1 && linear_model_flag == 1
        text(0.05,1.3,txt,'FontSize', 12)
    elseif linear_model_flag == 1
        text(0.05,2.1,txt,'FontSize', 12)
    end 
    
    subplot(2, l, j + l)
    
    plot(theta_p,theta_a,'r--','Linewidth',2.0)
    hold on
    plot(theta_p,theta_c,'b-','Linewidth',2.0)
    
    xlabel('Volumetric Water Content [m^{3} m^{-3}]','FontSize',12,'Color','k')
    ylabel('Capillary and Adsorbed Water [m^{3} m^{-3}]','FontSize',12)
    axis([0 1.2 * max(theta_p) 0 1.2 * max(theta_p)])
    title(soilname)
    pbaspect([1 1 1])
    set(gcf,'units','inches','position',[0,0,12,10])
    
    
    
    if j == 1
        legend({'Adsorbed Water','Capillary Water'},'fontsize',12,'location','northwest')
        legend boxoff
    end
    
    
end







% -------------------------------------------------------------------------
% Measured spectral reflectance from the model versus measured at specified
% water content
% -------------------------------------------------------------------------
% Specifiy the suggested degree of saturations in which you want the plot
% Note that final plots would be somewhat different because measured water
% contents do not necessarily match the desired values.

fin_Wc = zeros(size(sug_Sw)); % Final water content, Adjusted
moistures_indexes = zeros(size(sug_Sw)); % Final indexes for plots

% Finding column indexes of desired Sw's
sug_Wc = sug_Sw .* theta_sat;
for k = 1:max(size(sug_Wc))
    theta_want = sug_Wc(k);
    ind = 0;
    
    
    if theta_want == 0
        ind = max(size(theta_meas));
        
    elseif theta_want == theta_sat
        ind = 1;
    else
        
        
        % Finding the closest water content in measurements
        for l = 1:max(size(theta_meas))
            dev = theta_want - theta_meas(l);
            
            if dev < 0
                ind = ind + 1;
            end
        end
    end
    
    if ind == 0
        ind = 1;
    else
        fin_Wc(k) = theta_meas(ind);
        moistures_indexes(k) = ind;
    end
    
end

fin_Sw = fin_Wc ./ theta_sat;


figure('name','Spectral Reflectance at Known Moisture')
for j = 1:max(size(moistures_indexes))
    
    moist_index = moistures_indexes(j);
    moisture = theta_meas(moist_index);
    
    theta_cap_mois = interp1(theta_p,theta_c,moisture,'spline'); % Capillary water at measured water contents of evaporation experiment
    theta_ads_mois = interp1(theta_p,theta_a,moisture,'spline'); % Adsorptive water at measured water contents of evaporation experiment
    
    itr = 0;
    
    for i = i_start:i_end
        itr = itr + 1;
        
        r_d = result(i,2);
        c_a = result(i,3);
        c_c = result(i,4);
        p_a = result(i,5);
        p_c = result(i,6);
        
        
        r_pre_m(itr) = r_d + (c_a .* theta_ads_mois .^ p_a) + (c_c .* theta_cap_mois .^ p_c); %#ok<SAGROW> % Predicted transformed reflectance [Eq. (8)]
%         r_pre_m(itr) = r_d + c_a.*(moisture - theta_cap).^p_a + c_c.*theta_cap.^p_c;
        
        
        % Transform to regular R
        R(itr) = 1 + r_pre_m(itr) - (r_pre_m(itr).^2 + 2.*r_pre_m(itr)).^0.5;
        wave_SWIR(itr) = Wave(i); %#ok<SAGROW>
        
    end
    
    plot(Wave(SWIR_range),refl(SWIR_range,moist_index),'r-','Linewidth',2.5)
    hold on
    plot(wave_SWIR,R,'b--','Linewidth',1.0)
    xlabel('Wavelength [nm]','FontSize',14,'Color','k')
    ylabel('Reflectance [-]','FontSize',14,'Color','k')
    axis([1100 2600 0 0.6])
    
    pbaspect([1.3 1 1])
    xticks([1300 1700 2100 2500])
    

    
end

