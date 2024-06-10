% Code to fit the Lebeau and Konrad (2010) model to soil retention 
% measurements. Please feel free to contact me via email: sarem.nrz@gmail.com
clear
clc
close all


%--------------------------------------------------------------------------
% Reading the Measured values of the soil water retention curve
%--------------------------------------------------------------------------
flag_optimization = 1; % Put 1 if a new optimization is needed
soilname = 'AZ12';
version = 'Lebeau_Konrad';


load([soilname '_Measured_SWRC_Data'],'SWRC_Data') % .mat file containing the measured SWRC data

head_cm = - SWRC_Data(:,1);  % Measured  matric potential [cm]
head_meas = head_cm./100; % Measured  matric potential [m]

theta_meas = SWRC_Data(:,2); % Volumetric water content
theta_sat_vol = max(theta_meas); % Saturated water content


%--------------------------------------------------------------------------
% Optimization for finding the Lebeau and Konrad (2010) parameters
%--------------------------------------------------------------------------
if flag_optimization == 1 
    
    rng default 
    fun = @(x)lebeau_konrad(x, head_meas, theta_meas, theta_sat_vol);
    
    gs = GlobalSearch('FunctionTolerance',1e-20,'NumTrialPoints',2000,'XTolerance',1e-20,'Display','final');
    problem = createOptimProblem('fmincon','x0',[0.1, -1.0, 1],...
        'objective',fun,'lb',[eps, -10, eps],'ub',[0.46, -eps, 3]);
    
    
    ms = MultiStart(gs);
    %x = run(gs,problem);
    x = run(ms,problem,10);
    
    save(['All_AZ_Soils_' soilname '_' version '.mat'],'x')
else
    load(['All_AZ_Soils_' soilname '_' version '.mat'],'x')
    
end


% Solution
theta_0 = x(1);
h_median = x(2);
sigma = x(3);


% Refining the SWRC for plots
a = -2; b = 5; % Range of matric potential: 10^a to 10^b
head = -logspace(a,b,100); % A refined and Logarithmically distributed vector of matric potential 
h_dry = -10^5; % Matric potential at oven dryness

theta_c = theta_sat_vol.*(0.5.*erfc( log(head./h_median)./(sqrt(2).*sigma) )); % Capillary component
theta_a = theta_0.*(1 - log(abs(head))./log(abs(h_dry)) ) .* (1 - theta_c./theta_sat_vol); % Adsorptive component
theta_p = theta_c + theta_a; % Total predicted water content


%--------------------------------------------------------------------------
% Plot of SWRC and Water Components
%--------------------------------------------------------------------------
figure('name','Labeau-Konrad SWRC & Water Components')

% SWRC plot
subplot(2,1,1)
max_theta_label = 0.45; % Maximum of horizontal axis

semilogy(theta_meas,-head_meas,'ko','LineWidth',1.2,'MarkerSize',8); % Measured SWRC
hold on
semilogy(theta_p,-head,'k-','LineWidth',2.5); % Fitted SWRC 
semilogy(theta_c,-head,'b-','LineWidth',1.5); % Capillary component
semilogy(theta_a,-head,'r--','LineWidth',1.5); % Adsorptive component

xlabel('Volumetric Water Content [m^{3} m^{-3}]','FontSize',13,'Color','k')
ylabel('Matric Potential [-m] ','FontSize',13,'Color','k')
legend({'Measurements','Fitted SWRC','Capillary Component','Adsorptive Component'},'fontsize',12)
axis([0 max_theta_label 10^a 10^b])
xticks([0 .1 .2 .3 .4])
yticks([.01 1 100 10000])
title([soilname ' ' version])
pbaspect([1.5 1 1])
legend boxoff
text(0.405,10,'(a)','Color','k','FontSize',18,'FontWeight','bold','FontName','Times','HorizontalAlignment','left')

% Water components versus total water content
subplot(2,1,2)
plot(theta_p,theta_c,'b-','Linewidth',1.5) % Capillary water
hold on
plot(theta_p,theta_a,'r--','Linewidth',1.5) % Adsorbed water
xlabel('Volumetric Water Content [m^{3} m^{-3}]','FontSize',13,'Color','k')
ylabel('Capillary and Adsorbed Water [m^{3} m^{-3}]','FontSize',13)
legend({'Capillary Water', 'Adsorbed Water'},'fontsize',12,'Location','northwest')
axis([0 max_theta_label 0 max_theta_label])
xticks([0 .1 .2 .3 .4])
yticks([0 .1 .2 .3 .4])
%title(soilname,'FontSize',14)
pbaspect([1.5 1 1])
set(gcf,'units','inches','position',[6,2,6,9])
legend boxoff
text(0.405,0.2,'(b)','Color','k','FontSize',18,'FontWeight','bold','FontName','Times','HorizontalAlignment','left')


