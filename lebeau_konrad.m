function dev = lebeau_konrad(x, head, theta_meas, theta_sat_vol)

theta_0 = x(1);
h_median = x(2);
sigma = x(3); 

h_dry = -10^5;

theta_c = theta_sat_vol.*(0.5.*erfc( log(head./h_median)./(sqrt(2).*sigma) ));

theta_a = theta_0.*(1 - log(abs(head))./log(abs(h_dry)) ) .* (1 - theta_c./theta_sat_vol);


theta_p = theta_c + theta_a; 

dev = sum( (theta_p - theta_meas).^2 );