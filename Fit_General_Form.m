function dev = Fit_General_Form(x, r_d, rtran_meas, theta_ads, theta_cap)


c_a = x(1);
c_c = x(2);
p_a = x(3);
p_c = x(4);


r_pre = r_d + (c_a .* theta_ads.^ p_a) + (c_c .* theta_cap .^ p_c);

r_pre = real(r_pre);


dev = sum( (r_pre - rtran_meas).^2 );