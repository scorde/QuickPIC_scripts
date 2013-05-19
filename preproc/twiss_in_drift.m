

% Input parameters


s_waist = 0.;  % in m

sig_x = 30;  % in um
sig_y = 30;  % in um

eps_x = 50;  % in um
eps_y = 50;  % in um

gamma = 39139;



% Output parameters

beta_star_x = 1e-6*gamma*sig_x^2/eps_x;  % in m
beta_star_y = 1e-6*gamma*sig_y^2/eps_y;  % in m

alpha_x = s_waist/beta_star_x
beta_x = beta_star_x*(1+alpha_x^2)  % in m

alpha_y = s_waist/beta_star_y
beta_y = beta_star_y*(1+alpha_y^2)  % in m


