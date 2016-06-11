%{
This code graphs the V-r curve for two ions of elementary charge.
The force field used is the Coulombic interaction term, and a repulsive
inverse-tenth-power term.

Refer to Physical Chemistry for the Chemical and Biological Sciences
by Raymond Chang. Page 671. Equations 16.5, 16.6.
%}

% All values in atomic units
k = 1;              % Coulomb constant
n = 10;             % Integer between 8 and 12
r_e = 4.45975;      % Equilibrium bond length of the ion pair
q = 1;              % Elementary charge
b = (k*q^2*r_e^(n-1)/n);
syms r;
V(r) = (-k*q^2)/r + b/(r^n);
ezplot(V(r), [0, 50, -0.25, 0.25])