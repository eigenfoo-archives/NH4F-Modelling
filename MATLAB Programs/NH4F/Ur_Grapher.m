%{
This code graphs the V-r curves for all pairs of atoms in a NH4F cluster.

Original parameters from Topper-Freeman (1994)
D value revisions from Topper Errata (Sept 2014)
Estimates for F- C, D values from Mao and Pappu (2012)
F- q value from Chandrasekhar, Spellmeyer and Jorgensen (1983)
C and D values for H-H calculated from Oobatake and Ooi (1972)
%}

d = 0.65;        %   Partial charges of H, N and F
q_H = 1-d;
q_N = -3+4*d;
q_F = -1;

%   Matrix of parameters for Topper-Freeman force field
%   A           alpha        C          D           q1          q2
P = [0,         0,           3.3384,    6.4895e+04, q_H,        q_H;  % H-H
    104.74,     1.5611,      25.393,    40,         q_N,        q_N;  % N-N
    0,          0,           7.1717,    1.2481e+06, q_F,        q_F;  % F-F
    0,          0,           9.2072,    1.6111e+03, q_H,        q_N;  % H-N
    0,          0,           11,    300, q_H,        q_F;  % H-F
    0,          0,           13.4948,   7.0656e+03, q_N,        q_F]; % N-F

%   Potentials (given by force field) for all atomic pairs
%{
    Arguments of UHF in Utot are calculated through trigonometry: given
    that the N-F distance is r, the H-F distances are r - 1.912 and
    r^2 - 1.276*r + 3.656. This calculation assumes that one H is oriented
    towards the F- ion, which is the case, as this is the most stable
    configuration.
%}
syms r;
UHH(r) = P(1, 1)*exp(-P(1, 2)*r) + P(1, 4)/r^12 + P(1, 5)*P(1, 6)/r - P(1, 3)/r^6;
UNN(r) = P(2, 1)*exp(-P(2, 2)*r) + P(2, 4)/r^12 + P(2, 5)*P(2, 6)/r - P(1, 3)/r^6;
UFF(r) = P(3, 1)*exp(-P(3, 2)*r) + P(3, 4)/r^12 + P(3, 5)*P(3, 6)/r - P(3, 3)/r^6;
UHN(r) = P(4, 1)*exp(-P(4, 2)*r) + P(4, 4)/r^12 + P(4, 5)*P(4, 6)/r - P(4, 3)/r^6;
UHF(r) = P(5, 1)*exp(-P(5, 2)*r) + P(5, 4)/r^12 + P(5, 5)*P(1, 6)/r - P(5, 3)/r^6;
UNF(r) = P(6, 1)*exp(-P(6, 2)*r) + P(6, 4)/r^12 + P(6, 5)*P(6, 6)/r - P(6, 3)/r^6;
Utot(r) = UNF(r) + UHF(r - 1.912) + 3*UHF(sqrt(r^2 - 1.276*r + 3.656));

%   Plotting graphs
hold on
fplot(@(r) UHH(r), [1, 10, -0.5, 2], 'y--')
fplot(@(r) UNN(r), [1, 10, -0.5, 2], 'm--')
fplot(@(r) UFF(r), [1, 10, -0.5, 2], 'c--')
fplot(@(r) UHN(r), [1, 10, -0.5, 2], 'r--')
fplot(@(r) UHF(r), [1, 10, -0.5, 2], 'g--')
fplot(@(r) UNF(r), [1, 10, -0.5, 2], 'b--')
fplot(@(r) Utot(r), [1, 10, -0.5, 2], 'k-')
fplot(@(r) 0, [1, 10, -0.5, 2], 'k-')
hold off
title('U(r) for (NH4)F Potential Pairs')
xlabel('r /Bohr')
ylabel('U /Hartree')
legend('H-H', 'N-N', ' F-F', 'H-N', 'H-F', 'N-F', 'NH4-F', 'y = 0')