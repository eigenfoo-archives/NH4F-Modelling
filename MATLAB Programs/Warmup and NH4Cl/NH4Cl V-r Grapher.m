%{
This code prompts the user for a pair of atoms within a NH4Cl nanocluster,
and displays the V-r graph for the selected pair.

Force field parameters from Topper-Freeman (1994)
with revisions from Topper Errata (Sept 2014)
%}

%   A           alpha        C          D           q1          q2
P = [1.0162,    1.9950,      2.9973,    100,        0.35,       0.35;   % H-H
    104.74,     1.5611,      25.393,    40,         -0.40,      -0.40;  % N-N
    125.55,     1.7489,      113.68,    800,        -1.00,      -1.00;  % Cl-Cl
    10.318,     1.7780,      8.7229,    80,         0.35,       -0.40;  % H-N
    0,          0,           10.033,    43884.0,    0.35,       -1.00;  % H-Cl
    114.22,     1.6550,      53.736,    200,        -0.40,      -1.00]; % N-Cl

syms r;
UHH(r) = P(1, 1)*exp(-P(1, 2)*r) + P(1, 4)/r^12 + P(1, 5)*P(1, 6)/r - P(1, 3)/r^6;
UNN(r) = P(2, 1)*exp(-P(2, 2)*r) + P(2, 4)/r^12 + P(2, 5)*P(2, 6)/r - P(1, 3)/r^6;
UClCl(r) = P(3, 1)*exp(-P(3, 2)*r) + P(3, 4)/r^12 + P(3, 5)*P(3, 6)/r - P(3, 3)/r^6;
UHN(r) = P(4, 1)*exp(-P(4, 2)*r) + P(4, 4)/r^12 + P(4, 5)*P(4, 6)/r - P(4, 3)/r^6;
UHCl(r) = P(5, 1)*exp(-P(5, 2)*r) + P(5, 4)/r^12 + P(5, 5)*P(1, 6)/r - P(5, 3)/r^6;
UNCl(r) = P(6, 1)*exp(-P(6, 2)*r) + P(6, 4)/r^12 + P(6, 5)*P(6, 6)/r - P(6, 3)/r^6;

sprintf('1 = H-H\n2 = N-N\n3 = F-F\n4 = H-N\n5 = H-F\n6 = N-F\n0 = all graphs')
in = input('Input vector of desired graphs:\n')

if find(in == 0)
    in = [1,2,3,4,5,6];
end
if find(in == 1)
    ezplot(UHH(r), 0, 1)
    hold on
end
if find(in == 2)
    ezplot(UNN(r), 0, 1)
    hold on
end
if find(in == 3)
    ezplot(UClCl(r), 0, 1)
    hold on
end
if find(in == 4)
    ezplot(UHN(r), 0, 1)
    hold on
end
if find(in == 5)
    ezplot(UHCl(r), 0, 1)
    hold on
end
if find(in == 6)
    ezplot(UNCl(r), 0, 1)
    hold on
end