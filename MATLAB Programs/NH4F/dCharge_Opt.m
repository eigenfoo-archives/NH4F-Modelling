%{
This code calculates the optimum interatomic distance for a range of
ppossible partial charges. The distances are tabulated, along with the
partial charges.

Original parameters from Topper-Freeman (1994)
D value revisions from Topper Errata (Sept 2014)
Estimates for F- C, D values from Mao and Pappu (2012)
F- q value from Chandrasekhar, Spellmeyer and Jorgensen (1983)
C and D values for H-H calculated from Oobatake and Ooi (1972)
%}

dVal = [];
HCharge = [];
NCharge = [];
optLen = [];

q_F = -1;
syms r;

for d = 0:0.05:1
    %   Partial charges of H and N (partial charge of F is always -1)
    q_H = 1-d;
    q_N = -3+4*d;

    %   Matrix of parameters for Topper-Freeman force field
    %   C          D           q1          q2
    P = [4.8931,   2.8460e+05, q_H,        q_F;    % H-F
        13.4948,   7.0656e+03, q_N,        q_F];   % N-F

    %   Potentials (given by force field) for relevant atomic pairs
    UHF(r) = P(1, 2)/r^12 + P(1, 3)*P(1, 4)/r - P(1, 1)/r^6;
    UNF(r) = P(2, 2)/r^12 + P(2, 3)*P(2, 4)/r - P(2, 1)/r^6;
    Utot(r) = UNF(r) + UHF(r - 1.912) + 3*UHF(sqrt(r^2 - 1.276*r + 3.656));
    
    %   Append partial charges and corresponding N-F distance to appropriate array
    dVal = [dVal; d];
    HCharge = [HCharge; 1-d];
    NCharge = [NCharge; -3+4*d];
    optLen = [optLen; fminsearch(@(r) Utot(r), 5)];
end

%   Tabulate partial charges and N-F distance
T = table(dVal, HCharge, NCharge, optLen, 'VariableNames', {'dVal', 'HCharge', 'NCharge', 'optLength'})