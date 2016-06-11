%{
This code calculates the optimum interatomic distance and potential well
depth for a range of C-D (Lennard Jones parameters) pairs, and tabulates
the distances and depths, along with the parameters. The code also graphs
the data points on two separate graphs.

Original parameters from Topper-Freeman (1994)
D value revisions from Topper Errata (Sept 2014)
Estimates for F- C, D values from Mao and Pappu (2012)
F- q value from Chandrasekhar, Spellmeyer and Jorgensen (1983)
C and D values for H-H calculated from Oobatake and Ooi (1972)
%}

CVal = [];
DVal = [];
optLength = [];
wellDepth = [];
euclidean = [];

%   Partial charges of H, N and F
d = 0.5;
q_H = 1-d;
q_N = -3+4*d;
q_F = -1;

%   "Target" values of N-F distance and Utot potential well depth
targetLen = 2.205310229 + 1.012;
targetDepth = -0.23729574;

syms r;

for C = 1:5:71
    for D = 300:100:2000
        %   Matrix of parameters for Topper-Freeman force field
        %   C          D           q1          q2
        P = [C,        D,          q_H,        q_F;    % H-F
            13.4948,   7.0656e+03, q_N,        q_F];   % N-F

        %   Potentials (given by force field) for all atomic pairs
        UHF(r) = P(1, 2)/r^12 + P(1, 3)*P(1, 4)/r - P(1, 1)/r^6;
        UNF(r) = P(2, 2)/r^12 + P(2, 3)*P(2, 4)/r - P(2, 1)/r^6;
        Utot(r) = UNF(r) + UHF(r - 1.912) + 3*UHF(sqrt(r^2 - 1.276*r + 3.656));

        %   Calculate optimum N-F distance, well depth, and Euclidean distance between target values and calculated values.
        %   Append (C, D) values, N-F distance, well depth and Euclidean distance to appropriate array
        len = fminsearch(@(r) Utot(r), 5);
        depth = double(Utot(len));
        dist = sqrt((len-targetLen)^2 + (depth-targetDepth)^2);
        CVal = [CVal; C];
        DVal = [DVal; D];
        optLength = [optLength; len];
        wellDepth = [wellDepth; double(Utot(len))];
        euclidean = [euclidean; dist];
        fprintf('%d \t %d \t %d \t %d \t %d \n', C, D, len, depth, dist)
    end
end

%   Tabulate parameter values, N-F distance, well depth and Euclidean distance
T = table(CVal, DVal, optLength, wellDepth, euclidean, 'VariableNames', {'CVal', 'DVal', 'optLength', 'wellDepth', 'euclidean'})

%   Graph N-F distance, well depth and Euclidean distance all as functions of (C, D)
figure
scatter3(CVal, DVal, optLength, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c')
xlabel('C Value')
ylabel('D Value')
zlabel('Optimum Length /Bohr')

figure
scatter3(CVal, DVal, wellDepth, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
xlabel('C Value')
ylabel('D Value')
zlabel('Well Depth /Hartree')

figure
scatter3(CVal, DVal, euclidean, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'm')
xlabel('C Value')
ylabel('D Value')
zlabel('Euclidean Distance')