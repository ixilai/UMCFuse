function [MPC] = MonogenicPC(im, filterbankargs)
% ======================================================================
% Computes phase congruency using monogenic signals in an image.
%
% Reference:
% 
%    Luo, X.-G., Wang, H.-J., Wang, S.: Monogenic signal theory based
%    feature similarity index for image quality assessment. AEU-International
%    Journal of Electronics and Communications 69(1), 75-81 (2015)
%-----------------------------------------------------------------------

nscale = filterbankargs.nscale;
epsilon = 0.0001;
Rx = cell(1,nscale);
Ry = cell(1,nscale);
xi = 1.5;

[I, T, cutoff_gamma, cutoff_s] = loggaborconv2(im, filterbankargs);
% Get image components in different bands.

for s = 1:nscale
    [Rx{s}, Ry{s}] = RieszTransform(I{s});
end
An = cell(1,nscale); % Local amplitude
En = cell(1,nscale); % Local energy
A = 0; % Sum of local amplitude
E = 0; % Sum of local energy
for s = 1:nscale
    An{s} = sqrt(Rx{s}.^2 + Ry{s}.^2 + I{s}.^2);
    A = A + An{s};
    En{s} = (Rx{s}).^2 + (Ry{s}).^2 + (I{s}).^2;
    E = E + En{s};
end
E = sqrt(E);

N = 2;
c = A./(N*(max(max(A)) + epsilon));
W = 1./(1 + exp(cutoff_gamma*(cutoff_s - c)));
% Calculate the weight function with regard to the scale of the filter bank.

T_alpha = 2;
T = T.*T_alpha;
f1 = -(1 - xi*acos(E./A));
f2 = E - T;
f1(f1<0) = 0;
f2(f2<0) = 0;
% Noise compensation for sharpness of the MPC.

MPC = W.*f1.*f2./(A + epsilon);
MPC = MPC/sum(sum(MPC));
end