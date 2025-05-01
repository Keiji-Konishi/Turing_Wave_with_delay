%------
% Max eigenvalue for delayed systems: Extended Brusselator model 
%------
%
clear
%
%== parameters ==
%
%a=1.08;        %-- a --
%b=3.10;        %-- b --
a=0.96;         %-- a --
b=2.85;         %-- b --
D1=0.01;        %-- D --
D2=0.1;
D3=1.0;
%
tau=1.0;      %-- \tau --
ep=0.05;      %-- \ep --
%
N=5000;
%
% -- log space --
gamma_min=1e-2;
gamma_max=1e+2;
gamma=logspace(log10(gamma_min),log10(gamma_max),N);
% -- linear space --
%gamma_min=0;
%gamma_max=1e+2;
%gamma=linspace(gamma_min,gamma_max,N+1);
%
M = 20;%the discretization index (larger M=> smaller error
%
eig_data=zeros(0,3);
v_m = zeros(1,N);
v_p = zeros(1,N);
%
for ii=1:1:N             %-- varied gamma ‚ð•Ï‰» --
%
        par = [a,b,D1,D2,D3,gamma(ii),tau,ep];
        [lambda_m,MM_m] = eigAM('DDE_Wave_m',par,M);
        [lambda_p,MM_p] = eigAM('DDE_Wave_p',par,M);
%
        lambda_pm=[lambda_m;lambda_p];
        [v_pm,I]=max(real(lambda_pm));
%
        eig_data(ii,1)=gamma(ii);
        eig_data(ii,2)=v_pm;                     %-- Real part: Largest real
        eig_data(ii,3)=abs(imag(lambda_pm(I)));  %-- Imaginary part: Largest real
end
%
%
%save eig_data_a108_b310.dat eig_data -ascii;  
save eig_data_a096_b285_tau100_ep005.dat eig_data -ascii;  
%