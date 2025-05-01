%------
% Stability region in ep - tau plane: Coupled extended Brusselator model 
%------
%
clear
%
%== parameters ==
%
a=0.96;         %-- a --
b=2.85;         %-- b --
D1=0.01;        %-- D --
D2=0.1;
D3=1.0;
%
N_ep=100;
ep_min=0.0;         %-- minimum coupling strength --
ep_max=0.1;         %-- maximum coupling strength --
ep=linspace(ep_min,ep_max,N_ep+1);
%
N_tau=100;
tau_min=0.0;        %-- minimum delay time --
tau_max=6.0;        %-- maximum delay time --
tau=linspace(tau_min,tau_max,N_tau+1);
%
% -- log space --
N_gamma=1000;
gamma_min=1e-2;
gamma_max=1e+2;
gamma=logspace(log10(gamma_min),log10(gamma_max),N_gamma);
%
M = 20;%the discretization index (larger M=> smaller error
%
lambda_pm_max = zeros(1,N_gamma);
eig_data_ep_tau=zeros(0,4);
ii=1;
%
for i_ep=1:1:N_ep+1
    for i_tau=1:1:N_tau+1
        %
        for i_gamma=1:1:N_gamma             %-- varied gamma --
            %
            par = [a,b,D1,D2,D3,gamma(i_gamma),tau(i_tau),ep(i_ep)];
            lambda_m = eigAM('DDE_Wave_m',par,M);
            lambda_p = eigAM('DDE_Wave_p',par,M);
            %
            lambda_pm=[lambda_m;lambda_p];
            lambda_pm_max(i_gamma)=max(real(lambda_pm));
            %
        end
%
        [Real_max,I]=max(lambda_pm_max);
%
        eig_data_ep_tau(ii,1)=ep(i_ep);
        eig_data_ep_tau(ii,2)=tau(i_tau);
        eig_data_ep_tau(ii,3)=Real_max;
        eig_data_ep_tau(ii,4)=gamma(I);
%
        ii=ii+1;
%
    end
end
%
%save eig_data_ep_tau.dat eig_data_ep_tau -ascii;  
csvwrite('eig_data_ep_tau_a0.96_b285.csv',eig_data_ep_tau);
%