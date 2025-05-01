%------
% Max eigenvalue: Extended Brusselator model 
%------
%
clear
%
%== parameters ==
%
a=1.08;         %-- a --
b=3.08;         %-- b --
%a=0.96;        %-- a --
%b=2.85;        %-- b --
D1=0.01;        %-- D --
D2=0.1;
D3=1.0;
%
%
AA=[b-2 a^2 1; -b -a^2 0; 1 0 -1];  %-- A --
DD=[D1 0 0; 0 D2 0; 0 0 D3];        %-- D --
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
eig_data=zeros(0,3);
%
for ii=1:1:N             %-- varied gamma --
%
        Sys_A=AA-gamma(ii)*DD;
        eigen_Sys_A=eigs(Sys_A);
%
        [M,I]=maxk(real(eigen_Sys_A),3); %-- First, second, third largest real parts --
%
        eig_data(ii,1)=gamma(ii);
        eig_data(ii,2)=real(eigen_Sys_A(I(1)));       %-- First largest
        eig_data(ii,3)=abs(imag(eigen_Sys_A(I(1))));
        eig_data(ii,4)=real(eigen_Sys_A(I(2)));       %-- Second largest
        eig_data(ii,5)=abs(imag(eigen_Sys_A(I(2))));
        eig_data(ii,6)=real(eigen_Sys_A(I(3)));       %-- Third largest
        eig_data(ii,7)=abs(imag(eigen_Sys_A(I(3))));
%
end
%
%
save eig_data_a108_b308.dat eig_data -ascii;  
%save eig_data_a096_b285.dat eig_data -ascii;  
%