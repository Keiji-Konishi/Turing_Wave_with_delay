%
clear;
%
a=1.08; b=3.08;
D1=0.01; D2=0.1; D3=1.0;
%
L=16;
%
u1_eq=a; u2_eq=b/a; u3_eq=a; 
%
a11=-(2+b)+2*u1_eq*u2_eq;
a12=u1_eq*u1_eq;
a13=1;
a21=b-2*u1_eq*u2_eq;
a22=-u1_eq*u1_eq;
a23=0;
a31=1;
a32=0;
a33=-1;
%
ii=1;
%
for m=1:30
    gamma=m*m*(2*pi/L)*(2*pi/L);
    A=[a11 a12 a13; a21 a22 a23; a31 a32 a33];
    D=[D1 0 0; 0 D2 0; 0 0 D3];
    s_eig=eig(A-gamma*D);
    %
    Sig_eig(ii,1)=m;
    Sig_eig(ii,2)=gamma;
    Sig_eig(ii,3)=real(s_eig(1));
    Sig_eig(ii,4)=abs(imag(s_eig(1)));
    Sig_eig(ii,5)=real(s_eig(2));
    Sig_eig(ii,6)=abs(imag(s_eig(2)));
    Sig_eig(ii,7)=real(s_eig(3));
    Sig_eig(ii,8)=abs(imag(s_eig(3)));
    %
    ii=ii+1;
end
%
save Sig_eig_wave_num.dat Sig_eig -ascii
%

