function [ G_stretch ] = G_stretch(  fc,bp, T )
% ��֪force,critical force����DNA����ת�����force�������칦G_stretch,��λ��pN��nm��

% �������¹�ʽ����G4��������칦��x_G4( f ) = l_0*coth( f*l0/kBT) ? kBT/f. 
l_0 = 1.7;

kBT = 1.3806504e-2*( T+273.15 ); 
x_G4 = @(f)l_0*coth(f*l_0./kBT) - kBT./f;
G_G4 = integral(x_G4,0,fc);

%����ssDNA�ľ��鹫ʽ����ssDNA��������칦
% �ù�ʽ��ʹ��ʱ��Ҫ���ֶ�����bp��
h =0.34;
K=0.1;  %K+ Ũ��, ��λΪmol
a1=0.21;   a2=0.34;   a3=2.1*log(K/0.0025)/log(0.15/0.0025);
f1=0.0037;  f2=2.9;
f3=8000;
% bp = 21;
ssDNA_length =@(f) bp*h*(a1*log(f/f1)/(1+a3*exp(-f/f2))-a2-f/f3);
G_ss = integral(ssDNA_length,0,fc);
%��������
G_stretch =(G_ss - G_G4);



end

