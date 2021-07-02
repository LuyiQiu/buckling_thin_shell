clear;
clc;
%% For balloon
a0=[]; %initial deformation function
A=[]; 
S=struct([]);

% set constants for E.coli
E=2.5*10^7;  
Ec=6*10^7; 
h=5e-9;
nu=0.2;
R=5e-7;
D=Ec*h^3/(12*(1-nu^2)); 
alpha=sqrt(D/(E*h*R^2));

% Calculate the deformation of balloon
phat=0.1; %phat in balloon normalization, defined in the main text
khat=0.01; %phat in balloon normalization, defined in the main text
a0=[0,-khat^2/8/(alpha^2+phat/3),0,0,0];
s=E_minimization(phat,khat,a0,alpha);% change the energy function in E_minimization to balloon energy expression
% s describing the deformation is the result.

% Calculate the deformation of shell
pbar=2.5*p_i-2.5; %pbar in shell normalization
kbar=k_i*0.2; %kbar in shell normalization
a0=[0,-kbar^2/8/(1+pbar/3),0,0,0];
s=E_minimization(pbar,kbar,a0,alpha); % change the energy function in E_minimization to shell energy expression
% s describing the deformation is the result.

