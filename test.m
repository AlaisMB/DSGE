clear   
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		The business cycle, J.-O. Hairault, F. Portier 	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  Choix du Pays %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cntrymap = containers.Map([1,2],{'France','USA'}); % Map the selected number with country name (To display country name later)
prompt= 'Type the number of the country you want to simulate and press Enter: \n 1:France[default] \n 2:USA \n'; % Ask the user to choose a country

% Stock the selected country and check that the user has entered a legit value :
cntryid = input(prompt);
if isempty(cntryid) % test if empty
    cntryid = 1;
elseif cntryid==1 || cntryid==2 % test if != 1 or 2
else
    cntryid = 1;
end
country = cntrymap(cntryid); % Stock the country name
fprintf('The model will be simulated for %s\n ',country);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration and steady state %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C_us = num2cell([ 0.95 ; 0.95 ; 0.009  ; 0.008  ;  1.014 ; 0.377 ; 0.009   ; 0.2  ; 0.3   ;  0.025 ; - 1/3 ; 1/9 ; 0.197 ; 0.1  ; 1 ; 0.01  ; 0.988; 0]);
C_fr = num2cell([ 0.98 ; 0.99 ;  0.009 ; 0.0058 ; 1.028  ;  0.50 ; 0.0098  ; 0.2  ; 0.335 ; 0.0125 ; - 1/3 ; 1/9 ; 0.225 ; 0.25 ; 1 ; 0.007 ; 0.988; 0]);

if cntryid == 1 
	calib = [C_fr];
else
	calib = [C_us]; 
end	

% Declare variables
[rho_an rho_ap epsilon_an epsilon_ap g_barre rho_g epsilon_g H alpha delta a sigma gamma MPY phi MC beta zeta] = calib{:};


rho_a = rho_an; % We suppose this for the moment. We should let the user decide.

A_barre = 0.2;
P = 0.3; 

z = (1/beta) - 1 + delta; 

KY = (alpha / (1 + gamma )) * (1/z) ; 
CY = 1 - delta * KY  ; 
YH = (A_barre)^(1/(1-sigma)) * KY^(alpha/(1-alpha)) * H ;

Y = YH * H;
C = CY* Y ;
K = KY * Y ; 

w = ((1-alpha)/ (1+gamma))* YH ;
gamma_2 = - (1- g_barre * (1+z - delta)) * (MPY * (1/CY))^(1/sigma) ;


vP = gamma / (1 + gamma) ; 
v = vP * P ;
mu = log(vP); % ?????? Equation 30 et 31, d'o? sort mu ?!

MP = MPY * Y ;

% Change in variables
Avar = (1- sigma * ( a - 1)) / sigma;
Bvar = (1 - g_barre * (1 - z - delta ))* MP/C;
Cvar = (1 - sigma * ( a - 1))/ (sigma * (C^((sigma-1)/sigma) + gamma_2*(MP^((sigma-1)/sigma))));
Dvar = (g_barre / beta ) - 1;
Evar = (K*(delta - z) - w * H + C)^(-1);


%%%%%%%%%%%% Display parameters value %%%%%%%%%%%%

disp(' ');
param1=[alpha beta delta H rho_ap rho_g gamma sigma];
param2=[Y     C     K     w   v   z   MP ];
disp('    alpha     beta  	delta 	   H   	   rho_a    rhog     gamma      sigma'); 
disp(' ');
disp(param1);
disp(' ');
disp('    Y         C   	K  	  w         v  	      z	       MP ');
disp(param2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		   	Linearized Model in Matrix Form 				  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M1=zeros(7,7); 
M2=zeros(7,5);
M3=zeros(7,2);


M4I = zeros(5,5);
M4L = zeros(5,5);
M5I = zeros(5,7); 
M5L = zeros(5,7); 
M6I = zeros(5,2);
M6L = zeros(5,2);


%%%%%%%%%%%% STATIC EQUATIONS %%%%%%%%%%%%

% Equation 23
M1(1,1) = (Avar/(1-Bvar)) - 1/sigma;
M1(1,7) = (Avar*Bvar)/(1-Bvar);

M2(1,2) = (Avar*Bvar)/(1-Bvar);
M2(1,5) = 1;

% Equation 24
M1(2,2) = (zeta * H)/(1-H);
M1(2,3) = 1;

M2(2,5) = -1;

% Equation 25
M1(3,1) = -(Dvar*Cvar*C^((sigma-1)/sigma));
M1(3,4) = g_barre*z;
M1(3,7) = (Dvar*((sigma-1)/sigma) + 1 - Cvar);

M2(3,2) = (-1/sigma)*Dvar + Cvar;
M2(3,5) = -Dvar;

% Equation 30
M1(4,:) = [0 -1 -1 0 1 0 0];

M2(4,4) = mu/(1- mu);

% Equation 31
M1(5,:) = [0 0 0 -1 1 0 0];

M2(5,1) = 1;
M2(5,4) = mu/(1-mu);

% Equation 34
M1(6,2) = alpha - 1;
M1(6,5) = 1;

M2(6,1) = alpha;

M3(6,1) = 1;

% Equation 35
M1(7,2) = H*w;
M1(7,3) = H*w;
M1(7,4) = z*K;
M1(7,5) = -Y;
M1(7,6) = -w*H - z*K;

M2(7,1) = -z*K;

%%%%%%%%%%%% DYNAMICS EQUATIONS %%%%%%%%%%%%

% Equation 2 (Technology process)
%M6I(6,1) = 1;
%M6L(6,1) = -rho_a;

% Equation 8 (Money growth process)
%M6I(7,2) = 1;
%M6L(7,2) = -rho_g;

% Equation 26
M4I(1,5) = 1;
M4L(1,5) = -1;

M5I(1,4) = -(1-beta*(1-delta));

% Equation 27
M4I(2,1) = Evar * K;
M4L(2,1) = Evar*K*(delta-1-z);

M5L(2,:) = [-Evar Evar*w*H Evar*w*H Evar*w*K 0 1 0];

% Equation 32
M4L(3,4) = Y/(1+mu);

M5I(3,7) = phi*beta*g_barre^3;
M5L(3,7) = -phi*g_barre^2;

% Equation 40
M4I(4,3) = -1;
M4L(4,3) = 1;

M5I(4,7) = 1;
M6I(4,2) = -1;

% Equation 40 bis magique
M4I(5,2) = 1;
M4L(5,2) = -1;

M5L(5,7) = 1;
M6L(5,2) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      	Reduced Form            			  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M40 = M4I - M5I*inv(M1)*M2;
M41 = M4L - M5L*inv(M1)*M2;

M50 = M6I + M5I*inv(M1)*M3;
M51 = M6L + M5L*inv(M1)*M3;

W = -inv(M40)*M41;
R = inv(M40)*M50;
Q = inv(M40)*M51;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      	Eigenvalues             			  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[P,D]=eigs(W); % eigs() function gives eigenvalues in increasing order
EIG = eigs(W);

UR = numel(EIG(abs(EIG)==1));
back = numel(EIG(abs(EIG)<1));

if UR ~= 0,
    disp('Unit Root !')
else
    disp('No unit root :) ')
end

disp('Backward variables ')
disp(back)
disp('Forward variables  ')
disp(numel(EIG)-(back + UR))
pause

disp('Eigenvalues:')
disp(EIG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      	Saddle Path                			  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P1 = inv(P);
