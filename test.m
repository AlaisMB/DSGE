clear   
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Money, New-Keynesian macroeconomics            %
%                    and the business cycle                   %
%                                                             %
%                  J.-O. Hairault, F. Portier                 %
%               European Economic Review, 1993                %
%                                                             %
% ------------------------------------------------------------%
%                     Replication Code by                     %
%    Ala?s Martin-Baillon, Antoine Mayerowitz, Arthur Stalla  %
%                          Fall 2016                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     Country Choice                        %%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Calibration and steady state                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%      |  US  | |  FR  |
param=[ 0.95     0.98  ;    % rho_an
        0.95     0.99  ;    % rho_ap
        0.009    0.009 ;    % epsilon_an
        0.008    0.0058;    % epsilon_ap
        1.014    1.028 ;    % g_barre
        0.377    0.50  ;    % rho_g
        0.009    0.0098;    % epsilon_g
        0.2      0.2   ;    % H
        0.3      0.335 ;    % alpha 
        0.025    0.0125;    % delta
        -1/3     -1/3  ;    % a
        1/9      1/9   ;    % sigma
        0.197    0.225 ;    % gamma
        0.1      0.25  ;    % MPY
        1        1     ;    % phi
        0.1      0.007 ;    % MC
        0.988    0.988 ;    % beta
        0        0     ;];  % zeta

if cntryid == 1 
	calib = num2cell(param(:,2));
else
	calib = num2cell(param(:,1)); 
end	

% Declare variables
[rho_an rho_ap epsilon_an epsilon_ap g_barre rho_g epsilon_g H alpha delta a sigma gamma MPY phi MC beta zeta] = calib{:};


rho_a = rho_an; % We suppose this for the moment. We should let the user decide.
A_barre = 1;
P = 0.3; 

z = (1/beta) - 1 + delta; 

KY = (alpha / (1 + gamma )) * (1/z) ; 
CY = 1 - delta * KY  ; 
YH =(A_barre)^(1/(1-sigma)) *  KY^(alpha/(1-alpha))  ;

Y = YH * H;
C = CY* Y ;
K = KY * Y ; 

w = ((1-alpha)/ (1+gamma))* YH ;
gamma_2 = - (1- g_barre * (1+z - delta)) * (MPY * (1/CY))^(1/sigma) ;
vP = gamma / (1 + gamma) ; 
v = vP * P ;
mu = gamma;

MP = MPY * Y ;

% Change in variables
Avar = (1+ sigma * ( a - 1)) / sigma;
Bvar = (1 - g_barre * (1 + z - delta ))* MP/C;
Cvar = (1 + sigma * ( a - 1))/ (sigma * (C^((sigma-1)/sigma) + gamma_2*(MP^((sigma-1)/sigma))));
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
%%		   	Linearized Model in Matrix Form 				 %%
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

M2(1,2) = - (Avar*Bvar)/(1-Bvar);
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

M2(4,4) = mu/(1+ mu);

% Equation 31
M1(5,:) = [0 0 0 -1 1 0 0];

M2(5,1) = 1;
M2(5,4) = mu/(1+mu);

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

% Equation 26
M4I(1,5) = 1;
M4L(1,5) = -1;

M5I(1,4) = -(1-beta*(1-delta));

% Equation 27
M4I(2,1) = Evar * K;
M4L(2,1) = Evar*K*(delta-1-z);

M5L(2,:) = [-Evar Evar*w*H Evar*w*H Evar*z*K 0 1 0];

% Equation 32
M4L(3,4) = Y/(1+mu);

M5I(3,7) = phi*beta*g_barre^3;
M5L(3,7) = -phi*g_barre^2;

% Equation 40
M4I(4,2) = 1;
M4L(4,3) = -1;

% Equation 40 bis magique
M4I(5,3) = -1;
M4L(5,3) = 1;

M5I(5,7) = 1;
M6I(5,2) = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Reduced Form                 	     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M40 = M4I - M5I*inv(M1)*M2;
M41 = M4L - M5L*inv(M1)*M2;


M50 = M6I + M5I*inv(M1)*M3;
M51 = M6L + M5L*inv(M1)*M3;

% Matrices of equation A45 in KP&R :
W = -inv(M40)*M41;
R = inv(M40)*M50;
Q = inv(M40)*M51;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      	Eigenvalues            			  %
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

if back == 2 && 5-(back + UR) == 3
    disp('Blanchard-Kahn Conditions are met !')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     	Saddle Path                			  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P1 = inv(P);

RHO = [rho_a 0 ; 0 rho_g];

% We decompose P1 into 3 matrices in order to decompose predetermined and jump
% variables and we'll get P1 = [Pss Psa ; Pas Paa]
Pss = P1(1:2,1:2);
Psa = P1(1:2,3:5);
Pas = P1(3:5,1:2);
Paa = P1(3:5,3:5);

% We split W into 2 systems of equation, W1 that depends on predetermined
% and W2 that depends on jump such that W = [W1;W2]
W1 = W(1:2,:);
W2 = W(3:5,:);

% We split W1 into the matrix that multiply s and the one that multiply a:

W11 = W1(:,1:2);
W12 = W1(:,3:5);

% Then we have the new system for s_t+1

PIs = W11 - W12*inv(Paa)*Pas;

% We now want to express the forward variables as a linear combination of
% shocks.

mu_m = EIG(3);
mu_mu = EIG(4);
mu_l = EIG(5);

% CL
C = P1*R*RHO+P1*Q;
C = C(3:5,:);

CL = [ -C(1,1)*1/(mu_m - rho_a)  -C(1,2)*1/(mu_m - rho_g);
       -C(2,1)*1/(mu_mu - rho_a) -C(2,2)*1/(mu_mu - rho_g);
       -C(3,1)*1/(mu_l - rho_a)  -C(3,2)*1/(mu_l - rho_g)];

% Computation of PIa
PIa = W12*inv(Paa)*CL+R(1:2,:)*RHO+Q(1:2,:);

M = [PIs PIa ; zeros(2) RHO];


%% PId matrix:

V = -inv(Paa)*Pas;
U = inv(Paa)*CL;

mu = [V(2,1:2) U(2,1:2)];
lambda = [V(3,1:2) U(3,1:2)];

% last 2 lines of PI
PID = zeros(9,4);
PID(8,:) = mu;
PID(9,:) = lambda;

X = inv(M1)*M2;
Y = inv(M1)*M3;
Z = X(:,3:5)*inv(Paa)*Pas;

G1 = X(:,1)-Z(:,1);
G2 = X(:,2)-Z(:,2);
G3 = X(:,3:5)*inv(Paa)*CL+Y;

PID(1:7,:) = [G1 G2 G3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IRF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
nrep = 10;

%Monetary shock
shock=[0;0;0;1];

for i = 1:nrep+1
    Rd(:,i) = (PID*M^(i-1))*shock;
    Rstate(:,i) = (M^(i-1))*shock;
end

k = Rstate(1,1:nrep)';
C = Rd(1,1:nrep)';
H = Rd(2,1:nrep)';
w = Rd(3,1:nrep)';
z = Rd(4,2:(nrep+1))';
Y = Rd(5,1:nrep)';
pi = Rd(6,1:nrep)';
f = Rd(7,1:nrep)';

figure
plot(1:nrep,[w,C,z],'-s');
legend('w','C','r','Location','southoutside');
title('French calibration - 1% shock on g')
xlabel('Quarters');
ylabel('relative deviation (%)');
set(gcf, 'Color', [1,1,1]);
hline = refline([0 0]);
set(hline,'LineStyle','--','color','black');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stochastic Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

nsimul = 5

nlong=10;

% Drawing shocks
aleaa = normrnd(0,epsilon_ap,nlong,nsimul);
aleag = normrnd(0,epsilon_g,nlong,nsimul);

% Drawing autoregressive shocks
for i=2:nlong
    epsa(i,:) = rho_ap*aleaa(i-1,:) + aleaa(i);
    epsg(i,:) = rho_g*aleag(i-1,:) + aleag(i);
end
for j=2:nsimul
    eps = [epsa(:,j) epsg(:,j);
    MS(:,1)=[0;0];
    for i=1:nlong
        MS(:,1) = M(1:2,1:2)*MS(:,i-1)+M(1:2,3)*eps(i,1)+M(1:2,4)*eps(i,2);
    end
end
