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
%    Alais Martin-Baillon, Antoine Mayerowitz, Arthur Stalla  %
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


rho_a = rho_ap; % We suppose this for the moment. We should let the user decide.
epsilon_a = epsilon_ap;
A_barre = 1; 

z = (1/beta) - 1 + delta; 

KY = (alpha / (1 + gamma )) * (1/z) ; 
CY = 1 - delta * KY  ; 

Y = KY^(alpha/(1-alpha))*H;
C = CY* Y ;
K = KY * Y ; 

w = ((1-alpha)/ (1+gamma))* (Y/H) ;
gamma_2 = - ( 1 - g_barre * (1 + z - delta) ) * ( MPY / CY )^( 1/sigma ) ;
vP = gamma / (1 + gamma) ;  
mu = gamma;

MP = MPY * Y ;
MPY=0.25;

% Change in variables
Avar = (1+ sigma * ( a - 1)) / sigma;
Bvar = (1 - g_barre * (1 + z - delta ))* MP/C;
Cvar = (1 + sigma * ( a - 1))/ (sigma * (C^((sigma-1)/sigma) + gamma_2*(MP^((sigma-1)/sigma))));
Dvar = (g_barre / beta ) - 1;
Evar = (K*(delta - z) - w * H + C)^(-1);


%%%%%%%%%%%% Display parameters value %%%%%%%%%%%%

disp(' ');
param1=[alpha beta delta H rho_ap rho_g gamma sigma];
param2=[Y     C     K     w      z   MP ];
disp('    alpha     beta  	delta 	   H   	   rho_a    rhog     gamma      sigma'); 
disp(' ');
disp(param1);
disp(' ');
disp('    Y         C   	K  	  w       z	       MP ');
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

M2(1,2) = (Avar*Bvar)/(1-Bvar);
M2(1,5) = 1;

% Equation 24 
M1(2,2) = (zeta * H)/(1-H);
M1(2,3) = 1;

M2(2,5) = -1;

%Equation 40
M1(3,7) = 1;

M2(3,:) = [0 1 -1 0 0];

M3(3,2) = 1;

% Equation 30
M1(4,:) = [0 1 1 0 -1 0 0];

M2(4,4) = -mu/(1+ mu);

% Equation 31
M1(5,:) = [0 0 0 1 -1 0 0];

M2(5,1) = -1;
M2(5,4) = - mu/(1+mu);

% Equation 34
M1(6,2) = alpha - 1;
M1(6,5) = 1;

M2(6,1) = alpha;

M3(6,1) = 1;

% % Equation 35
M1(7,2) = H*w;
M1(7,3) = H*w;
M1(7,4) = z*K;
M1(7,5) = -Y;
M1(7,6) =Y -w*H - z*K;

M2(7,1) = -z*K;


%%%%%%%%%%%% DYNAMICS EQUATIONS %%%%%%%%%%%%

% Equation 26 -> OK
M4I(1,5) = 1;
M4L(1,5) = -1;

M5I(1,4) = -1 - beta * ( delta - 1);

% Equation 27

M4I(2,1) = Evar * K;
M4L(2,1) = Evar*K*(delta-1-z);

M5L(2,:) = [-Evar*C Evar*w*H Evar*w*H Evar*z*K 0 1 0];

% Equation 32
M4L(3,4) = Y/(1+mu);

M5I(3,7) = phi*beta*g_barre^3;
M5L(3,7) = -phi*g_barre^2;

%Equation 25 
M5I(5,1) = -( Dvar*Cvar*C^( (sigma-1)/sigma ) );
M5I(5,4) = g_barre*z;
M5I(5,7) = (Dvar + 1 +Dvar*( (-1/sigma) + gamma_2 * Cvar*(MP)^((sigma-1)/sigma)) );

M4I(5,2) = Dvar* ( (-1/sigma) + gamma_2 * Cvar*(MP)^((sigma-1)/sigma) );
M4I(5,5) = -Dvar;


% Equation 40 bis magique
M4I(4,2) = 1;
M4L(4,3) = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Reduced Form                 	     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M40 = M4I - M5I*inv(M1)*M2;
M41 = M4L - M5L*inv(M1)*M2;


M50 = M6I + M5I*inv(M1)*M3;
M51 = M6L + M5L*inv(M1)*M3;

% Reduces form model:
% S(t+1) = ( W * S(t) ) + ( R * E(t+1) ) + ( Q * E(t) )
% with S the vector of state variables and E the matrix of exogenous
% variables.

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
    disp('Unit Root(s) !')
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

% The purpose of this section is to substitute the equations for which
% |eigenvalue| > 1.
% In order to achieve this we must solve forward the equations of m(t+1),
% mu(t), and lambda(t).
% The following code follows the method discussed by King, Plosser & Rebelo and
% Cochrane and is based on the seminal work of Blanchard and Kahn.
%
% We try to keep the same notation as in the original article to keep the
% reading understandable.

% We are looking for a solution of the form S(t+1) = M*S(t) <=>
%
% | k(t+1) |   |    ?_s        ?_e      |   | k(t) |
% | m(t+1) |   |   (2x2)      (2x2)     |   | m(t) |
% |        | = |                        | * |      |
% | A(t+1) |   |  0     0   ?(A)    0   |   | A(t) |
% | g(t+1) |   |  0     0   0     ?(g)  |   | g(t) |

P1 = inv(P);

RHO = [rho_a 0 ; 0 rho_g];

% -------------------------------------------
% PIs computation

% We decompose P1 into 3 matrices in order to decompose predetermined and jump
% variables and we'll get
%             |             |
%             |  Pss   Psa  |
%             | (1x2) (1x2) |
%   inv(P)=   |             | 
%             |  Pas   Paa  |
%             | (1x2) (1x2) |
%             |             |
Pss = P1(1:2,1:2);
Psa = P1(1:2,3:5);
Pas = P1(3:5,1:2);
Paa = P1(3:5,3:5);

% We split W into 2 systems of equation. W1 is the equations describing the
% dynamics of [k(t+1) , m(t+1)] and W2 that describes the dynamics of the costate variables.
%         |   W1  |
%         | (2x5) |
%    W  = |       |
%         |   W2  |
%         | (3,5) |
W1 = W(1:2,:);
W2 = W(3:5,:);

% We split W1 to isolate the dynamics of [k(t+1) ; m(t+1)] that depends
% on its past values (W11) and (W12) that depends on exogenous variables.
%                 |   W11    W12  |   | k(t)  |
%                 |  (2x2)  (2x3) |   | m(t)  |
%    W * ?(t) =   |               | * | m(t+1)|
%                 |       W2      |   | A(t)  |
%                 |     (3,5)     |   | g(t)  |
W11 = W1(:,1:2);
W12 = W1(:,3:5);

% We can now compute ?_s :
PIs = W11 - W12*inv(Paa)*Pas;

% -------------------------------------------
% PIe computation

% We now want to express the forward variables as a linear combination of
% shocks.

% We save the eigenvalue corresponding to each forward variable.
mu_m = EIG(3);
mu_mu = EIG(4);
mu_l = EIG(5);

% Here we compute a critical matrix to solve the model.
% CL expresses the forward variables in term of shocks.
% It is a direct application of the transversality condition to the model
% in order to keep it converging.

% To simplify our resolution we can express R*E(t+1) + Q*E(t) = C*E(t).
% This is possible because we suppose AR() process on exogenous variables
C = P1*R*RHO+P1*Q;
C = C(3:5,:); % We just want to keep the dynamic of our forward variables (i.e. last 3 equations)


% The computation of CL is given by :
CL = [ -C(1,1)*1/(mu_m - rho_a)  -C(1,2)*1/(mu_m - rho_g);
       -C(2,1)*1/(mu_mu - rho_a) -C(2,2)*1/(mu_mu - rho_g);
       -C(3,1)*1/(mu_l - rho_a)  -C(3,2)*1/(mu_l - rho_g)];

% We can now substitute CL and find our ?_e matrix
PIe = (W12 * inv(Paa)*CL) + ( R(1:2,:) * RHO ) + Q(1:2,:);

% We can fully identify our state variables dynamic
M = [PIs PIe ; zeros(2) RHO];


% -------------------------------------------
% Control variables dynamics

%% PId matrix:

V = -inv(Paa)*Pas;
U = inv(Paa)*CL;

mu = [V(2,1:2) U(2,1:2)];
lambda = [V(3,1:2) U(3,1:2)];

% last 2 lines of PI
PID = zeros(9,4);
PID(8,:) = mu;
PID(9,:) = lambda;

XX = inv(M1)*M2;
YY = inv(M1)*M3;
ZZ = XX(:,3:5)*inv(Paa)*Pas;

G1 = XX(:,1)-ZZ(:,1);
G2 = XX(:,2)-ZZ(:,2);
G3 = XX(:,3:5)*inv(Paa)*CL+YY;

PID(1:7,:) = [G1 G2 G3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Impulse Response Function         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% -------------------------------------------
% Initialization

nrep = 12; % Periods

% 1% shock
%       k   m   A   g
shock=[ 0 ; 0 ; 0 ; 1 ];

% --------------------------------------------
% Compute state and control variables

for i = 1:nrep+1
    Rd(:,i) = (PID*M^(i-1))*shock; % Control variables response to shock.
    Rstate(:,i) = (M^(i-1))*shock; % State variables response to shock.
end

% --------------------------------------------
% Extract series from previous matrices

series = {'C' 'H' 'w' 'z' 'Y' 'pi' 'f' 'mu'};

for i = 1:length(series)
    irf.( series{i} ) = Rd( i , 1:nrep )';
end
irf.k  = Rstate(1,1:nrep)';
irf.m  = Rstate(2,1:nrep)';
irf.g  = Rstate(4,1:nrep)';
irf.z  = Rd(4,2:(nrep+1))';
irf.I  = (irf.Y-CY*irf.C)/CY ;

% --------------------------------------------
% Plot figures

figure
plot(1:nrep,[irf.w,irf.C,irf.z],'-o');
legend('w','C','r','Location','southoutside');
xlabel('Quarters');
ylabel('relative deviation (%)');
set(gcf, 'Color', [1,1,1]);
hline = refline([0 0]);
set(hline,'LineStyle','--','color','black');

figure
plot(1:nrep,[irf.pi,irf.mu],'-o');
legend('pi','mu','Location','southoutside');
xlabel('Quarters');
ylabel('relative deviation (%)');
set(gcf, 'Color', [1,1,1]);

figure
plot(1:nrep,[irf.H,irf.Y,irf.I],'-o');
legend('H','Y','I','Location','southoutside');
xlabel('Quarters');
ylabel('relative deviation (%)');
set(gcf, 'Color', [1,1,1]);
hline = refline([0 0]);
set(hline,'LineStyle','--','color','black');

figure
plot(1:nrep,[irf.f,irf.z],'-o');
legend('inflation','interest rates','Location','southoutside');
xlabel('Quarters');
ylabel('relative deviation (%)');
set(gcf, 'Color', [1,1,1]);
hline = refline([0 0]);
set(hline,'LineStyle','--','color','black');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Stochastic Simulation            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsimul = 100; % Number of simulations

nlong=100; % Number of periods in each simulations

% --------------------------------------------
% Technology and money growth process

% Drawing innovations
aleaa = normrnd(0,epsilon_a,nlong,nsimul);
aleag = normrnd(0,epsilon_g,nlong,nsimul);

% Drawing autoregressive processes
epsa(1,:) = aleaa(1,:);
epsg(1,:) = aleag(1,:);

for t=2:nlong
    epsa(t,:) = rho_a*epsa(t-1,:) + aleaa(t);
    epsg(t,:) = rho_g*epsg(t-1,:) + aleag(t);
end

% --------------------------------------------
% Simulation of the model

h = waitbar(0,'Simulations are running, please wait...');
for j=1:nsimul
    waitbar(j / nsimul)
    
    S(:,1)=[0;0]; % State variables initialization
    for i=2:nlong;
        S(:,i)=M(1:2,1:2)*S(:,i-1)+M(1:2,3)*epsa(i-1,j)+M(1:2,4)*epsg(i-1,j); % We compute the dynamics of the state variables
    end;
    
    for i=1:nlong;
        CC(i,j) = PID(1,1:2)*S(:,i)+PID(1,3)*epsa(i,j)+PID(1,4)*epsg(i,j);
        HC(i,j) = PID(2,1:2)*S(:,i)+PID(2,3)*epsa(i,j)+PID(2,4)*epsg(i,j);
        WC(i,j) = PID(3,1:2)*S(:,i)+PID(3,3)*epsa(i,j)+PID(3,4)*epsg(i,j);
        RC(i,j) = PID(4,1:2)*S(:,i)+PID(4,3)*epsa(i,j)+PID(4,4)*epsg(i,j);
        YC(i,j) = PID(5,1:2)*S(:,i)+PID(5,3)*epsa(i,j)+PID(5,4)*epsg(i,j);
        PC(i,j) = PID(6,1:2)*S(:,i)+PID(6,3)*epsa(i,j)+PID(6,4)*epsg(i,j);
        FC(i,j) = PID(7,1:2)*S(:,i)+PID(7,3)*epsa(i,j)+PID(7,4)*epsg(i,j);
        IC(i,j) = (YC(i,j)-CY*CC(i,j))/CY;
    end
end
close(h)


%%%%%% Moments computation %%%%%%%
% HP filter
[~,CCe] = hpfilter(CC,1600);
[~,HCe] = hpfilter(HC,1600);
[~,WCe] = hpfilter(WC,1600);
[~,RCe] = hpfilter(RC,1600);
[~,YCe] = hpfilter(YC,1600);
[~,PCe] = hpfilter(PC,1600);
[~,FCe] = hpfilter(FC,1600);
[~,ICe] = hpfilter(IC,1600);

% Standard deviation
CCsd = mean(std(CCe));
HCsd = mean(std(HCe));
WCsd = mean(std(WCe));
RCsd = mean(std(RCe));
YCsd = mean(std(YCe));
PCsd = mean(std(PCe));
FCsd = mean(std(FCe));
ICsd = mean(std(ICe));


%% Plot 1 of the simulations

figure
subplot(221),plot(1:nlong,[CC(:,1),YC(:,1)]);
legend('Consumption','Output','Location','southoutside');
title('Revenu and consumption')
xlabel('Quarters')
ylabel('% Dev.')
set(gcf, 'Color', [1,1,1]);

subplot(222),plot(YC(:,1));
title('Output')
xlabel('Quarters')
ylabel('% Dev.')
set(gcf, 'Color', [1,1,1]);

subplot(223),plot(WC(:,1));
title('Wages')
xlabel('Quarters')
ylabel('% Dev.')
set(gcf, 'Color', [1,1,1]);

subplot(224),plot(HC(:,1));
title('Hours')
xlabel('Quarters')
ylabel('% Dev.')
set(gcf, 'Color', [1,1,1]);

figure
subplot(221),plot(aleaa(:,1));
title('White noise, Technology')
xlabel('Quarters')
ylabel('% Dev.')
set(gcf, 'Color', [1,1,1]);

subplot(222),plot(aleag(:,1));
title('White noise, Money growth')
xlabel('Quarters')
ylabel('% Dev.')
set(gcf, 'Color', [1,1,1]);

subplot(223),plot(epsa(:,1));
title('Total factor productivity')
xlabel('Quarters')
ylabel('% Dev.')
set(gcf, 'Color', [1,1,1]);

subplot(224),plot(epsg(:,1));
title('Money growth rate')
xlabel('Quarters')
ylabel('% Dev.')
set(gcf, 'Color', [1,1,1]);

figure
plot(1:nlong,[FC(:,1), RC(:,1), epsg(:,1)]);
legend('inflation','interest rate','money growth rate','Location','southoutside');
xlabel('Quarters');
ylabel('relative deviation (%)');
set(gcf, 'Color', [1,1,1]);