clear;
close all;

% set random_seed
random_seed = 123;

% define Flockwork parameters
edges = [];
N = 100;
Q = 0.5;
k = 1/(1-Q);
use_random_rewiring = false;

% define epidemic parameters
t_run_total = 1234.56;
R0 = 1.5;
recovery_rate = 1;
susceptible_rate = 1;
rewiring_rate = 1;
number_of_vaccinated = 10;
number_of_initially_infected = N/2;

infection_rate = R0 * recovery_rate / k;

% =========================== SIS ======================

[I,SI,R0,edgelist] = FlockworkSIS(edges,...
                                  N,...
                                  Q,...
                                  t_run_total,...
                                  infection_rate,...
                                  recovery_rate,...
                                  rewiring_rate,...
                                  number_of_vaccinated,...
                                  number_of_initially_infected,...
                                  use_random_rewiring,...
                                  random_seed...
                                 );

figure;

t = I(:,1);
i = I(:,2) / N;
plot(t,1-i); hold on;
plot(t,i); 

title('SIS')
legend('S','I')
xlabel('time t')
ylabel('population ratios')

figure;

%% ========================== SIRS ============================

[I,R,SI,R0,edgelist] = FlockworkSIRS(edges,...
                                  N,...
                                  Q,...
                                  t_run_total,...
                                  infection_rate,...
                                  recovery_rate,...
                                  susceptible_rate,...
                                  rewiring_rate,...
                                  number_of_vaccinated,...
                                  number_of_initially_infected,...
                                  use_random_rewiring,...
                                  random_seed...
                                 );
                             
t = R(:,1);
r = R(:,2) / N; 
plot(t,r); hold on;

t = I(:,1);
i = I(:,2) / N;
plot(t,i); 

title('SIRS')
legend('R','I')
xlabel('time t')
ylabel('population ratios')

figure;

%% ===================== SIR ==================

% let SIR run until all infected are gone
t_run_total = 0;

[I,R,SI,R0,edgelist] = FlockworkSIR(edges,...
                                  N,...
                                  Q,...
                                  infection_rate,...
                                  recovery_rate,...
                                  rewiring_rate,...
                                  t_run_total,...
                                  number_of_vaccinated,...
                                  number_of_initially_infected,...
                                  use_random_rewiring,...
                                  random_seed...
                                 );
                             
t = R(:,1);
r = R(:,2) / N; 
plot(t,r); hold on;

t = I(:,1);
i = I(:,2) / N;
plot(t,i); 


title('SIR')
legend('R','I')
xlabel('time t')
ylabel('population ratios')
