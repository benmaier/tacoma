clear;
close all;

% set random_seed
random_seed = 12634;

% define Flockwork parameters
edges = [];
N = 500;
Q = 0.5;
k = 1/(1-Q);
use_random_rewiring = 0; %Careful! "false" doesn't work.

% define epidemic parameters
t_run_total = 123.56;
R0 = 4;
recovery_rate = 1;
susceptible_rate = 1;
rewiring_rate = 1;
number_of_vaccinated = 10;
number_of_initially_infected = 10;

infection_rate = R0 * recovery_rate / k;

%% ========================= Equilibration ======================

num_timesteps = 0; % if this is zero, the function determines this number
                   % by itself, using N and Q
edges = FlockworkEq(edges,N,Q,random_seed,num_timesteps);
random_seed = random_seed + 1;

%% ========================= Simulation ======================

num_timesteps = N; 
edges = FlockworkEq(edges,N,Q,random_seed,num_timesteps);
random_seed = random_seed + 1;

%% =========================== SIS ======================
% start this simulation as a fully disconnected graph and equilibrate
% within simulation
equilibrate_flockwork = 1; % true doesn't work

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
                                  equilibrate_flockwork,...
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
ylim([0,1])

figure;

%% ========================== SIRS ============================
% start this simulation with the equilibrated graph from above
equilibrate_flockwork = 0; % false doesn't work

random_seed = random_seed + 1;

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
                                  equilibrate_flockwork,...
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
ylim([0,1])

figure;

%% ===================== SIR ==================
random_seed = random_seed + 1;

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
                                  equilibrate_flockwork,...
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
ylim([0,1])
