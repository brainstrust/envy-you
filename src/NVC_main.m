clean
tic

% global variables
global CASE J_PLC startpulse lengthpulse 
CASE        = 2; 
J_PLC 		= 0.2;
startpulse  = 200; %200
lengthpulse = 50; %[s] Ostby: 30  

% Time span for which the simulation runs:
t_start = 0;
t_end = 300;

%% load the constants for the fluxes and pointers:
all_indices();
all_constants();
%% load the initial conditions of the system:
state0 = InitCond();
%% Ensure single filenames for the writing of data in other files
global csvfilename
csvfilename = 'Data_simulation.csv';
try
delete(csvfilename) % remove file, if present from older simulation.
end
%% Solve the proces from initial position tot Steady State:
options = odeset('OutputFcn',@odeprogWD,'Events',@odeabort,'Stats','on','RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1); 
[t,state] = ode15s(@DEsyst,[t_start t_end],state0,options);

%% Write output and info to file/cmd
output.info.completiontime = toc;
fprintf('ODE solution time: %.3f seconds\n', output.info.completiontime)

%% Plot statement:
plot_all()

% figure; plot(time,DATA(:,smcoff+flu.M)+state(:,ind.AMp)+state(:,ind.AM)+state(:,ind.Mp)); hold on;
% plot(time,DATA(:,smcoff+flu.M),'r'); plot(time,state(:,ind.Mp),'g'); plot(time,state(:,ind.AMp),'b');plot(time,state(:,ind.AM),'k');
% legend('Total Myosin','[M]','[Mp]','[AMp]','[AM]')
% title('New')

% to plot a single flux, type in plot(time,DATA(:,flu.(name))     
% to plot a single state variable, type in plot(time,DATA(:,ind.(name))
%(don't forget to put the offset!! e.g. smcoff+flu.1_c)