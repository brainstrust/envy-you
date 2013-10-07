clean
tic

% global variables
global CASE J_PLC startpulse lengthpulse C_Hillmann stretch_ch only_Koenig

%% Parameters to adjust the model:
t_start = 0;
t_end = 500;
startpulse  = 200;  % (s) 
lengthpulse = 200;  % (s) 
CASE        = 2;    % (see all_constants.m for details)
J_PLC 		= 0.18;  % (muM s-1) EC agonist concentration  
C_Hillmann  = 1;    % scaling factor for the Hai&Murphy rate constants (see all_constants.m for details)
stretch_ch  = 'ON'; % choose 'ON'/'OFF' to activate/deactivate stretch-activated channels in EC and SMC
only_Koenig = 'OFF';% choose 'ON'/'OFF' to simulate only the Koenigsberger model (other sub-models will still be considered, but the KIR channel is set to 0)

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
hold all

%% save figures & parameters
save_all()


% to create .tikz figures:
% matlab2tikz('test.tikz', 'height', '\figureheight', 'width', '\figurewidth');


% figure; plot(time,DATA(:,smcoff+flu.M)+state(:,ind.AMp)+state(:,ind.AM)+state(:,ind.Mp)); hold on;
% plot(time,DATA(:,smcoff+flu.M),'r'); plot(time,state(:,ind.Mp),'g'); plot(time,state(:,ind.AMp),'b');plot(time,state(:,ind.AM),'k');
% legend('Total Myosin','[M]','[Mp]','[AMp]','[AM]')
% title('New')

% to plot a single flux, type in plot(time,DATA(:,flu.(name))     
% to plot a single state variable, type in plot(time,DATA(:,ind.(name))
%(don't forget to put the offset!! e.g. smcoff+flu.1_c)
