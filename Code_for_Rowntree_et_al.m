% The following code was used in the paper "Bilateral feedback in 
%oscillator model is required to explain the coupling dynamics of Hes1 with the cell cycle."
% By Rowntree et. al. 
% Submitted to MDPI mathematics

clc, clear all, close all

%% Load in data

load Hes1_expression_data_1
load Hes1_expression_data_2
load Hes1_expression_data_3
load Phase_reconstruction_of_data_1
load Phase_reconstruction_of_data_2
load Phase_reconstruction_of_data_3

% Vectors of normalised Hes1 expression in pseudo-time
Hes1_expression_data = {Hes1_expression_data_1, Hes1_expression_data_2, Hes1_expression_data_3};

% Vectors of phase reconstruction of Hes1 expression, values between 0 and 2pi.
Phase_reconstruction_of_data = {Phase_reconstruction_of_data_1,Phase_reconstruction_of_data_2,Phase_reconstruction_of_data_3};

%% Parameters and preallocation

n=2;   % Number of particles.

% Taken from mean Hes1 periodicity found in Sabherwal et. al. (2001)

Hes1_Oscillator_intrinsic_periodicity = ; 
Cell_cycle_Oscillator_intrinsic_periodicity = ; 

% Set up vector with these two intrinsic frequencies denoted as $omega_{Hes1}$,
% $omega_{C.C.}$ in paper

Omega=[1/Hes1_Oscillator_intrinsic_periodicity;...
       1/Cell_cycle_Oscillator_intrinsic_periodicity]; 
   
% Step Size for numerical scheme
h=0.1; 

% Max number of interations to run numerical scheme
iter = 10000;

% Construct time vector
t = 0:h:h*iter; 

% Set initial conditions
Initial_conditions_1 = Phase_reconstruction_of_data{1}(1);
Initial_conditions_2 = Phase_reconstruction_of_data{2}(1);
%.....
Initial_conditions_n = Phase_reconstruction_of_data{3}(1);

Initial_conditions = [Initial_conditions_1,Initial_conditions_2,Initial_conditions_n]; %
                 
% Set up array for the solution vectors, one cell for each intial condition
% and each cell will contain Hes1 and cell cycle phase solutions from the
% model simulations.

theta = cell(1,length(Initial_conditions));

% Set up array for the Residual error_matrices

Residual_error_matix = cell(1,length(Initial_conditions));
Residual_error_matix_transpose = cell(1,length(Initial_conditions));

% Set up array for the best coupling strength parameters

Best_K_position = cell(1,length(Initial_conditions));
Best_K_parameter_value = cell(1,length(Initial_conditions));

% Set up grid mesh for coupling parameters to search through
kappa_scale = 0.01:0.01:1;

%% Run the model

% Loop through initial conditions

for this_initial_condition = 1:length(Initial_conditions)

    % Set up matrix grid to contain residual error
    Residual_error_matix{this_initial_condition} = zeros(length(kappa_scale),length(kappa_scale)); 
    Residual_error_matix_transpose{this_initial_condition} = zeros(length(kappa_scale),length(kappa_scale)); 
    
    % Loop through two dimensional parameter space
    
    for  kappa_index_cell_cycle_onto_Hes1 = 1:length(kappa_scale)
        for  kappa_index_Hes1_onto_cell_cycle = 1:length(kappa_scale)
            
            % Set up vector of K (omega*kappa) values to be fed into the
            % function
            K=[kappa_scale(kappa_index_cell_cycle_onto_Hes1)*Omega(1);...
               kappa_scale(kappa_index_Hes1_onto_cell_cycle)*Omega(2)];

            % Set up solution vectors as zeros except first value which is
            % our initial conditions for Hes1 at 0 for cell cycle
            theta{this_initial_condition}=zeros(n,iter); % Phase vectors
            theta{this_initial_condition}(:,1)=[Initial_conditions(this_initial_condition);0]; % Initial values
     
            % Run the numerical scheme (4th order Runge Kutta)
            for j=1:iter
                k1=Coupled_oscillator_model_function(theta{this_initial_condition}(:,j),K,n,Omega);
                k2=Coupled_oscillator_model_function(theta{this_initial_condition}(:,j)+0.5*h*k1,K,n,Omega);
                k3=Coupled_oscillator_model_function(theta{this_initial_condition}(:,j)+0.5*h*k2,K,n,Omega);           
                k4=Coupled_oscillator_model_function(theta{this_initial_condition}(:,j)+h*k3,K,n,Omega);
                theta{this_initial_condition}(:,j+1)=theta{this_initial_condition}(:,j)+(h/6)*(k1+2*k2+2*k3+k4);
                
                % Halt the simulation once the cell cycle oscillator
                % surpasses 2pi
                
                if theta{this_initial_condition}(2,j)>2*pi
                    break
                end

            end
            
            % Truncate the solution vector to remove zeros after simulation
            % ended
            theta{this_initial_condition} = theta{this_initial_condition}(:,1:j-1);
            
            % In order to compare Hes1 phase dynamics, we make the solution
            % vector the same length as our data vector which is in
            % pseudo-time
            Hes1_phase_solution_interpolation =interp1(1:length(theta{this_initial_condition}(1,:)),...
                                                            theta{this_initial_condition}(1,:)',...
                                                            linspace(1,length(theta{this_initial_condition}),...
                                                            length(Phase_reconstruction_of_data{this_initial_condition}(1))));

            Residual_error_matix{this_initial_condition}(kappa_index_cell_cycle_onto_Hes1,kappa_index_Hes1_onto_cell_cycle)=...
                sqrt(sum((Phase_reconstruction_of_data{this_initial_condition}(1)-...
                          Hes1_phase_solution_interpolation).^2));

        end
    end

    where_is_the_min_1 = find(Residual_error_matix{this_initial_condition}...
        ==min(min(Residual_error_matix{this_initial_condition})));

    if numel(where_is_the_min_1)>1
        where_is_the_min_1=where_is_the_min_1(1);
    end

    Best_K_position{this_initial_condition}(1) = mod(where_is_the_min_1,length(kappa_scale));

    Residual_error_matix_transpose{this_initial_condition}= Residual_error_matix{this_initial_condition}';

    where_is_the_min_2 = find(Residual_error_matix_transpose{this_initial_condition}...
        ==min(min(Residual_error_matix_transpose{this_initial_condition})));

    if numel(where_is_the_min_2)>1
        where_is_the_min_2=where_is_the_min_2(1);
    end

    Best_K_position{this_initial_condition}(2) = mod(where_is_the_min_2,length(kappa_scale));

    if Best_K_position{this_initial_condition}(1)==0
        Best_K_position{this_initial_condition}(1)=11;
    end
    if Best_K_position{this_initial_condition}(2)==0
        Best_K_position{this_initial_condition}(2)=11;
    end

    Best_K_parameter_value{this_initial_condition} = kappa_scale(Best_K_position{this_initial_condition});

end

%% Run model one more time but with the best parameters to save 
%the cell cycle length simulations and the optimal solution vectors


CC_duration = zeros(1,length(Initial_conditions));
Hes1_phase_solution_interpolation = cell(1,length(Initial_conditions));
Hes1_dynamics_rebuilt = cell(1,length(Initial_conditions));
Cell_cycle_phase = cell(1,length(Initial_conditions));

for this_initial_condition = 1:length(Initial_conditions)
    
    K=Best_K_parameter_value{this_initial_condition}*max([Omega(1),Omega(2)]);
    K=K';

    theta{this_initial_condition}=zeros(n,iter); % Phase vectors
    theta{this_initial_condition}(:,1)=[Initial_conditions(this_initial_condition);0]; % Initial values
     
    for j=1:iter
        
        k1=Coupled_oscillator_model_function(theta{this_initial_condition}(:,j),K,n,Omega);
        k2=Coupled_oscillator_model_function(theta{this_initial_condition}(:,j)+0.5*h*k1,K,n,Omega);
        k3=Coupled_oscillator_model_function(theta{this_initial_condition}(:,j)+0.5*h*k2,K,n,Omega);           %4-th order Runge-Kutta method.
        k4=Coupled_oscillator_model_function(theta{this_initial_condition}(:,j)+h*k3,K,n,Omega);
        theta{this_initial_condition}(:,j+1)=theta{this_initial_condition}(:,j)+(h/6)*(k1+2*k2+2*k3+k4);
        
        if theta{this_initial_condition}(2,j)>2*pi
            break
        end
    end
        
    % Truncate the solution vector to remove zeros after simulation
    % ended
    theta{this_initial_condition} = theta{this_initial_condition}(:,1:j-1);
    
    % In order to compare Hes1 phase dynamics, we make the solution
    % vector the same length as our data vector which is in
    % pseudo-time
    Hes1_phase_solution_interpolation{this_initial_condition} =interp1(1:length(theta{this_initial_condition}(1,:)),...
        theta{this_initial_condition}(1,:)',...
        linspace(1,length(theta{this_initial_condition}),...
        length(Phase_reconstruction_of_data{this_initial_condition})));
    
    Cell_cycle_phase{this_initial_condition} = interp1(1:length(theta{this_initial_condition}(2,:)),...
        theta{this_initial_condition}(2,:)',...
        linspace(1,length(theta{this_initial_condition}),...
        length(Phase_reconstruction_of_data{this_initial_condition})));
       
    Hes1_dynamics_rebuilt{this_initial_condition} = zscore(cos(Hes1_phase_solution_interpolation{this_initial_condition}));

    CC_duration(this_initial_condition) = j;
    
end

colour = {'r','b','g'};
colour_dot = {'r.','b.','g.'};

figure
for this_initial_condition = 1:length(Initial_conditions)
    hold on
    plot(Cell_cycle_phase{this_initial_condition},...
        mod(Hes1_phase_solution_interpolation{this_initial_condition},2*pi),...
        colour_dot{this_initial_condition},'MarkerSize',15)

ylabel('Hes1 phase')
xlabel('Cell cycle phase')
xlim([0 2*pi])
ylim([0 2*pi])
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
yticks(0:pi/2:2*pi)
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
end

for this_initial_condition = 1:3
    figure
    hold on
    plot(zscore(Hes1_expression_data{this_initial_condition}),colour{this_initial_condition},'LineWidth',2)
    plot(zscore(Hes1_dynamics_rebuilt{this_initial_condition}),'k','LineWidth',2)
    ylim([-2.5 2.5])
    xticks([])
    yticks([])
    ylabel({'Normalised Hes1 intensity'})
    xlabel('Pseudo-time')
    legend({'Mean cluster trace','Model reconstruction'})
end

figure('DefaultAxesFontSize',12)
b=bar(CC_duration);
xlabel('Cluster')
ylabel('Time steps for full cell cycle')
b.FaceColor = 'flat';
b.CData(1,:) = [1 0 0];
b.CData(2,:) = [0 0 1];
b.CData(3,:) = [0 1 0];
b.FaceAlpha = 0.5;
title('Time steps required for each model simulation')

%% Oscillator model function

function f=Coupled_oscillator_model_function(x,K,n,Omega)

f=Omega+(K).*sum(sin(x*ones(1,n)-(ones(n,1)*x')))';

end