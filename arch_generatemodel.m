%thein: modified script for model generation in ARCH
clear
close all
sbioreset

%% Create the model
immune_oncology_model_TNBC

% number of tumors
ntumors=length(TumorArray);
disp(ntumors)

%% Define dosing
dose_schedule = schedule_dosing({'pembrolizumab'});

% total number of virtual patients
n_PSA = 1000; 

% number of VPs per job array
npatients_subset=50;

% suffix for tumors. this is only used for growth rate calculation
organs=["primary","lung","others"];
suffix=["","_Ln1","_other"];

% variables to store growth rates
growth_rates=NaN(n_PSA,length(organs));

% thein: 1=read parameters from already generated VPs;
% 0=randomly sample VP parameters again
readfromfile = 0;

% Set distributions of the selected parameters for random sampling
% change parameter file name here
params_in  = PSA_param_in_TNBC_w60_vct(TumorArray,cancer_clones,n_T_specs);
% Set boundaries for model species (only include those added in the present model)
params_out = PSA_param_out(model);
% Add selected model outputs to sensitivity analysis
params_in  = PSA_param_obs(params_in);
% Randomly generate parameter sets using Latin-Hypercube Sampling
params_in  = PSA_setup(model,params_in,n_PSA,readfromfile);

% save the workspace
save model_workspace_w60.mat