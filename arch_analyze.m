%thein: modified script to analyze results from ARCH
% some functions might not work in the QSP model with mets but are not
% required

clear
close all
sbioreset

load('model_workspace_w60.mat')

%npatients_subset=3; % comment later

growth_rates=NaN(n_PSA,length(organs));
%simDataPSA=struct;

% number of tasks
ntasks=n_PSA/npatients_subset;

for i=1:ntasks

    % concatenate output of all tasks

    disp(i)

    tic

    loaded_ws=load(['PSA_output' num2str(i) '.mat']);

    start_patient=npatients_subset*(i-1)+1;
    end_patient=start_patient+npatients_subset-1;

    growth_rates(start_patient:end_patient,:)=loaded_ws.growth_rates(start_patient:end_patient,:);
    % growth_rates=[growth_rates,growth_rates_temp];

    % put together simDataPSA

%     if(i==1)
%         simDataPSA.simData=load(['PSA_output' num2str(i) '.mat']).simDataPSA.simData;
%     else
%         simDataPSA.simData=[load(['PSA_output' num2str(i) '.mat']).simDataPSA.simData];
% 
%     end

    for k=start_patient:end_patient

        simDataPSA(k).simData=loaded_ws.simDataPSA(k).simData;

    end

    % put together params.out

    params_out.patient(start_patient:end_patient)=loaded_ws.params_out.patient(start_patient:end_patient);
     
    params_out.ICs(1,start_patient:end_patient)=loaded_ws.params_out.ICs(1,start_patient:end_patient);

%     if(i==1)
% 
%         params_out.iPatient(start_patient:end_patient)=load(['PSA_output' num2str(i) '.mat']).params_out.iPatient(start_patient:end_patient);
%     else
% 
%         params_out.iPatient(start_patient:end_patient)=load(['PSA_output' num2str(i) '.mat']).params_out.iPatient(start_patient:end_patient);
% 
%     end

tic

end

ncount_patient=1;
ncount_healthy=1;
params_out.iPatient=[];
params_out.iHealthy=[];

for j=1:n_PSA
    
    if(params_out.patient(j)==1) % patient

        params_out.iPatient(ncount_patient)=j;
        ncount_patient=ncount_patient+1;

    else % healthy if 0

        params_out.iHealthy(ncount_healthy)=j;
        ncount_healthy=ncount_healthy+1;

    end

end


%% Postprocess
growth_rates = array2table(growth_rates,'VariableNames',organs);

% thein: write growth rates in a file
writetable(growth_rates,'growth_rates.csv');

% Postprocess Data -> Calculate Clonality, Percentages and ...
simDataPSApost = PSA_post(simDataPSA,params_in,params_out,TumorArray);

% Add pre-treatment observables to the params_in
params_in = PSA_preObs(simDataPSA,simDataPSApost,params_in,params_out);

% Prepare the data for the rest of the analysis
params_out = PSA_prep(simDataPSA,simDataPSApost,params_out);

% Save and print data of interest (by assigning a unique code name for the trial)
%sprint_data(simDataPSA, simDataPSApost, params_in, params_out, 'test_thein')
%sprint_data(simDataPSA, simDataPSApost, params_in, params_out, 'thein') 

%% Perform and plot different types of analysis

% Partial Rank Correlation Coefficients
% PSA_PRCC(params_in,params_out,'plausible')
% PSA_PRCC(params_in,params_out)

% t-SNE Analysis
% PSA_tSNE(params_in,params_out,'plausible')
% PSA_tSNE(params_in,params_out,'patient')

% eFAST

% Principle Component Analysis

%% Plot Results
% close all
% Plot percent change in size and RECIST
% PSA_plot_RECIST(simDataPSA,simDataPSApost,params_out)

% Tumor Size
% PSA_plot_TumSize(simDataPSA,simDataPSApost,params_out)

% Kaplan-Meier Progression-free survival (PFS)
% PSA_plot_KaplanMeier(simDataPSA,simDataPSApost,params_out)

%% Waterfall plots
% Plot percent change in size using waterfall plots for a parameter.
% Waterfall plots can be plotted using either end tumor sizes or best overall responses.
%
% PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'n_T1_clones')
% PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'k_P1_d1')
 PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'k_C1_growth')
 PSA_plot_Waterfall_mets(simDataPSApost,model,params_in,params_out,'k_C1_growth',TumorArray)

 save('analyzed_results_w60.mat')
