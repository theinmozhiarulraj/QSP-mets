% function to run VP simulations in ARCH
function [simDataPSA, params_out, growth_rates] = arch_simulate(taskno)
    
    % load the model 
    load('model_workspace_w60.mat')

    %npatients_subset=3; % comment later
    
    start_patient=npatients_subset*(taskno-1)+1;
    end_patient=start_patient+npatients_subset-1;
    
    sbioaccelerate(model, dose_schedule)
    
    tic
    [simDataPSA, params_out, growth_rates] = simbio_PSA_arch(model,params_in,params_out,TumorArray,growth_rates,suffix,start_patient,end_patient,dose_schedule);
    toc
    
    stringname=append('PSA_output',num2str(taskno),'.mat');
    save(stringname)

end
