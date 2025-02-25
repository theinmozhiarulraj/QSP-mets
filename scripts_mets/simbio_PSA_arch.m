% Function to run the model with parameter snesitivity analysis cases
%
% Inputs: model         -- simbio model object
%         params_in     -- object containing LHS values of model parameters
%         params_out    -- object containing a list of output parameters
%         dose_schedule -- dose object
%
% Outputs: simDataPSA  -- results from the LHS simulations
%          params_out  -- updated out parameters object containing patients
%                         status (1:patient,0:healthy,-1:simulation failed to converge)


function [simDataPSA, params_out, growth_rates] = simbio_PSA_arch(model,params_in,params_out,TumorArray,growth_rates,suffix,start_patient,end_patient,dose_schedule)
config = getconfigset(model);
time = get(config.SolverOptions,'OutputTimes');

n_dose = 1;
idx_nab = 0;
for k = 1:length(dose_schedule)
    if ~isempty(strfind(dose_schedule(k).Name, 'nabp'))
        idx_nab(n_dose) = k;
        n_dose = n_dose + 1;
        dose = dose_schedule(k).Amount / get_params(model,'BSA'); % default BSA
    end
end

% Check if ICs are set
new_sim = ~isfield(params_out,'ICs');

n_PSA = length(params_in.(params_in.names{1}).LHS);

for i = start_patient:end_patient
    display(['Sample ',num2str(i),'/',num2str(n_PSA)]);
    % Set the new parameters
    if i > 1
      %delete(model_PSA)
    end
    model_PSA = copyobj(model);
    variantObj = addvariant(model_PSA, ['v',num2str(i,'%5.5i')]);
    for j = 1:length(params_in.names)
        if ~isempty(sbioselect (model, 'Type', 'parameter', 'Name', params_in.names{j}))
            addcontent(variantObj, {'parameter', params_in.names{j}, 'Value', params_in.(params_in.names{j}).LHS(i)});
        elseif ~isempty(sbioselect (model, 'Type', 'compartment', 'Name', params_in.names{j}))
            addcontent(variantObj, {'compartment', params_in.names{j}, 'Capacity', params_in.(params_in.names{j}).LHS(i)});
        else
            disp(['Unable to identify parameter/compartment named ', params_in.names{j}, ' in the model'])
        end
    end

    %thein: add fold changes here
    %addcontent(variantObj, {'parameter','k_C1_growth_Ln1', 'Value', params_in.('foldchange_k_C1_growth_Ln1').LHS(i)*params_in.('k_C1_growth').LHS(i)});

    % Set Initial Conditions
    if (new_sim || isempty(params_out.ICs(i).Values))
      %addcontent(variantObj,{'parameter','seeding','Value',0}) %thein
      [model_PSA,growth_rates,success,~] = initial_conditions(model_PSA,TumorArray,i,growth_rates,suffix,'Variant',variantObj);
      params_out.patient(i) = double(success);
    else
      success = logical(params_out.patient(i));
      model_PSA = set_ICs(model_PSA,params_out.ICs(i).Values);
    end

    % Run the model with drugs
    if (success)
          if idx_nab ~= 0
              for k = 1:length(idx_nab)
                  dose_schedule(idx_nab(k)).Amount = dose*params_in.BSA.LHS(i);
              end
          end

          %thein
          %addcontent(variantObj,{'parameter','do_surgery','Value',1})
          addcontent(variantObj,{'parameter','seeding','Value',0})
          disp("seeding was switched off- thein!!!!!!!");

          try
              simData = sbiosimulate(model_PSA,[],variantObj,dose_schedule);
              % Remove simulations that reached the time limit
              if size(simData.Data,1) < length(time)
                  simData = [];
                  params_out.patient(i) = -1;
                  disp('Simulation takes longer than the preset time limit');
              end
          catch
              disp('Integration tolerance not met');
              simData = [];
              params_out.patient(i) = -1;
          end

    else
        disp('Initial conditions not reached');
        simData = [];
    end

    % save model output struture within an array of structures
    simDataPSA(i).simData = simData;

    % save model ICs
    params_out.ICs(i).Values = get_ICs(simData);
end

params_out.iPatient = find(params_out.patient == 1);
params_out.iHealthy = find(params_out.patient == 0);

%% parfor implementation gives error on MJ's mac machine
% If had "Not enough input arguments" error from parfor
% Enter matlab command window 'prefdir' -->
% You can see C:\Users\xxx\AppData\Roaming\MathWorks\MATLAB\R2016a (or wahtever you see on your machine) -->
% Delete local_cluster_jobs and restart matlab -->
% I hope this will fix the issue.