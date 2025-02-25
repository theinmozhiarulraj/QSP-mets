% Function to add RECIST and response status to the postprocessed parameters
%
% Inputs: simDataPSA      -- Object containing the simbiology model outputs
%                            for all batch simulations
%         simDataPSApost  -- Object containing the postprocessed simbiology
%                           model outputs for all batch simulations
%         params_outPlaus -- Object containing organized outputs
%
% Outputs: simDataPSAout  -- Updated Object containing organized outputs
%                            with response status and RECIST added


function params_outUpdate = responseStatus_mets(simDataPSA,simDataPSApost,params_out)

n_PSA = length(params_out.patient);
params_outUpdate = params_out;

% thein: testing frequency in weeks
test_freq=9;
duration_sd=24; % weeks
duration_sd_days=duration_sd*7;

time = simDataPSA(params_out.iPatient(1)).simData.time;
time=time(1:test_freq*7:end);

% find the index of postprocessed percent change in tumor size
k = find(strcmp(simDataPSApost(params_out.iPatient(1)).simData.DataNames, 'D_T_all_perc' ));
h = find(strcmp(simDataPSApost(params_out.iPatient(1)).simData.DataNames, 'D_T_all' ));
for j = 1:n_PSA
    % see if the initial tumor diameter has reached
    if (params_out.patient(j)==1)
        % Find the final perentage change in the diameter
        D_T_perc = simDataPSApost(j).simData.Data(:,k);
        D_T = simDataPSApost(j).simData.Data(:,h);

        % thein: testing at intervals
        D_T_perc=D_T_perc(1:test_freq*7:end);
        D_T=D_T(1:test_freq*7:end);

        idx_min = find(D_T_perc==min(D_T_perc), 1);

        % Add RECIST status, TTR, and DOR of the patient
        if (min(D_T_perc) <= -30)
            idx_res = find(D_T_perc<=-30, 1);
            params_outUpdate.RECIST_mets(j,1) = {'CR/PR'};

            if(min(D_T)<0.2)
                params_outUpdate.RECIST_mets_ref(j,1) = {'CR'};
            else
                params_outUpdate.RECIST_mets_ref(j,1) = {'PR'};
            end

            params_outUpdate.Response_mets(j,1) = {'Responder'};
            params_outUpdate.TTR_mets(j,1) = time(idx_res) - time(1);
            if (max(D_T(idx_min:end)) >= 1.2*min(D_T(idx_min))) && (max(D_T(idx_min:end)) - min(D_T) >= 0.5)
                params_outUpdate.DOR_mets(j,1) = time(find(D_T(idx_min:end) >= 1.2*min(D_T(idx_min)), 1)+idx_min-1)-time(idx_res);
            else
                params_outUpdate.DOR_mets(j,1) = time(end)-time(idx_res);
            end
        elseif (max(D_T(idx_min:end)) >= 1.2*min(D_T(idx_min))) && (max(D_T(idx_min:end)) - min(D_T) >= 0.5) ...
            && (time(find(D_T(idx_min:end) >= 1.2*min(D_T(idx_min)), 1)+idx_min-1) - time(1) <= duration_sd_days)
            % assuming a minimum duration of stable disease of 8 weeks (56 days)
            params_outUpdate.RECIST_mets(j,1) = {'PD'};
            params_outUpdate.RECIST_mets_ref(j,1) = {'PD'};
            params_outUpdate.Response_mets(j,1) = {'Non-responder'};
            params_outUpdate.TTR_mets(j,1) = 0; % thein changed this from time(end) -> 0
            params_outUpdate.DOR_mets(j,1) = 0;
        else
            params_outUpdate.RECIST_mets(j,1) = {'SD'};
            params_outUpdate.RECIST_mets_ref(j,1) = {'SD'};
            params_outUpdate.Response_mets(j,1) = {'Non-responder'};
            params_outUpdate.TTR_mets(j,1) = 0; % thein changed this from time(end) -> 0
            params_outUpdate.DOR_mets(j,1) = 0;
        end

        % If the patient progressed add time of the event for Kaplan-Meier plots (time to progression)
        % PMID:29765247
        D_T_perc_max = max(D_T_perc);
        if (D_T_perc_max > (20))
            idx_max = find(D_T_perc>=20, 1);
            params_outUpdate.tPFS_mets(j,1) = time(idx_max);
        else
            params_outUpdate.tPFS_mets(j,1) = time(end);
        end

        % Assuming efficacy endpoints are evaluated every 8 weeks
        % params_outUpdate.TTR(j,1) = 56*ceil(params_outUpdate.TTR(j,1)/56);
        % params_outUpdate.DOR(j,1) = 56*ceil(params_outUpdate.DOR(j,1)/56);
        % params_outUpdate.tPFS(j,1) = 56*ceil(params_outUpdate.tPFS(j,1)/56);

    else
        params_outUpdate.RECIST_mets(j,1)   = {'NP'};
        params_outUpdate.RECIST_mets_ref(j,1) = {'NP'};
        params_outUpdate.Response_mets(j,1) = {'Non-patient'};
    end
end
