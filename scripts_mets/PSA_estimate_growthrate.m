function growth_rates = PSA_estimate_growthrate(simData,patient_id,growth_rates,idx,suffix)
 
% script to calculate tumor growth rate in the simulations        

% organs=["primary","lung","others"];
% suffix=["","_Ln1","_other"];

% set the target initial and final tumor volumes for growth rate
% calculation

idx_final=idx;
tumor_target_initial = 0.2; % 0.2 cm or 2 mm

    for i=1:length(suffix)
    
        % get the tumor volume (in ml)
        pos = find(ismember(simData.DataNames,append("V_T",suffix(i))));
        tumor_volume = simData.Data(:,pos);
        
        % calculate tumor diameter from volume (in cm) 
        tumor_dia = 2*(3/(4*pi)*tumor_volume).^(1/3);
        time=simData.time; % in days
        
        if(tumor_dia(end)>0.2) % 0.2 cm

            % index where tumor diameter reaches target values
            idx_initial = find((tumor_target_initial-tumor_dia)<0,1); 
            
            % calculate the growth rate
            growth_rates(patient_id,i)=(log(tumor_volume(idx_final))-log(tumor_volume(idx_initial)))/((time(idx_final)-time(idx_initial))*24);

        end

    end

end