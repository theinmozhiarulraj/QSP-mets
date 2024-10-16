% script to calculate tumor growth rate in the simulations        

organs=["primary","lung","others"];
suffix=["","_Ln1","_other"];

% set the target initial and final tumor volumes for growth rate
% calculation

tumor_target_initial = 0.2; % 0.2 cm or 2 mm
tumor_target_final = 2; % 2 cm or 20 mm

growth_rate_temp=NaN(1,length(organs));

for i=1:length(organs)

    % get the tumor volume (in ml)
    pos = find(ismember(simData.DataNames,append("V_T",suffix(i))));
    tumor_volume = simData.Data(:,pos);
    
    % calculate tumor diameter from volume (in cm)
    
    tumor_dia = 2*(3/(4*pi)*tumor_volume).^(1/3);
    time=simData.time; % in days
    
    % index where tumor diameter reaches target values
    idx_initial = find((tumor_target_initial-tumor_dia)<0,1); 
    idx_final = find((tumor_target_final-tumor_dia)<0,1); 
    
    % calculate the growth rate
    growth_rate_temp(i)=(log(tumor_dia(idx_final))-log(tumor_dia(idx_initial)))/(time(idx_final)-time(idx_initial));

end
