% Initial Conditions Module
%
% Generates initial conditions
%
% Inputs: model           -- pre-initialized SimBiology model object
%       TumorArray        -- Array of suffixes for tumor comp. names
%         varargin        -- variant object
%                         -- true: debug mode, false: catches the error
%
% Outputs: model   -- initialized SimBiology model
%          success -- shows if it reached the initial tumor diameter
%          simData -- initial condition data
%
% future release: name/value pair inputs, IC output time resolution

function [model_out,growth_rates,success,simData] = initial_conditions(model_in,TumorArray,patient_id,growth_rates,suffix,varargin)

% number of tumors
ntumors=length(TumorArray);

% Optional Inputs
p = inputParser;
addParameter(p,'Variant','');
addParameter(p,'Debug','');
parse(p,varargin{:});
Variant = p.Results.Variant;
Debug = p.Results.Debug;
if (isempty(Debug))
    Debug = false;
end

model = copyobj(model_in);

if (~isempty(Variant))
    variantWithInitialTumourDiameter = false;
    % Check if 'Variant Object' Contains Initial Tumour Diameter
    for i = 1:length(Variant.content)
        variantParam = Variant.content{i};
        if strcmp(variantParam(2),'initial_tumour_diameter')
            % Get Initial Tumour Diameter from 'variant'
            D_Ti = sbioselect(model,'Name','initial_tumour_diameter');
            tumour_diameter.Value = variantParam{4};
            tumour_diameter.Units = D_Ti.ValueUnits;
            tumour_diameter.Notes = D_Ti.Notes;

            variantWithInitialTumourDiameter = true;
        end
    end
    if ~variantWithInitialTumourDiameter
        % Get Initial Tumour Diameter from 'model'
        D_Ti = sbioselect(model,'Name','initial_tumour_diameter');
        tumour_diameter.Value = D_Ti.Value;
        tumour_diameter.Units = D_Ti.ValueUnits;
        tumour_diameter.Notes = D_Ti.Notes;
    end
else
    % Get Initial Tumour Diameter from 'model'
    D_Ti = sbioselect(model,'Name','initial_tumour_diameter');
    tumour_diameter.Value = D_Ti.Value;
    tumour_diameter.Units = D_Ti.ValueUnits;
    tumour_diameter.Notes = D_Ti.Notes;
end

if(ntumors>1)
    % thein: met diameter
    if (~isempty(Variant))
        variantWithInitialMetDiameter = false;
        % Check if 'Variant Object' Contains Initial Met Diameter
        for i = 1:length(Variant.content)
            variantParam = Variant.content{i};
            if strcmp(variantParam(2),'initial_met_diameter')
                % Get Initial Met Diameter from 'variant'
                D_Ti_met = sbioselect(model,'Name','initial_met_diameter');
                met_diameter.Value = variantParam{4};
                met_diameter.Units = D_Ti_met.ValueUnits;
                met_diameter.Notes = D_Ti_met.Notes;
    
                variantWithInitialMetDiameter = true;
            end
        end
        if ~variantWithInitialMetDiameter
            % Get Initial Met Diameter from 'model'
            D_Ti_met = sbioselect(model,'Name','initial_met_diameter');
            met_diameter.Value = D_Ti_met.Value;
            met_diameter.Units = D_Ti_met.ValueUnits;
            met_diameter.Notes = D_Ti_met.Notes;
        end
    else
    % Get Initial Met Diameter from 'model'
    D_Ti_met = sbioselect(model,'Name','initial_met_diameter');
    
    met_diameter.Value = D_Ti_met.Value;
    met_diameter.Units = D_Ti_met.ValueUnits;
    met_diameter.Notes = D_Ti_met.Notes;
    
    end
end

% Calculate Target Tumour Volume
tumour_volume = 4/3*pi*(tumour_diameter/2)^3;

% Calculate Target Met Volume
if(ntumors>1)
    met_volume = 4/3*pi*(met_diameter/2)^3;
end

% Get User-Defined Confiuration Settings
config = getconfigset(model);
user_output_times = config.SolverOptions.OutputTimes;
user_abs_tol = get(config.SolverOptions,'AbsoluteTolerance');
user_rel_tol = get(config.SolverOptions,'RelativeTolerance');

% Reset Output Times for IC Simulation (assumes time in days)
set(config,'StopTime',2000); % thein: changed 5000 to 2000
set(config.SolverOptions,'OutputTimes',[]);
set(config.SolverOptions,'MaxStep',1);

% IC Simulation
if ((Debug)&&(~isempty(Variant)))
    simData = sbiosimulate(model,Variant);
elseif ((Debug)&&(isempty(Variant)))
    simData = sbiosimulate(model);
elseif ((~Debug)&&(isempty(Variant)))
    try
        simData = sbiosimulate(model);
    catch
        % disp('There was an error while finding the initial conditions');
        model_out = copyobj(model);
        success = false;
        simData = [];
        return
    end
elseif ((~Debug)&&(~isempty(Variant)))
    try
        simData = sbiosimulate(model,Variant);
    catch
        % disp('There was an error while finding the initial conditions');
        model_out = copyobj(model);
        success = false;
        simData = [];
        return
    end
end

% Get Time to Reach Target Tumour Size
% Get Tumour Volume from Simulation
for i = 1:size(simData.Data,2)
    if (strcmp(simData.DataInfo{i}.Name,'V_T'))
        idx = i;
        break;
    end
end

V_T = simData.Data(:,idx);
V_T_units = simData.DataInfo{idx}.Units;
% Convert Units of Target Volume to Units in Simulation
zero.Value = 0; % define 'zero' parameter object with same units as simulation
zero.Units = V_T_units;

%thein - for primary tumor
if(tumour_volume.Value==0)
    target_V_T = 0;
else
    tumour_volume = zero + tumour_volume; % operator converts to units of first input
    target_V_T = tumour_volume.Value; % tumour volume in same units as simulation
end

if(ntumors>1)

    % thein: identify the timepoints when mets reach initial conditions 
    % index of all mets, index of primary will be NaN here
    idx_met=NaN(1,ntumors); % index of mets volume
    idx_met_time=NaN(1,ntumors); % Timepoint at which mets volume cross target met volume
    
    for k=2:ntumors
    
        for i = 1:size(simData.Data,2)
    
            if (strcmp(simData.DataInfo{i}.Name,['V_T' TumorArray{k}]))
                idx_met(k) = i;
                break;
            end
    
        end
    
        V_T_met = simData.Data(:,idx_met(k));
        V_T_units_met = simData.DataInfo{idx_met(k)}.Units;
        % Convert Units of Target Volume to Units in Simulation
        zero.Value = 0; % define 'zero' parameter object with same units as simulation
        zero.Units = V_T_units_met;
        
        %thein: target volume for mets
        if(met_volume.Value==0)
            target_met = 0;
        else
            met_volume = zero + met_volume; % operator converts to units of first input
            target_met = met_volume.Value; % met volume in same units as simulation
        end
    
        % Get Time of Target Tumour Size
        if(~isempty(find((target_met-V_T_met)<0,1)))
            idx_met_time(k) = find((target_met-V_T_met)<0,1); % difference should cross zero at target time
        end
    
    end
    
    [value,index]=min(idx_met_time); % index of met tumor that reaches target size first
    
    % find the first timepoint where both primary and selected met each target
    % size first
    
    for i = 1:size(simData.Data,2)
        if (strcmp(simData.DataInfo{i}.Name,['V_T' TumorArray{index}]))
            idx_sel = i;
            break;
        end
    end
    
    V_T_met_selected = simData.Data(:,idx_sel);

end

if(ntumors>1)
    % Get Time of Target Tumour Size
    idx = find((target_V_T-V_T)<0 & (target_met-V_T_met_selected)<0,1); % difference should cross zero at target time
%     disp("ntumors>1");
%     disp(target_V_T);
%     disp(V_T(1));
%     disp(target_met);
%     disp(V_T_met_selected(1));
else
    idx = find((target_V_T-V_T)<0,1); 
%     disp("ntumors=1");
%     disp(target_V_T);
%     disp(V_T(1));
end

% Re-Initialize Model with New ICs
if isempty(idx) % tumour did not reach target size
    success = false;
else
    success = true;
    % Reset User Defined Configuration Settings
    set(config,'StopTime',user_output_times(end));
    set(config.SolverOptions,'OutputTimes',user_output_times);

    %thein: calculate growth rate here
    growth_rates=PSA_estimate_growthrate(simData,patient_id,growth_rates,idx,suffix);

    % Set New ICs for Species
    for i = 1:length(model.Species)
        try
            model.Species(i).InitialAmount = simData.Data(idx,i);
        catch
            model.Species(i).InitialAmount = 0;
            % disp(['Initial Amount of Species ' num2str(i) ' is negative, Inf, or NaN'])
        end
    end
    % Set New ICs for varying Parameters
    for i = 1:length(model.Parameters)
        if ~(isempty(find(strcmp(simData.DataNames,model.Parameters(i).Name ))))
            k = find(strcmp(simData.DataNames, model.Parameters(i).Name));
            model.Parameters(i).Value =  simData.Data(idx,k);
        end
    end
    % Set New ICs for size-changing Compartment
    for i = 1:length(model.Compartments)
        if ~(isempty(find(strcmp(simData.DataNames,model.Compartments(i).Name ))))
            k = find(strcmp(simData.DataNames, model.Compartments(i).Name));
            model.Compartments(i).Capacity =  simData.Data(idx,k);
        end
    end
    % Reset Initial Assignment Rules
    for k = 1:length(model.Rules)
        if strcmp(model.Rules(k).RuleType,'initialAssignment')
            model.Rules(k).Active = false;
        end
    end
end

model_out = copyobj(model);
delete(model)

end