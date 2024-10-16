% Cancer Module
%
% Models cancer cell growth under baseline conditions assuming logistic
% growth and first-order death
%
% Inputs: model        -- simbio model object with tumour compartment, T
%         species_name -- name of cancer cells [must be unique]
%         params       -- object containing model parameter Values, Units, and Notes:
%                         - k_C_growth--cancer growth rate
%                         - C_max--cancer cell capacity
%                         - k_C_death--cancer death rate
%         TumorArray -- Array of suffixes for tumor comp. name
%         LNArray    -- Array of suffixes for LN comp. names
% Outputs: model -- SimBiology model object with new cancer module


%function model = cancer_module(model,species_name,params,TumorArray,LNArray,varargin)
function model = cancer_module(model,species_name,params,TumorArray,LNArray)

% number of tumors, LNs and cancer clones
ntumors = length(TumorArray); 
uniqLNArray = unique(LNArray);
nLNs = length(uniqLNArray);
ncancer_clones=length(species_name);

first_call = true;
%in = inputParser;
%addOptional(in,'initial_amount',1); % initial cancer cell number
%parse(in,varargin{:});
%init = in.Results.initial_amount;

for nc=1:ncancer_clones
    for array=1:ntumors

        % number of cells belonging to each cancer clone in a given tumor
        ncells_clone = addparameter(model,['ncells_' species_name{nc} TumorArray{array}],params.(['ncells_' species_name{nc} TumorArray{array}]).Value,'ValueUnits',params.(['ncells_' species_name{nc} TumorArray{array}]).Units);
            set(ncells_clone,'Notes','numer of cells');

        % Add Species
        if(array==1) % only set the initial value for primary tumor
            Csp = addspecies(sbioselect(model, 'Name', ['V_T' ,TumorArray{array}]),species_name{nc},params.(['ncells_' species_name{nc} TumorArray{array}]).Value,'InitialAmountUnits',params.(['ncells_' species_name{nc} TumorArray{array}]).Units); 
                set(Csp,'Notes','Number of cancer cells in tumour');
        else % Metastatic tumors are seeded later
            Csp = addspecies(sbioselect(model, 'Name', ['V_T' ,TumorArray{array}]),species_name{nc},1e-6,'InitialAmountUnits','cell');
                set(Csp,'Notes','Number of cancer cells in tumour');
        end
        % Add Parameters
        % Growth
        k_C_growth = addparameter(model,['k_' species_name{nc} '_growth' TumorArray{array}],params.(['k_' species_name{nc} '_growth' TumorArray{array}]).Value,'ValueUnits',params.(['k_' species_name{nc} '_growth' TumorArray{array}]).Units,'ConstantValue',false);
            set(k_C_growth,'Notes',['Cancer cell growth rate ' params.(['k_' species_name{nc} '_growth' TumorArray{array}]).Notes]);
        % Death
        k_C_death = addparameter(model,['k_' species_name{nc} '_death' TumorArray{array}],params.(['k_' species_name{nc} '_death' TumorArray{array}]).Value,'ValueUnits',params.(['k_' species_name{nc} '_death' TumorArray{array}]).Units,'ConstantValue',false);
           set(k_C_death,'Notes',['Cancer cell death rate from innate immune cells ' params.(['k_' species_name{nc} '_death' TumorArray{array}]).Notes]);
    end
end

try % only add C_max and vaculature parameter/species once

    % Vasculature parameter
    p = addparameter(model,'k_K_g',params.k_K_g.Value,'ValueUnits',params.k_K_g.Units,'ConstantValue',false); % 5.33
        set(p,'Notes',['Tumour vasculature growth rate ' params.k_K_g.Notes]);
    p = addparameter(model,'k_K_d',params.k_K_d.Value,'ValueUnits',params.k_K_d.Units,'ConstantValue',false); % 7.9e-3
        set(p,'Notes',['Tumour vasculature inhibition rate ' params.k_K_d.Notes]);
    p = addparameter(model,'k_vas_Csec',params.k_vas_Csec.Value,'ValueUnits',params.k_vas_Csec.Units,'ConstantValue',false);
        set(p,'Notes',['Secretion rate of angiogenic factors by cancer cells ' params.k_vas_Csec.Notes]);
    p = addparameter(model,'k_vas_deg',params.k_vas_deg.Value,'ValueUnits',params.k_vas_deg.Units,'ConstantValue',false);
        set(p,'Notes',['Degradation rate of angiogenic factors ' params.k_vas_deg.Notes]);
    p = addparameter(model,'c_vas_50',params.c_vas_50.Value,'ValueUnits',params.c_vas_50.Units,'ConstantValue',false);
        set(p,'Notes',['Half-maximal conc. of angiogenic factor on tumor capacity growth ' params.c_vas_50.Notes]);
    
    for array=1:ntumors    
        genvarname(['C_max' TumorArray{array}]) = addparameter(model,['C_max' TumorArray{array}],params.C_max.Value,'ValueUnits',params.C_max.Units,'ConstantValue',false);
            set(genvarname(['C_max' TumorArray{array}]),'Notes',['Cancer cell capacity ' params.C_max.Notes]);
        % Vasculature species
        s = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'K',params.K0.Value,'InitialAmountUnits',params.K0.Units); % 2.4e8 
            set(s,'Notes',['Maximal tumor capacity ' params.K0.Notes]);
        s = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'c_vas',0,'InitialAmountUnits','picogram/milliliter'); 
            set(s,'Notes','Angiogenic factors ');
       
        addrule(model,['C_max' TumorArray{array} ' = V_T' TumorArray{array} '.K'],'repeatedAssignment');
    end

catch
    first_call = false;
end

% Initial Tumour Diameter
try % only add once
    p = addparameter(model,'initial_tumour_diameter',params.initial_tumour_diameter.Value,'ValueUnits',params.initial_tumour_diameter.Units);
        set(p,'Notes','Pre-treatment tumor diameter');
    p = addparameter(model,'initial_met_diameter',params.initial_met_diameter.Value,'ValueUnits',params.initial_met_diameter.Units);
        set(p,'Notes','Pre-treatment met diameter');
catch
end

for array=1:ntumors
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.c_vas']);
        set(reaction,'ReactionRate',['k_vas_Csec*C_total' TumorArray{array}]);
        set(reaction,'Notes','Secretion of angiogenic factors by cancer cells');
    reaction = addreaction(model,['V_T' TumorArray{array} '.c_vas -> null']);
        set(reaction,'ReactionRate',['k_vas_deg*V_T' TumorArray{array} '.c_vas']);
        set(reaction,'Notes','Degradation of tumor angiogenic factors');
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.K']);
        set(reaction,'ReactionRate',['k_K_g*C_total' TumorArray{array} '*V_T' TumorArray{array} '.c_vas/(V_T' TumorArray{array} '.c_vas+c_vas_50)']);
        set(reaction,'Notes','Growth of tumor carrying capacity');
    reaction = addreaction(model,['V_T' TumorArray{array} '.K -> null']);
        set(reaction,'ReactionRate',['k_K_d*V_T' TumorArray{array} '.K*(nthroot(C_total' TumorArray{array} '/cell*2.57e-6,3))^2']); % 2.57e-6 mm^3/cell is the volume of a cancer cell
        set(reaction,'Notes','Endogenous inhibition of previously generated vasculature');

    addevent(model,['C_total' TumorArray{array} ' < 0.5*cell'],['V_T' TumorArray{array} '.K = 0.01*cell']);

    for nc=1:ncancer_clones
        % Therapy
        param = addparameter(model,['k_' species_name{nc} '_therapy' TumorArray{array}],0,'ValueUnits','1/day','ConstantValue',false);
        set(param,'Notes',['Rate of ' species_name{nc} ' killing by therapy (see Rules)']);
        % Add Reactions
        % Growth
        reaction = addreaction(model,['null -> V_T' TumorArray{array} '.' species_name{nc}]);
            %set(reaction,'ReactionRate',['k_C_growth*V_T' TumorArray{array} '.C*log(max(C_max' TumorArray{array} '/(C_total' TumorArray{array} '+cell), 1))*(1-stp' TumorArray{array} ')']);
           set(reaction,'ReactionRate',['start' TumorArray{array} '*k_' species_name{nc} '_growth' TumorArray{array} '*V_T' TumorArray{array} '.' species_name{nc} '*log(max(C_max' TumorArray{array} '/(C_total' TumorArray{array} '+cell), 1))']);
            set(reaction,'Notes','Cancer cell growth');
        % Death
        reaction = addreaction(model,['V_T' TumorArray{array} '.' species_name{nc} ' -> V_T' TumorArray{array} '.C_x']);
            set(reaction,'ReactionRate',['k_' species_name{nc} '_death' TumorArray{array} '*V_T' TumorArray{array} '.' species_name{nc}]);
            set(reaction,'Notes','Cancer cell death');

        % Tumour Eradication 
        addevent(model,['V_T' TumorArray{array} '.' species_name{nc} ' < 0.5*cell'],['V_T' TumorArray{array} '.' species_name{nc} ' = 0.01*cell']);
        
        % Get Model Rules for Updating
        %model_rules = get(model,'Rules'); %thein: commented
        
        % Update Total Number of Cancer Cells (Rule 2)
        %total_cancer_rule = model_rules(2); %thein: commented
        total_cancer_rule = sbioselect(model,'Where','Rule', 'regexp', ['^C_total' TumorArray{array} ' = ']); 
        rule = get(total_cancer_rule,'Rule');
        set(total_cancer_rule,'Rule',[rule '+ V_T' TumorArray{array} '.' species_name{nc}]);

    end
end 