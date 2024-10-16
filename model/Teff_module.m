% Teff Module
%
% Models Teff transport and activation by APCs [Use before antigen corresponding module]
%
% Inputs: model        -- SimBiology model object with four compartments
%         ID           -- T cell-antigen ID number [must be unique]
%         params       -- object containing model parameter Values, Units, and Notes:
%         cancer_types -- cell array of strings containing names of cancer types Teff kill
%         n_T_specs    -- number of T cell specificities
%         TumorArray   -- Array of suffixes for tumor comp. names
%         LNArray      -- Array of suffixes for LN comp. names
% Outputs: model -- SimBiology model object with new Tcell module

function model = Teff_module(model,ID,params,cancer_types,n_T_specs,TumorArray,LNArray)

% number of tumors
ntumors = length(TumorArray);
uniqLNArray = unique(LNArray);
nLNs = length(uniqLNArray);

% thein: ID is used in one place
% Species Names
species_name = ['T' ID];
%antigen = ['P' ID];

% Add Species
% Naive T cells
nCD8_C  = addspecies(sbioselect(model, 'Name', 'V_C'),'nT',params.nCD8_C.Value/params.nCD8_div.Value,'InitialAmountUnits','cell'); 
    set(nCD8_C,'Notes',['Number of naive ' species_name ' cells in the central compartment']);
nCD8_P  = addspecies(sbioselect(model, 'Name', 'V_P'),'nT',params.nCD8_P.Value/params.nCD8_div.Value,'InitialAmountUnits','cell'); 
    set(nCD8_P,'Notes',['Number of naive ' species_name ' cells in the peripheral compartment']);

for nt=1:n_T_specs
    % Mature T cells
    T_C = addspecies(sbioselect(model, 'Name', 'V_C'),['T' num2str(nt)],0,'InitialAmountUnits','cell');
        set(T_C,'Notes',['Number of ' species_name ' cells in the central compartment']);
    T_P = addspecies(sbioselect(model, 'Name', 'V_P'),['T' num2str(nt)],0,'InitialAmountUnits','cell'); 
        set(T_P,'Notes',['Number of ' species_name ' cells in the peripheral compartment']);
end

for array=1:nLNs
    
    nCD8_LN = addspecies(sbioselect(model, 'Name', ['V_LN' uniqLNArray{array}]),'nT',params.nCD8_LN.Value/params.nCD8_div.Value,'InitialAmountUnits','cell'); 
        set(nCD8_LN,'Notes',['Number of naive ' species_name ' cells in the lymph node']);

    for nt=1:n_T_specs
        % Activated T cells
        aT = addspecies(sbioselect(model, 'Name', ['V_LN' uniqLNArray{array}]),['aT' num2str(nt)],0,'InitialAmountUnits','cell'); 
            set(aT,'Notes',['Number of activated ' species_name ' cells in the lymph node']);
        T_LN = addspecies(sbioselect(model, 'Name', ['V_LN' uniqLNArray{array}]),['T' num2str(nt)],0,'InitialAmountUnits','cell'); 
            set(T_LN,'Notes',['Number of ' species_name ' cells in the lymph node compartment']);
    end
end

for array=1:ntumors
    for nt=1:n_T_specs
        T_T = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),['T' num2str(nt)],0,'InitialAmountUnits','cell'); 
            set(T_T,'Notes',['Number of ' species_name ' cells in the tumor compartment']);
    end
end

% Determine if first call
first_call = true;
try % add IL2 if it does not exist yet
    for array=1:nLNs
        IL2 = addspecies(sbioselect(model, 'Name', ['V_LN' uniqLNArray{array}]),'IL2',1.9e-4,'InitialAmountUnits','nanomolarity'); % PMID: 21774806 
            set(IL2,'Notes','Concentration of IL2 in the lymph node compartment');
    end
catch
    first_call = false;
end

% Add Hill Functions for APC/mAPC
H_APC = 'H_mAPC';

% Add Parameters
nCD8_div = addparameter(model,'nCD8_div',params.nCD8_div.Value,'ValueUnits',params.nCD8_div.Units);
    set(nCD8_div,'Notes',['T Cell Diversity ',params.nCD8_div.Notes]);

for nt=1:n_T_specs
    n_clones_tum = addparameter(model,['n_clones_tum' num2str(nt)],params.(['n_clones_tum' num2str(nt)]).Value,'ValueUnits',params.(['n_clones_tum' num2str(nt)]).Units);
        set(n_clones_tum,'Notes',['Number of T cell clones ' params.(['n_clones_tum' num2str(nt)]).Notes]);
end

q_nCD8_LN_in = addparameter(model,'q_nCD8_LN_in',params.q_nCD8_LN_in.Value,'ValueUnits',params.q_nCD8_LN_in.Units);
    set(q_nCD8_LN_in,'Notes',['Rate of naive T cell transport into the LN ' params.q_nCD8_LN_in.Notes]);
q_CD8_LN_out = addparameter(model,'q_CD8_LN_out',params.q_CD8_LN_out.Value,'ValueUnits',params.q_CD8_LN_out.Units);
    set(q_CD8_LN_out,'Notes',['Rate of activated T cell transport out of the LN ' params.q_CD8_LN_out.Notes]);
q_nCD8_LN_out = addparameter(model,'q_nCD8_LN_out',params.q_nCD8_LN_out.Value,'ValueUnits',params.q_nCD8_LN_out.Units);
    set(q_nCD8_LN_out,'Notes',['Rate of naive T cell transport out of the LN ' params.q_nCD8_LN_out.Notes]);
k_nCD8_act = addparameter(model,'k_nCD8_act',params.k_nCD8_act.Value,'ValueUnits',params.k_nCD8_act.Units);
    set(k_nCD8_act,'Notes',[species_name ' activation rate ' params.k_nCD8_act.Notes]);
k_CD8_pro = addparameter(model,'k_CD8_pro',params.k_CD8_pro.Value,'ValueUnits',params.k_CD8_pro.Units);
    set(k_CD8_pro,'Notes',[species_name ' proliferation rate ' params.k_CD8_pro.Notes]);
k_CD8_death = addparameter(model,'k_CD8_death',params.k_CD8_death.Value,'ValueUnits',params.k_CD8_death.Units);
    set(k_CD8_death,'Notes',[species_name ' death rate ' params.k_CD8_death.Notes]);

for array=1:ntumors

    q_CD8_T_in_tmp = addparameter(model,['q_CD8_T_in' TumorArray{array}],params.(['q_CD8_T_in' TumorArray{array}]).Value,'ValueUnits',params.(['q_CD8_T_in' TumorArray{array}]).Units);
        set(q_CD8_T_in_tmp,'Notes',['rate of ' species_name ' transport into the tumor compartment ' params.(['q_CD8_T_in' TumorArray{array}]).Notes]);
end

q_CD8_P_in = addparameter(model,'q_CD8_P_in',params.q_CD8_P_in.Value,'ValueUnits',params.q_CD8_P_in.Units);
        set(q_CD8_P_in,'Notes',['rate of ' species_name ' transport into the peripheral compartment ' params.q_CD8_P_in.Notes]);
q_CD8_P_out = addparameter(model,'q_CD8_P_out',params.q_CD8_P_out.Value,'ValueUnits',params.q_CD8_P_out.Units);
    set(q_CD8_P_out,'Notes',['rate of ' species_name ' transport out of the peripheral compartment ' params.q_CD8_P_out.Notes]);
q_nCD8_P_in = addparameter(model,'q_nCD8_P_in',params.q_nCD8_P_in.Value,'ValueUnits',params.q_nCD8_P_in.Units);
    set(q_nCD8_P_in,'Notes',['rate of n' species_name ' transport into the peripheral compartment ' params.q_nCD8_P_in.Notes]);
q_nCD8_P_out = addparameter(model,'q_nCD8_P_out',params.q_nCD8_P_out.Value,'ValueUnits',params.q_nCD8_P_out.Units);
    set(q_nCD8_P_out,'Notes',['rate of n' species_name ' transport out of the peripheral compartment ' params.q_nCD8_P_out.Notes]);
Q_nCD8_thym = addparameter(model,'Q_nCD8_thym',params.Q_nCD8_thym.Value,'ValueUnits',params.Q_nCD8_thym.Units);
    set(Q_nCD8_thym,'Notes',['Thymic output of n',species_name,' ', params.Q_nCD8_thym.Notes]);
k_nCD8_pro = addparameter(model,'k_nCD8_pro',params.k_nCD8_pro.Value,'ValueUnits',params.k_nCD8_pro.Units);
    set(k_nCD8_pro,'Notes',['rate of n',species_name,' proliferation ', params.k_nCD8_pro.Notes]);
K_nT_pro = addparameter(model,'K_nT_pro',params.K_nT_pro.Value,'ValueUnits',params.K_nT_pro.Units);
    set(K_nT_pro,'Notes',['half-maximal peripheral proliferation of of n',species_name,' ', params.K_nT_pro.Notes]);
k_nT_death = addparameter(model,'k_nT_death',params.k_nT_death.Value,'ValueUnits',params.k_nT_death.Units);
    set(k_nT_death,'Notes',['rate of n',species_name,' death ', params.k_nT_death.Notes]);

% IL2 Parameters
if (first_call)
    p = addparameter(model,'k_IL2_deg',params.k_IL2_deg.Value,'ValueUnits',params.k_IL2_deg.Units);
        set(p,'Notes',['rate of IL2 degradation ' params.k_IL2_deg.Notes]);
    p = addparameter(model,'k_IL2_cons',params.k_IL2_cons.Value,'ValueUnits',params.k_IL2_cons.Units);
        set(p,'Notes',['rate of IL2 consumption by T cells ' params.k_IL2_cons.Notes]);
    p = addparameter(model,'k_IL2_sec',params.k_IL2_sec.Value,'ValueUnits',params.k_IL2_sec.Units);
        set(p,'Notes',['rate of IL2 secretion from T cells ' params.k_IL2_sec.Notes]);
    p = addparameter(model,'IL2_50',params.IL2_50.Value,'ValueUnits',params.IL2_50.Units);
        set(p,'Notes',['T cell activation half-maximal IL2 concentration ' params.IL2_50.Notes]);
    p = addparameter(model,'IL2_50_Treg',params.IL2_50_Treg.Value,'ValueUnits',params.IL2_50_Treg.Units);
        set(p,'Notes',['Treg activation half-maximal IL2 concentration ' params.IL2_50_Treg.Notes]);
    for array=1:nLNs
        p = addparameter(model,['N0' uniqLNArray{array}],params.N0.Value,'ValueUnits',params.N0.Units);
            set(p,'Notes',['numer of activated T cell generation by TCR signaling only' params.N0.Notes]);
        p = addparameter(model,['N_costim' uniqLNArray{array}],params.N_costim.Value,'ValueUnits',params.N_costim.Units);
            set(p,'Notes',['numer of activated T cell generation by co-stimulatory signaling only ' params.N_costim.Notes]);
        p = addparameter(model,['N_IL2_CD8' uniqLNArray{array}],params.N_IL2_CD8.Value,'ValueUnits',params.N_IL2_CD8.Units);
            set(p,'Notes',['maximum number of activated T cell generations due to IL2 ' params.N_IL2_CD8.Notes]);
        p = addparameter(model,['N_IL2_CD4' uniqLNArray{array}],params.N_IL2_CD4.Value,'ValueUnits',params.N_IL2_CD4.Units);
            set(p,'Notes',['maximum number of activated T cell generations due to IL2 ' params.N_IL2_CD4.Notes]);
    end
end

k_Tcell = addparameter(model,'k_Tcell',params.k_Tcell.Value,'ValueUnits',params.k_Tcell.Units);
    set(k_Tcell,'Notes',['Rate of T cell exhaustion by cancer cells ' params.k_Tcell.Notes]);
k_C_Tcell = addparameter(model,'k_C_Tcell',params.k_C_Tcell.Value,'ValueUnits',params.k_C_Tcell.Units);
    set(k_C_Tcell,'Notes',['Rate of cancer cell death by T cells ' params.k_C_Tcell.Notes]);

try % only add once
    k_Treg = addparameter(model,'k_Treg',params.k_Treg.Value,'ValueUnits',params.k_Treg.Units);
        set(k_Treg,'Notes',['Rate of T cell death by Tregs ' params.k_Treg.Notes]);
catch
end

for array=1:ntumors
    % Antigen Default Hill Function
    for nt=1:n_T_specs
        p = addparameter(model,['H_P' num2str(nt) TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
            set(p,'Notes',['Hill function of tumor antigen P' num2str(nt)]);
    end

end

K_T_C = addparameter(model,'K_T_C',params.K_T_C.Value,'ValueUnits',params.K_T_C.Units);
    set(K_T_C,'Notes',['Dependence of Teff killing rate on Teff/C ratio ',params.K_T_C.Notes]);
K_T_Treg = addparameter(model,'K_T_Treg',params.K_T_Treg.Value,'ValueUnits',params.K_T_Treg.Units);
    set(K_T_Treg,'Notes',['Dependence of Teff killing rate on Teff/Treg ratio ',params.K_T_Treg.Notes]);

% Add Reactions
% Thymic output of naive T cell
reaction = addreaction(model,'null -> V_C.nT');
    set(reaction,'ReactionRate','Q_nCD8_thym/nCD8_div');
    set(reaction,'Notes','Thymic output of naive T cell to blood');
% Naive T cell proliferation in the peripheral compartment and the LN compartment
reaction = addreaction(model,'null -> V_P.nT');
    set(reaction,'ReactionRate','k_nCD8_pro/nCD8_div*V_P.nT/(K_nT_pro/nCD8_div+V_P.nT)');
    set(reaction,'Notes','Naive T cell proliferation in the peripheral compartment');
    % Naive T cell death in the peripheral compartment and the central compartment
reaction = addreaction(model,'V_P.nT -> null');
    set(reaction,'ReactionRate','k_nT_death*V_P.nT');
    set(reaction,'Notes','Naive T cell death in the peripheral compartment');
reaction = addreaction(model,'V_C.nT -> null');
    set(reaction,'ReactionRate','k_nT_death*V_C.nT');
    set(reaction,'Notes','Naive T cell death in the central compartment');
% Naive T cell transport into and out of the peripheral compartment
reaction = addreaction(model,'V_C.nT -> V_P.nT');
    set(reaction,'ReactionRate','q_nCD8_P_in*V_C.nT');
    set(reaction,'Notes','Naive T cell entry into the peripheral compartment');
reaction = addreaction(model,'V_P.nT -> V_C.nT');
    set(reaction,'ReactionRate','q_nCD8_P_out*V_P.nT');
    set(reaction,'Notes','Naive T cell exit from the peripheral compartment');

for nt=1:n_T_specs
    % T cell Death
    reaction = addreaction(model,['V_C.T' num2str(nt) ' -> null']);
        set(reaction,'ReactionRate',['k_CD8_death*V_C.T' num2str(nt)]);
        set(reaction,'Notes','T cell death in the central compartment');
    reaction = addreaction(model,['V_P.T' num2str(nt) ' -> null']);
        set(reaction,'ReactionRate',['k_CD8_death*V_P.T' num2str(nt)]);
        set(reaction,'Notes','T cell death in the peripheral compartment');   
    % T cell transport
    % Central & Peripheral
    reaction = addreaction(model,['V_C.T' num2str(nt) ' -> V_P.T' num2str(nt)]);
        set(reaction,'ReactionRate',['q_CD8_P_in*V_C.T' num2str(nt)]);
        set(reaction,'Notes','T cell transport into the peripheral compartment');
    reaction = addreaction(model,['V_P.T' num2str(nt) ' -> V_C.T' num2str(nt)]);
        set(reaction,'ReactionRate',['q_CD8_P_out*V_P.T' num2str(nt)]);
        set(reaction,'Notes','T cell transport out of the peripheral compartment');
end

for array=1:ntumors 

    for nt=1:n_T_specs
        % Naive T cell activation
        reaction = addreaction(model,['V_LN' LNArray{array} '.nT -> null']);
            set(reaction,'ReactionRate',['k_nCD8_act*' H_APC TumorArray{array} '*H_P' num2str(nt) TumorArray{array} '*V_LN' LNArray{array} '.nT']);
            set(reaction,'Notes','Naive T cell activation');
        reaction = addreaction(model,['null -> V_LN' LNArray{array} '.aT' num2str(nt)]);
            set(reaction,'ReactionRate',['k_nCD8_act*' H_APC TumorArray{array} '*H_P' num2str(nt) TumorArray{array} '*V_LN' LNArray{array} '.nT*n_clones_tum' num2str(nt)]);
            set(reaction,'Notes','Naive T cell activation');
    end
end

for array=1:nLNs
    
    reaction = addreaction(model,['null -> V_LN' uniqLNArray{array} '.nT']);
        set(reaction,'ReactionRate',['k_nCD8_pro/nCD8_div*V_LN' uniqLNArray{array} '.nT/(K_nT_pro/nCD8_div+V_LN' uniqLNArray{array} '.nT)']);
        set(reaction,'Notes','Naive T cell proliferation in the TDLN compartment'); 
    reaction = addreaction(model,['V_LN' uniqLNArray{array} '.nT -> null']);
        set(reaction,'ReactionRate',['k_nT_death*V_LN' uniqLNArray{array} '.nT']);
        set(reaction,'Notes','Naive T cell death in the TDLN compartment');
    % Naive T cell transport into and out of the lymph node
    reaction = addreaction(model,['V_C.nT -> V_LN' uniqLNArray{array} '.nT']);
        set(reaction,'ReactionRate',['q_nCD8_LN_in*V_C.nT']);
        set(reaction,'Notes','Naive T cell entry into the lymph node');
    reaction = addreaction(model,['V_LN' uniqLNArray{array} '.nT -> V_C.nT']);
        set(reaction,'ReactionRate',['q_nCD8_LN_out*V_LN' uniqLNArray{array} '.nT']);
        set(reaction,'Notes','Naive T cell exit from the lymph node');

    for nt=1:n_T_specs
        % Activated T cell proliferation
        reaction = addreaction(model,['V_LN' uniqLNArray{array} '.aT' num2str(nt) ' -> null']);
            set(reaction,'ReactionRate',['k_CD8_pro/N_aT' uniqLNArray{array} '*V_LN' uniqLNArray{array} '.aT' num2str(nt)]);
            set(reaction,'Notes',['a' species_name ' cell proliferation']);
        reaction = addreaction(model,['null -> V_LN' uniqLNArray{array} '.T' num2str(nt)]);
            set(reaction,'ReactionRate',['k_CD8_pro/N_aT' uniqLNArray{array} '*2^(N_aT' uniqLNArray{array} ')*V_LN' uniqLNArray{array} '.aT' num2str(nt)]); % *(1-H_PD1_APC)
            set(reaction,'Notes',['a' species_name ' cell proliferation']);  
        reaction = addreaction(model,['V_LN' uniqLNArray{array} '.T' num2str(nt) ' -> null']);
            set(reaction,'ReactionRate',['k_CD8_death*V_LN' uniqLNArray{array} '.T' num2str(nt)]);
            set(reaction,'Notes','T cell death in the lymph node compartment');
         % Central & LN
        reaction = addreaction(model,['V_LN' uniqLNArray{array} '.T' num2str(nt) ' -> V_C.T' num2str(nt)]);
            set(reaction,'ReactionRate',['q_CD8_LN_out*V_LN' uniqLNArray{array} '.T' num2str(nt)]);
            set(reaction,'Notes','T cell transport out of the lymph node compartment');

    end

    % IL2 Reactions
    if (first_call)
        % IL2 Degradation
        reaction = addreaction(model,['V_LN' uniqLNArray{array} '.IL2 -> null']);
            set(reaction,'ReactionRate',['k_IL2_deg*V_LN' uniqLNArray{array} '.IL2*V_LN' uniqLNArray{array}]);
            set(reaction,'Notes','IL2 degradation');

        for nt=1:n_T_specs
            % IL2 Consumption
            % not used?
            reaction = addreaction(model,['V_LN' uniqLNArray{array} '.IL2 -> null']);
                set(reaction,'ReactionRate',['k_IL2_cons*V_LN' uniqLNArray{array} '.T1' num2str(nt) '*V_LN' uniqLNArray{array} '.IL2/(IL2_50+V_LN' uniqLNArray{array} '.IL2)']);
                set(reaction,'Notes','IL2 consumption by T cells');
        end
    end

    for nt=1:n_T_specs
        % IL2 Secretion by Activated T Cells
        reaction = addreaction(model,['null -> V_LN' uniqLNArray{array} '.IL2']);
            set(reaction,'ReactionRate',['k_IL2_sec*V_LN' uniqLNArray{array} '.aT' num2str(nt)]);
            set(reaction,'Notes','IL2 secretion from activated T cells');
    end

end


for array=1:ntumors

    for nt=1:n_T_specs

        reaction = addreaction(model,['V_T' TumorArray{array} '.T' num2str(nt) ' -> V_T' TumorArray{array} '.T1_exh']);
            set(reaction,'ReactionRate',['k_CD8_death*V_T' TumorArray{array} '.T' num2str(nt)]);
            set(reaction,'Notes','T cell death in the tumor compartment');
        % T cell clearance upon Ag clearance
        reaction = addreaction(model,['V_T' TumorArray{array} '.T' num2str(nt) ' -> V_T' TumorArray{array} '.T1_exh']);
            set(reaction,'ReactionRate',['k_cell_clear*V_T' TumorArray{array} '.T' num2str(nt) '*(Kc_rec/(C_total' TumorArray{array} '^2 + Kc_rec))']);
            set(reaction,'Notes','T cell clearance upon antigen clearance');       
        % T cell death from Treg
        reaction = addreaction(model,['V_T' TumorArray{array} '.T' num2str(nt) ' -> V_T' TumorArray{array} '.T1_exh']);
            set(reaction,'ReactionRate',['k_Treg*V_T' TumorArray{array} '.T' num2str(nt) '*Tregs_' TumorArray{array} '/(V_T' TumorArray{array} '.T' num2str(nt) '+Tregs_' TumorArray{array} '+cell)']);
            set(reaction,'Notes','T cell death from Tregs');       
        % T cell exhaustion from cancer
        reaction = addreaction(model,['V_T' TumorArray{array} '.T' num2str(nt) ' -> V_T' TumorArray{array} '.T1_exh']);
            set(reaction,'ReactionRate',['k_Tcell*V_T' TumorArray{array} '.T' num2str(nt) '*C_total' TumorArray{array} '/(C_total' TumorArray{array} '+V_T' TumorArray{array} '.T' num2str(nt) '+cell)*H_PD1_C1' TumorArray{array}]);
            set(reaction,'Notes','T cell death from cancer');     
        % Central & tumor
        reaction = addreaction(model,['V_C.T' num2str(nt) ' -> V_T' TumorArray{array} '.T' num2str(nt)]);
            set(reaction,'ReactionRate',['q_CD8_T_in' TumorArray{array} '*V_T' TumorArray{array} '*V_C.T' num2str(nt) '*(C_total' TumorArray{array} '^2/(C_total' TumorArray{array} '^2 + Kc_rec))']);
            set(reaction,'Notes','T cell transport into the tumor compartment');
    end
    
end

for array=1:nLNs
    
    % Add Rules
    if (first_call)
        % Set Number of Activated T Cell Generations
        p = addparameter(model,['N_aT' uniqLNArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
            set(p,'Notes','Number of Activated CD8+ T Cell Generations (see Rules)');
    
        addrule(model,['N_aT' uniqLNArray{array} ' = N0' uniqLNArray{array} ' + N_costim' uniqLNArray{array} '*H_CD28_APC' uniqLNArray{array} ' + N_IL2_CD8' uniqLNArray{array} '*V_LN' uniqLNArray{array} '.IL2/(IL2_50+V_LN' uniqLNArray{array} '.IL2)'],'repeatedAssignment');
    
        p = addparameter(model,['N_aT0' uniqLNArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
            set(p,'Notes','Number of Activated Treg Generations (see Rules)');
    
        addrule(model,['N_aT0' uniqLNArray{array} ' = N0' uniqLNArray{array} ' + N_costim' uniqLNArray{array} '*H_CD28_APC' uniqLNArray{array} ' + N_IL2_CD4' uniqLNArray{array} '*V_LN' uniqLNArray{array} '.IL2/(IL2_50+V_LN' uniqLNArray{array} '.IL2)'],'repeatedAssignment');
    
    end

    % Update Total T Cells in LN (Rule 4)
    %Tcell_rule = model_rules(4); %thein
    Tcell_rule = sbioselect(model,'Where','Rule', 'regexp', ['^T_total_LN' uniqLNArray{array} ' = ']); 

    for nt=1:n_T_specs
        rule = get(Tcell_rule,'Rule');
        set(Tcell_rule,'Rule',[rule '+V_LN' uniqLNArray{array} '.T' num2str(nt)]);
    end

end
    
for array=1:ntumors
    % Get Model Rules for Updating
    %model_rules = get(model,'Rules'); %thein
    
    % Update Total T Cells in tumor (Rule 3)
    %Tcell_rule = model_rules(3); %thein
    Tcell_rule = sbioselect(model,'Where','Rule', 'regexp', ['^T_total' TumorArray{array} ' = ']); 

    for nt=1:n_T_specs
        rule = get(Tcell_rule,'Rule');
        set(Tcell_rule,'Rule',[rule '+V_T' TumorArray{array} '.T' num2str(nt)]);
    end
    
    T1cell_rule = sbioselect(model,'Where','Rule', 'regexp', ['^Tcyt_total' TumorArray{array} ' = ']); 

    for nt=1:n_T_specs
        rule = get(T1cell_rule,'Rule');
        set(T1cell_rule,'Rule',[rule '+V_T' TumorArray{array} '.T' num2str(nt)]);
    end


    ruleobj=sbioselect(model,'Where','Rule', 'regexp', ['^R_Tcell' TumorArray{array} ' = ']); 
    
    % Update Cancer Killing by T Cells (Rule 5)
    if (exist('cancer_types','var'))
        for i = 1:length(cancer_types)
            reaction = addreaction(model,['V_T' TumorArray{array} '.' cancer_types{i} ' -> V_T' TumorArray{array} '.C_x']);
                set(reaction,'ReactionRate',['k_C_Tcell*V_T' TumorArray{array} '.' cancer_types{i} '*' cancer_types{i} '_Tcyt_total' TumorArray{array} '/(K_T_C*V_T' TumorArray{array} '.' cancer_types{i} '+' cancer_types{i} '_Tcyt_total' TumorArray{array} '+cell)*' cancer_types{i} '_Tcyt_total' TumorArray{array} '/(' cancer_types{i} '_Tcyt_total' TumorArray{array} '+K_T_Treg*Tregs_' TumorArray{array} '+cell)*(1-H_TGFb_Teff' TumorArray{array} ')*(1-H_PD1_C1' TumorArray{array} ')']);                
                set(reaction,'Notes','Cancer cell killing by T cells');
               % rule = get(model_rules(5),'Rule'); %thein
                rule = get(ruleobj,'Rule'); 
                set(ruleobj,'Rule',[rule,'+(k_' cancer_types{i} '_death' TumorArray{array} '+k_' cancer_types{i} '_therapy' TumorArray{array} ')*V_T' TumorArray{array} '.' cancer_types{i} '+k_C_Tcell*V_T' TumorArray{array} '.' cancer_types{i} '*' cancer_types{i} '_Tcyt_total' TumorArray{array} '/(K_T_C*V_T' TumorArray{array} '.' cancer_types{i} '+' cancer_types{i} '_Tcyt_total' TumorArray{array} '+cell)*' cancer_types{i} '_Tcyt_total' TumorArray{array} '/(' cancer_types{i} '_Tcyt_total' TumorArray{array} '+K_T_Treg*Tregs_' TumorArray{array} '+cell)*(1-H_TGFb_Teff' TumorArray{array} ')*(1-H_PD1_C1' TumorArray{array} ')']);
                           
        end
    end

    comp_tumor=sbioselect(model,'Name',['V_T' TumorArray{array}]);

end

 for array=1:nLNs

    comp_LN=sbioselect(model,'Name',['V_LN' uniqLNArray{array}]);
    rename(sbioselect(comp_LN,'Name','nT'),['n' species_name]);

end
% Rename Objects with 'species_name'
rename(nCD8_div,['div_' species_name]);
%rename(n_clones_tum,['n_' species_name '_clones']);

for nt=1:n_T_specs
    rename(sbioselect(model,'Name',['n_clones_tum' num2str(nt)]),['n_T' num2str(nt) '_clones']);
end

for array=1:ntumors
    rename(sbioselect(model,'Name',['q_CD8_T_in' TumorArray{array}]),['q_' species_name '_T_in' TumorArray{array}]);
end

rename(nCD8_C,['n' species_name]);
rename(nCD8_P,['n' species_name]);
rename(q_CD8_LN_out,['q_' species_name '_LN_out']);
rename(k_nCD8_act,['k_' species_name '_act']);
rename(k_CD8_pro,['k_' species_name '_pro']);
rename(k_CD8_death,['k_' species_name '_death']);
rename(q_CD8_P_in,['q_' species_name '_P_in']);
rename(q_CD8_P_out,['q_' species_name '_P_out']);
%rename(q_CD8_T_in,['q_' species_name '_T_in']);
rename(q_nCD8_P_in,['q_n' species_name '_P_in']);
% rename(q_nCD8_T_in,['q_n' species_name '_T_in']);
rename(q_nCD8_LN_in,['q_n' species_name '_LN_in']);
rename(q_nCD8_P_out,['q_n' species_name '_P_out']);
rename(q_nCD8_LN_out,['q_n' species_name '_LN_out']);
rename(Q_nCD8_thym,['Q_n' species_name '_thym']);
rename(k_nCD8_pro,['k_n' species_name '_pro']);
rename(K_nT_pro,['K_n' species_name '_pro']);
rename(k_nT_death,['k_n' species_name '_death']);
rename(k_Tcell,['k_' species_name]);
rename(k_C_Tcell,['k_C_' species_name]);

warning('off','SimBiology:DimAnalysisNotDone_MatlabFcn_Dimensionless');
