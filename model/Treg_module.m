% Treg Module
%
% Models Treg transport and activation by APCs [Use before antigen corresponding module]
%
% Inputs: model        -- SimBiology model object with four compartments
%         params       -- object containing model parameter Values, Units, and Notes:
%         TumorArray   -- Array of suffixes for tumor comp. names
%         LNArray      -- Array of suffixes for LN comp. names
% Outputs: model -- SimBiology model object with new Tcell module

function model = Treg_module(model,params,n_T_specs,TumorArray,LNArray)

% number of tumors and LNs
ntumors = length(TumorArray);
uniqLNArray = unique(LNArray);
nLNs = length(uniqLNArray);

ID = '0';
% Species Names
species_name = ['T' ID];
antigen = ['P' ID];

% Add Species
% Naive T cells
nCD4_C  = addspecies(sbioselect(model, 'Name', 'V_C'),'nT',params.nCD4_C.Value/params.nCD4_div.Value,'InitialAmountUnits','cell'); 
    set(nCD4_C,'Notes',['Number of naive ' species_name ' cells in the central compartment']);
nCD4_P  = addspecies(sbioselect(model, 'Name', 'V_P'),'nT',params.nCD4_P.Value/params.nCD4_div.Value,'InitialAmountUnits','cell'); 
    set(nCD4_P,'Notes',['Number of naive ' species_name ' cells in the peripheral compartment']);

for array=1:ntumors
    % Mature T cells
    T_T = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'T',0,'InitialAmountUnits','cell'); 
        set(T_T,'Notes',['Number of ' species_name ' cells in the tumor compartment']);
end

for array=1:nLNs

    nCD4_LN = addspecies(sbioselect(model, 'Name', ['V_LN' uniqLNArray{array}]),'nT',params.nCD4_LN.Value/params.nCD4_div.Value,'InitialAmountUnits','cell'); 
        set(nCD4_LN,'Notes',['Number of naive ' species_name ' cells in the lymph node']);
    % Activated T cells
    aT = addspecies(sbioselect(model, 'Name', ['V_LN' uniqLNArray{array}]),'aT',0,'InitialAmountUnits','cell'); 
        set(aT,'Notes',['Number of activated ' species_name ' cells in the lymph node']);
    T_LN = addspecies(sbioselect(model, 'Name', ['V_LN' uniqLNArray{array}]),'T',0,'InitialAmountUnits','cell'); 
        set(T_LN,'Notes',['Number of ' species_name ' cells in the lymph node compartment']);

end

T_C = addspecies(sbioselect(model, 'Name', 'V_C'),'T',0,'InitialAmountUnits','cell'); % params.nTreg_C.Value 
    set(T_C,'Notes',['Number of ' species_name ' cells in the central compartment']);
T_P = addspecies(sbioselect(model, 'Name', 'V_P'),'T',0,'InitialAmountUnits','cell'); % params.nTreg_P.Value 
    set(T_P,'Notes',['Number of ' species_name ' cells in the peripheral compartment']);

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
% H_APC = 'H_APC'; % thein commented this line
for array=1:ntumors
% Add Treg Variable
    addrule(model,['Tregs_' TumorArray{array} ' = V_T' TumorArray{array} '.' species_name],'repeatedAssignment');
    % Antigen Default Hill Function
    p = addparameter(model,['H_' antigen TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes',['Hill function of tumor antigen' antigen]);
end

% Add Parameters
nCD4_div = addparameter(model,'nCD4_div',params.nCD4_div.Value,'ValueUnits',params.nCD4_div.Units);
    set(nCD4_div,'Notes',['T Cell Diversity ',params.nCD4_div.Notes]);
n_clones_slf = addparameter(model,'n_clones_slf',params.n_clones_slf.Value,'ValueUnits',params.n_clones_slf.Units);
    set(n_clones_slf,'Notes',['Number of T cell clones ' params.n_clones_slf.Notes]);
q_nCD4_LN_in = addparameter(model,'q_nCD4_LN_in',params.q_nCD4_LN_in.Value,'ValueUnits',params.q_nCD4_LN_in.Units);
    set(q_nCD4_LN_in,'Notes',['Rate of naive T cell transport into the LN ' params.q_nCD4_LN_in.Notes]);
q_CD4_LN_out = addparameter(model,'q_CD4_LN_out',params.q_CD4_LN_out.Value,'ValueUnits',params.q_CD4_LN_out.Units);
    set(q_CD4_LN_out,'Notes',['Rate of activated T cell transport out of the LN ' params.q_CD4_LN_out.Notes]);
q_nCD4_LN_out = addparameter(model,'q_nCD4_LN_out',params.q_nCD4_LN_out.Value,'ValueUnits',params.q_nCD4_LN_out.Units);
    set(q_nCD4_LN_out,'Notes',['Rate of naive T cell transport out of the LN ' params.q_nCD4_LN_out.Notes]);
k_nCD4_act = addparameter(model,'k_nCD4_act',params.k_nCD4_act.Value,'ValueUnits',params.k_nCD4_act.Units);
    set(k_nCD4_act,'Notes',[species_name ' activation rate ' params.k_nCD4_act.Notes]);
k_CD4_pro = addparameter(model,'k_CD4_pro',params.k_CD4_pro.Value,'ValueUnits',params.k_CD4_pro.Units);
    set(k_CD4_pro,'Notes',[species_name ' proliferation rate ' params.k_CD4_pro.Notes]);
k_CD4_death = addparameter(model,'k_CD4_death',params.k_CD4_death.Value,'ValueUnits',params.k_CD4_death.Units);
    set(k_CD4_death,'Notes',[species_name ' death rate ' params.k_CD4_death.Notes]);
q_CD4_P_in = addparameter(model,'q_CD4_P_in',params.q_CD4_P_in.Value,'ValueUnits',params.q_CD4_P_in.Units);
    set(q_CD4_P_in,'Notes',['rate of ' species_name ' transport into the peripheral compartment ' params.q_CD4_P_in.Notes]);
q_CD4_P_out = addparameter(model,'q_CD4_P_out',params.q_CD4_P_out.Value,'ValueUnits',params.q_CD4_P_out.Units);
    set(q_CD4_P_out,'Notes',['rate of ' species_name ' transport out of the peripheral compartment ' params.q_CD4_P_out.Notes]);
q_CD4_T_in = addparameter(model,'q_CD4_T_in',params.q_CD4_T_in.Value,'ValueUnits',params.q_CD4_T_in.Units);
    set(q_CD4_T_in,'Notes',['rate of ' species_name ' transport into the tumor compartment ' params.q_CD4_T_in.Notes]);
q_nCD4_P_in = addparameter(model,'q_nCD4_P_in',params.q_nCD4_P_in.Value,'ValueUnits',params.q_nCD4_P_in.Units);
    set(q_nCD4_P_in,'Notes',['rate of n' species_name ' transport into the peripheral compartment ' params.q_nCD4_P_in.Notes]);
q_nCD4_P_out = addparameter(model,'q_nCD4_P_out',params.q_nCD4_P_out.Value,'ValueUnits',params.q_nCD4_P_out.Units);
    set(q_nCD4_P_out,'Notes',['rate of n' species_name ' transport out of the peripheral compartment ' params.q_nCD4_P_out.Notes]);
Q_nCD4_thym = addparameter(model,'Q_nCD4_thym',params.Q_nCD4_thym.Value,'ValueUnits',params.Q_nCD4_thym.Units);
    set(Q_nCD4_thym,'Notes',['Thymic output of n',species_name,' ', params.Q_nCD4_thym.Notes]);
k_nCD4_pro = addparameter(model,'k_nCD4_pro',params.k_nCD4_pro.Value,'ValueUnits',params.k_nCD4_pro.Units);
    set(k_nCD4_pro,'Notes',['rate of n',species_name,' proliferation ', params.k_nCD4_pro.Notes]);
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

try % only add once
    k_Treg = addparameter(model,'k_Treg',params.k_Treg.Value,'ValueUnits',params.k_Treg.Units);
        set(k_Treg,'Notes',['Rate of T cell death by Tregs ' params.k_Treg.Notes]);
catch
end

% Add Reactions
% Thymic output of naive T cell
reaction = addreaction(model,'null -> V_C.nT');
    set(reaction,'ReactionRate','Q_nCD4_thym/nCD4_div');
    set(reaction,'Notes','Thymic output of naive T cell to blood');
% Naive T cell proliferation in the peripheral compartment and the LN compartment
reaction = addreaction(model,'null -> V_P.nT');
    set(reaction,'ReactionRate','k_nCD4_pro/nCD4_div*V_P.nT/(K_nT_pro/nCD4_div+V_P.nT)');
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
    set(reaction,'ReactionRate','q_nCD4_P_in*V_C.nT');
    set(reaction,'Notes','Naive T cell entry into the peripheral compartment');
reaction = addreaction(model,'V_P.nT -> V_C.nT');
    set(reaction,'ReactionRate','q_nCD4_P_out*V_P.nT');
    set(reaction,'Notes','Naive T cell exit from the peripheral compartment');
% T cell Death
reaction = addreaction(model,'V_C.T -> null');
    set(reaction,'ReactionRate','k_CD4_death*V_C.T');
    set(reaction,'Notes','T cell death in the central compartment');
reaction = addreaction(model,'V_P.T -> null');
    set(reaction,'ReactionRate','k_CD4_death*V_P.T');
    set(reaction,'Notes','T cell death in the peripheral compartment');
% T cell transport
% Central & Peripheral
reaction = addreaction(model,'V_C.T -> V_P.T');
    set(reaction,'ReactionRate','q_CD4_P_in*V_C.T');
    set(reaction,'Notes','T cell transport into the peripheral compartment');
reaction = addreaction(model,'V_P.T -> V_C.T');
    set(reaction,'ReactionRate','q_CD4_P_out*V_P.T');
    set(reaction,'Notes','T cell transport out of the peripheral compartment');

for array=1:ntumors
    % Naive T cell activation
    reaction = addreaction(model,['V_LN' LNArray{array} '.nT -> null']);
        set(reaction,'ReactionRate',['k_nCD4_act*' 'H_APC' TumorArray{array} '*H_' antigen TumorArray{array} '*V_LN' LNArray{array} '.nT']);
        set(reaction,'Notes','Naive T cell activation');
    reaction = addreaction(model,['null -> V_LN' LNArray{array} '.aT']);
        set(reaction,'ReactionRate',['k_nCD4_act*' 'H_APC' TumorArray{array} '*H_' antigen TumorArray{array} '*V_LN' LNArray{array} '.nT*n_clones_slf']);
        set(reaction,'Notes','Naive T cell activation');
end

for array=1:nLNs   
    reaction = addreaction(model,['V_LN' uniqLNArray{array} '.nT -> null']);
        set(reaction,'ReactionRate',['k_nT_death*V_LN' uniqLNArray{array} '.nT']);
        set(reaction,'Notes','Naive T cell death in the TDLN compartment');
    reaction = addreaction(model,['null -> V_LN' uniqLNArray{array} '.nT']);
        set(reaction,'ReactionRate',['k_nCD4_pro/nCD4_div*V_LN' uniqLNArray{array} '.nT/(K_nT_pro/nCD4_div+V_LN' uniqLNArray{array} '.nT)']);
        set(reaction,'Notes','Naive T cell proliferation in the TDLN compartment');
    % Naive T cell transport into and out of the lymph node
    reaction = addreaction(model,['V_C.nT -> V_LN' uniqLNArray{array} '.nT']);
        set(reaction,'ReactionRate',['q_nCD4_LN_in*V_C.nT']);
        set(reaction,'Notes','Naive T cell entry into the lymph node');
    reaction = addreaction(model,['V_LN' uniqLNArray{array} '.nT -> V_C.nT']);
        set(reaction,'ReactionRate',['q_nCD4_LN_out*V_LN' uniqLNArray{array} '.nT']);
        set(reaction,'Notes','Naive T cell exit from the lymph node');
    % Activated T cell proliferation
    reaction = addreaction(model,['V_LN' uniqLNArray{array} '.aT -> null']);
        set(reaction,'ReactionRate',['k_CD4_pro/N_aT0' uniqLNArray{array} '*V_LN' uniqLNArray{array} '.aT']);
        set(reaction,'Notes',['a' species_name ' cell proliferation']);    
    reaction = addreaction(model,['null -> V_LN' uniqLNArray{array} '.T']);
        set(reaction,'ReactionRate',['k_CD4_pro/N_aT0' uniqLNArray{array} '*2^(N_aT0' uniqLNArray{array} ')*V_LN' uniqLNArray{array} '.aT']); % *(1-H_PD1_APC)
        set(reaction,'Notes',['a' species_name ' cell proliferation']);    
    reaction = addreaction(model,['V_LN' uniqLNArray{array} '.T -> null']);
        set(reaction,'ReactionRate',['k_CD4_death*V_LN' uniqLNArray{array} '.T']);
        set(reaction,'Notes','T cell death in the lymph node compartment');
    % Central & LN
    reaction = addreaction(model,['V_LN' uniqLNArray{array} '.T -> V_C.T']);
        set(reaction,'ReactionRate',['q_CD4_LN_out*V_LN' uniqLNArray{array} '.T']);
        set(reaction,'Notes','T cell transport out of the lymph node compartment');

    % IL2 Reactions
    if (first_call)
        % IL2 Degradation
        reaction = addreaction(model,['V_LN' uniqLNArray{array} '.IL2 -> null']);
            set(reaction,'ReactionRate',['k_IL2_deg*V_LN' uniqLNArray{array} '.IL2*V_LN' uniqLNArray{array}]);
            set(reaction,'Notes','IL2 degradation');
    
         for nt=1:n_T_specs
            % IL2 Consumption
            reaction = addreaction(model,['V_LN' uniqLNArray{array} '.IL2 -> null']);
                set(reaction,'ReactionRate',['k_IL2_cons*V_LN' uniqLNArray{array} '.T' num2str(nt) '*V_LN' uniqLNArray{array} '.IL2/(IL2_50+V_LN' uniqLNArray{array} '.IL2)']);
                set(reaction,'Notes','IL2 consumption by T cells');
         end
    
    end
    
    % IL2 Consumption by Tregs
    reaction = addreaction(model,['V_LN' uniqLNArray{array} '.IL2 -> null']);
        set(reaction,'ReactionRate',['k_IL2_cons*V_LN' uniqLNArray{array} '.T*V_LN' uniqLNArray{array} '.IL2/(IL2_50_Treg+V_LN' uniqLNArray{array} '.IL2)']);
        set(reaction,'Notes','IL2 consumption by Tregs');

end
    
for array=1:ntumors

    reaction = addreaction(model,['V_T' TumorArray{array} '.T -> null']);
        set(reaction,'ReactionRate',['k_CD4_death*V_T' TumorArray{array} '.T']);
        set(reaction,'Notes','T cell death in the tumor compartment');
    % T cell clearance upon Ag clearance
    reaction = addreaction(model,['V_T' TumorArray{array} '.T -> null']);
        set(reaction,'ReactionRate',['k_cell_clear*V_T' TumorArray{array} '.T*(Kc_rec/(C_total' TumorArray{array} '^2 + Kc_rec))']);
        set(reaction,'Notes','T cell clearance upon antigen clearance');
    % Central & tumor
    reaction = addreaction(model,['V_C.T -> V_T' TumorArray{array} '.T']);
        set(reaction,'ReactionRate',['q_CD4_T_in*V_T' TumorArray{array} '*V_C.T*(C_total' TumorArray{array} '^2/(C_total' TumorArray{array} '^2 + Kc_rec))']);
        set(reaction,'Notes','T cell transport into the tumor compartment');
end

for array=1:nLNs
    % Add Rules
    if (first_call)
        % Set Number of Activated T Cell Generations
        p = addparameter(model,['N_aT' uniqLNArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
            set(p,'Notes',['Number of Activated CD8+ T Cell Generations']);
    
        addrule(model,['N_aT' uniqLNArray{array} ' = N0' uniqLNArray{array} ' + N_costim' uniqLNArray{array} '*H_CD28_APC' uniqLNArray{array} ' + N_IL2_CD8' uniqLNArray{array} '*V_LN' uniqLNArray{array} '.IL2/(IL2_50+V_LN' uniqLNArray{array} '.IL2)'],'repeatedAssignment');
    
        p = addparameter(model,['N_aT0' uniqLNArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
            set(p,'Notes',['Number of Activated Treg Generations']);
        addrule(model,['N_aT0' uniqLNArray{array} ' = N0' uniqLNArray{array} ' + N_costim' uniqLNArray{array} '*H_CD28_APC' uniqLNArray{array} ' + N_IL2_CD4' uniqLNArray{array} '*V_LN' uniqLNArray{array} '.IL2/(IL2_50+V_LN' uniqLNArray{array} '.IL2)'],'repeatedAssignment');
    
    end

    % Update Total T Cells in LN (Rule 4)
    %Tcell_rule = model_rules(4); %thein
    Tcell_rule = sbioselect(model,'Where','Rule', 'regexp', ['^T_total_LN' uniqLNArray{array} ' = ']); 
    rule = get(Tcell_rule,'Rule');
    set(Tcell_rule,'Rule',[rule '+V_LN' uniqLNArray{array} '.' species_name]);

end

for array=1:ntumors
    % Get Model Rules for Updating
    % model_rules = get(model,'Rules'); %thein
    
    % Update Total T Cells in tumor (Rule 3)
    %Tcell_rule = model_rules(3); %thein
    Tcell_rule = sbioselect(model,'Where','Rule', 'regexp', ['^T_total' TumorArray{array} ' = ']); 
    rule = get(Tcell_rule,'Rule');
    set(Tcell_rule,'Rule',[rule '+V_T' TumorArray{array} '.' species_name]);

    comp_tumor=sbioselect(model,'Name',['V_T' TumorArray{array}]);
    rename(sbioselect(comp_tumor,'Name','T'),species_name);

end 

for array=1:nLNs

    comp_LN=sbioselect(model,'Name',['V_LN' uniqLNArray{array}]);
    rename(sbioselect(comp_LN,'Name','nT'),['n' species_name]);
    rename(sbioselect(comp_LN,'Name','aT'),['a' species_name]);
    rename(sbioselect(comp_LN,'Name','T'),species_name);

end

% Rename Objects with 'species_name'
rename(nCD4_div,['div_' species_name]);
rename(n_clones_slf,['n_' species_name '_clones']);
rename(nCD4_C,['n' species_name]);
rename(nCD4_P,['n' species_name]);
% rename(nCD4_T,['n' species_name]);
% rename(aT_T,['a' species_name]);
rename(T_C,species_name);
rename(T_P,species_name);
rename(q_CD4_LN_out,['q_' species_name '_LN_out']);
rename(k_nCD4_act,['k_' species_name '_act']);
rename(k_CD4_pro,['k_' species_name '_pro']);
rename(k_CD4_death,['k_' species_name '_death']);
rename(q_CD4_P_in,['q_' species_name '_P_in']);
rename(q_CD4_P_out,['q_' species_name '_P_out']);
rename(q_CD4_T_in,['q_' species_name '_T_in']);
rename(q_nCD4_P_in,['q_n' species_name '_P_in']);
% rename(q_nCD4_T_in,['q_n' species_name '_T_in']);
rename(q_nCD4_LN_in,['q_n' species_name '_LN_in']);
rename(q_nCD4_P_out,['q_n' species_name '_P_out']);
rename(q_nCD4_LN_out,['q_n' species_name '_LN_out']);
rename(Q_nCD4_thym,['Q_n' species_name '_thym']);
rename(k_nCD4_pro,['k_n' species_name '_pro']);
rename(K_nT_pro,['K_n' species_name '_pro']);
rename(k_nT_death,['k_n' species_name '_death']);

warning('off','SimBiology:DimAnalysisNotDone_MatlabFcn_Dimensionless');
