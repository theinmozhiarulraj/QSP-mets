% Helper T Cell Module
%
% Models Th activation and transport
%
% Inputs: model        -- SimBiology model object with four compartments
%         params       -- object containing the default parameters
%         TumorArray   -- Array of suffixes for tumor comp. names
%         LNArray      -- Array of suffixes for LN comp. names
% Outputs: model -- SimBiology model object with new Th module

function model = Th_module(model,params,TumorArray,LNArray)

% number of tumors
ntumors=length(TumorArray);
uniqLNArray = unique(LNArray);
nLNs = length(uniqLNArray);

% nT0 is naive CD4+, T0 is CD4+/Treg
species_name = ['Th'];
n_clones = 'n_T1_clones'; % number of tumor neoantigen clones

% Add Species
% Mature T cells
T_C = addspecies(sbioselect(model, 'Name', 'V_C'),'T',0,'InitialAmountUnits','cell'); 
    set(T_C,'Notes',['Number of ' species_name ' cells in the central compartment']);
T_P = addspecies(sbioselect(model, 'Name', 'V_P'),'T',0,'InitialAmountUnits','cell'); 
    set(T_P,'Notes',['Number of ' species_name ' cells in the peripheral compartment']);

for array=1:ntumors

    T_T = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'T',0,'InitialAmountUnits','cell'); 
        set(T_T,'Notes',['Number of ' species_name ' cells in the tumor compartment']);
    s = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'TGFb',0,'InitialAmountUnits','nanomolarity'); 
        set(s,'Notes',['TGFb in the tumor comparment']);
    s = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'IFNg',0,'InitialAmountUnits','nanomolarity'); 
        set(s,'Notes',['IFNg in the tumor comparment']);

end

for array=1:nLNs
    % Activated T cells
    aT = addspecies(sbioselect(model, 'Name', ['V_LN' uniqLNArray{array}]),'aT',0,'InitialAmountUnits','cell'); 
        set(aT,'Notes',['Number of activated ' species_name ' cells in the lymph node']);

    T_LN = addspecies(sbioselect(model, 'Name', ['V_LN' uniqLNArray{array}]),'T',0,'InitialAmountUnits','cell'); 
        set(T_LN,'Notes',['Number of ' species_name ' cells in the lymph node compartment']);
end    

% Add Hill Functions for APC/mAPC
H_APC = 'H_APCh';

% Add Parameters
k_Th_act = addparameter(model,'k_Th_act',params.k_Th_act.Value,'ValueUnits',params.k_Th_act.Units);
    set(k_Th_act,'Notes',['Rate of T helper cell activation ' params.k_Th_act.Notes]);
k_Th_Treg = addparameter(model,'k_Th_Treg',params.k_Th_Treg.Value,'ValueUnits',params.k_Th_Treg.Units);
    set(k_Th_Treg,'Notes',[species_name ' differentiation rate to Treg ' params.k_Th_Treg.Notes]);
k_TGFb_Tsec = addparameter(model,'k_TGFb_Tsec',params.k_TGFb_Tsec.Value,'ValueUnits',params.k_TGFb_Tsec.Units); % 1.2e-10 12193750; 23393376 27589056 9697990
    set(k_TGFb_Tsec,'Notes',['TGF secretion rate by Treg ' params.k_TGFb_Tsec.Notes]);
k_TGFb_deg = addparameter(model,'k_TGFb_deg',params.k_TGFb_deg.Value,'ValueUnits',params.k_TGFb_deg.Units);
    set(k_TGFb_deg,'Notes',['TGF degradtion rate ' params.k_TGFb_deg.Notes]);
TGFb_50 = addparameter(model,'TGFb_50',params.TGFb_50.Value,'ValueUnits',params.TGFb_50.Units);
    set(TGFb_50,'Notes',['Half-Maximal TGFb level for Th-to-Treg differentiation / chemoresistance development / M1-to-M2 polarization ' params.TGFb_50.Notes]);
TGFb_50_Teff = addparameter(model,'TGFb_50_Teff',params.TGFb_50_Teff.Value,'ValueUnits',params.TGFb_50_Teff.Units);
    set(TGFb_50_Teff,'Notes',['Half-Maximal TGFb level for CD8 T cell inhibition ' params.TGFb_50_Teff.Notes]);
Kc_rec = addparameter(model,'Kc_rec',params.Kc_rec.Value,'ValueUnits',params.Kc_rec.Units);
    set(Kc_rec,'Notes',['Half-Maximal cancer cell number for T cell recruitment ' params.Kc_rec.Notes]);
TGFbase = addparameter(model,'TGFbase',params.TGFbase.Value,'ValueUnits',params.TGFbase.Units);
    set(TGFbase,'Notes',['Baseline TGFb level in breast tumor ' params.TGFbase.Notes]);
k_IFNg_sec = addparameter(model,'k_IFNg_sec',params.k_IFNg_sec.Value,'ValueUnits',params.k_IFNg_sec.Units);
    set(k_IFNg_sec,'Notes',['IFNg secretion rate by T helper cell ' params.k_IFNg_sec.Notes]);
k_IFNg_deg = addparameter(model,'k_IFNg_deg',params.k_IFNg_deg.Value,'ValueUnits',params.k_IFNg_deg.Units);
    set(k_IFNg_deg,'Notes',['IFNg degradation rate ' params.k_IFNg_deg.Notes]);
IFNg_50_ind = addparameter(model,'IFNg_50_ind',params.IFNg_50_ind.Value,'ValueUnits',params.IFNg_50_ind.Units);
    set(IFNg_50_ind,'Notes',['Half-Maximal IFNg level for PD-L1 induction ' params.IFNg_50_ind.Notes]);

for array=1:ntumors    
    p = addparameter(model,['H_TGFb' TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes','Hill function of TGFb for Th-to-Treg differentiation / chemoresistance development / M1-to-M2 polarization');
    addrule(model,['H_TGFb' TumorArray{array} ' = V_T' TumorArray{array} '.TGFb/(V_T' TumorArray{array} '.TGFb+TGFb_50)'],'repeatedAssignment');
    p = addparameter(model,['H_TGFb_Teff' TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes','Hill function of TGFb for Teff inhibition and chemoresistance development');
    addrule(model,['H_TGFb_Teff' TumorArray{array} ' = V_T' TumorArray{array} '.TGFb/(V_T' TumorArray{array} '.TGFb+TGFb_50_Teff)'],'repeatedAssignment');
end

% Add Reactions
% T cell Death
reaction = addreaction(model,'V_C.T -> null');
    set(reaction,'ReactionRate','k_T0_death*V_C.T');
    set(reaction,'Notes','T cell death in the central compartment');
reaction = addreaction(model,'V_P.T -> null');
    set(reaction,'ReactionRate','k_T0_death*V_P.T');
    set(reaction,'Notes','T cell death in the peripheral compartment');
% T cell transport
% Central & Peripheral
reaction = addreaction(model,'V_C.T -> V_P.T');
    set(reaction,'ReactionRate','q_T0_P_in*V_C.T');
    set(reaction,'Notes','T cell transport into the peripheral compartment');
reaction = addreaction(model,'V_P.T -> V_C.T');
    set(reaction,'ReactionRate','q_T0_P_out*V_P.T');
    set(reaction,'Notes','T cell transport out of the peripheral compartment');

for array=1:ntumors
        % Naive T cell activation
    reaction = addreaction(model,['V_LN' LNArray{array} '.nT0 -> null']);
        set(reaction,'ReactionRate',['k_Th_act*' H_APC TumorArray{array} '*H_P1' TumorArray{array} '*V_LN' LNArray{array} '.nT0']);
        set(reaction,'Notes','Naive T cell activation');
    reaction = addreaction(model,['null -> V_LN' LNArray{array} '.aT']);
        set(reaction,'ReactionRate',['k_Th_act*' H_APC TumorArray{array} '*H_P1' TumorArray{array} '*V_LN' LNArray{array} '.nT0*' n_clones]);
        set(reaction,'Notes','Naive T cell activation');
end

for array=1:nLNs
    % Activated T cell proliferation
    reaction = addreaction(model,['V_LN' uniqLNArray{array} '.aT -> null']);
        set(reaction,'ReactionRate',['(k_T0_pro/N_aTh' uniqLNArray{array} ')*V_LN' uniqLNArray{array} '.aT']);
        set(reaction,'Notes',['a' species_name ' cell proliferation']);
    reaction = addreaction(model,['null -> V_LN' uniqLNArray{array} '.T']);
        set(reaction,'ReactionRate',['(k_T0_pro/N_aTh' uniqLNArray{array} ')*2^N_aTh' uniqLNArray{array} '*V_LN' uniqLNArray{array} '.aT']); % *(1-H_PD1_APC)
        set(reaction,'Notes',['a' species_name ' cell proliferation']);
    reaction = addreaction(model,['V_LN' uniqLNArray{array} '.T -> null']);
        set(reaction,'ReactionRate',['k_T0_death*V_LN' uniqLNArray{array} '.T']);
        set(reaction,'Notes','T cell death in the lymph node compartment');
    % Central & LN
    reaction = addreaction(model,['V_LN' uniqLNArray{array} '.T -> V_C.T']);
        set(reaction,'ReactionRate',['q_T0_LN_out*V_LN' uniqLNArray{array} '.T']);
        set(reaction,'Notes','T cell transport out of the lymph node compartment');
    % IL2 Reactions
    % IL2 Secretion by Activated T Cells
    reaction = addreaction(model,['null -> V_LN' uniqLNArray{array} '.IL2']);
       set(reaction,'ReactionRate',['k_IL2_sec*V_LN' uniqLNArray{array} '.aT']);
       set(reaction,'Notes','IL2 secretion from activated T cells');
end
    
for array=1:ntumors
    % Differentiation between Treg and Th Cells
    reaction = addreaction(model,['V_T' TumorArray{array} '.T -> V_T' TumorArray{array} '.T0']);
        set(reaction,'ReactionRate',['k_Th_Treg*V_T' TumorArray{array} '.T*H_TGFb' TumorArray{array}]);
        set(reaction,'Notes','Differentiation between Treg and Th Cells'); 
    reaction = addreaction(model,['V_T' TumorArray{array} '.T -> V_T' TumorArray{array} '.Th_exh']);
        set(reaction,'ReactionRate',['k_T0_death*V_T' TumorArray{array} '.T']);
        set(reaction,'Notes','T cell death in the tumor compartment');    
    reaction = addreaction(model,['V_T' TumorArray{array} '.T -> V_T' TumorArray{array} '.Th_exh']);
        set(reaction,'ReactionRate',['k_cell_clear*V_T' TumorArray{array} '.T*(Kc_rec/(C_total' TumorArray{array} '^2 + Kc_rec))']);
        set(reaction,'Notes','T cell clearance upon antigen clearance');    
    % Central & tumor
    reaction = addreaction(model,['V_C.T -> V_T' TumorArray{array} '.T']);
        set(reaction,'ReactionRate',['q_T0_T_in*V_T' TumorArray{array} '*V_C.T*(C_total' TumorArray{array} '^2/(C_total' TumorArray{array} '^2 + Kc_rec))']);
        set(reaction,'Notes','T cell transport into the tumor compartment');   
    % TGFb Reactions
    % TGFb secretion and consumption by cancer cells
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.TGFb']);
        set(reaction,'ReactionRate',['k_TGFb_deg*(TGFbase - V_T' TumorArray{array} '.TGFb)*V_T' TumorArray{array}]);
        set(reaction,'Notes','TGFb secretion by triple-negative breast cancer cells');
    % TGFb Secretion by Activated Treg Cells
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.TGFb']);
       set(reaction,'ReactionRate',['k_TGFb_Tsec*V_T' TumorArray{array} '.T0']);
       set(reaction,'Notes','TGFb secretion from activated T cells in tumor');   
    % IFNg Secretion by CD4 T helper Cells
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.IFNg']);
        set(reaction,'ReactionRate',['k_IFNg_sec*V_T' TumorArray{array} '.T']);
        set(reaction,'Notes','IFNg secretion from T helper cells in tumor');
    % IFNg Degradation
    reaction = addreaction(model,['V_T' TumorArray{array} '.IFNg -> null']);
        set(reaction,'ReactionRate',['k_IFNg_deg*V_T' TumorArray{array} '.IFNg']);
        set(reaction,'Notes','IFNg degradtion in tumor');
end
    
for array=1:nLNs

    % Add Rules
    % Set Number of Activated T Cell Generations
    p = addparameter(model,['N_aTh' uniqLNArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes',['Number of Activated T Helper Cell Generations (see Rules)']);
    addrule(model,['N_aTh' uniqLNArray{array} ' = N0' uniqLNArray{array} ' + N_costim' uniqLNArray{array} '*H_CD28_APC' uniqLNArray{array} ' + N_IL2_CD4' uniqLNArray{array} '*V_LN' uniqLNArray{array} '.IL2/(IL2_50+V_LN' uniqLNArray{array} '.IL2)'],'repeatedAssignment');
    
    % Update Total T Cells in LN (Rule 4)
    %Tcell_rule = model_rules(4); %thein
    Tcell_rule = sbioselect(model,'Where','Rule', 'regexp', ['^T_total_LN' uniqLNArray{array} ' = ']); 
    rule = get(Tcell_rule,'Rule');
    set(Tcell_rule,'Rule',[rule '+V_LN' uniqLNArray{array} '.' species_name]);

end

for array=1:ntumors
    % Get Model Rules for Updating
    %model_rules = get(model,'Rules'); %thein
    
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
    rename(sbioselect(comp_LN,'Name','aT'),['a' species_name]);
    rename(sbioselect(comp_LN,'Name','T'),species_name);

end

% rename(aT_T,['a' species_name]);
rename(T_C,species_name);
rename(T_P,species_name);

warning('off','SimBiology:DimAnalysisNotDone_MatlabFcn_Dimensionless');
