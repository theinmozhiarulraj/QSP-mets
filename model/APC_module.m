% APC Module
%
% Models APC activation and transport
%
% Inputs: model        -- simbio model object with four compartments
%         species_name -- name of Tcells [must be unique]
%         params       -- object containing model parameter Values, Units, and Notes:
%                         - k_APC_mat--rate of APC maturation
%                         - k_APC_mig--rate of APC migration
%                         - k_APC_dth--rate of APC death
%                         - k_mAPC_dth--rate of mAPC death
%                         - APC0_T--baseline APC density in tumour
%                         - APC0_LN--baseline APC density in LN
%                         - k_c--cytokine time constant
%                         - c0--baseline cytokine concentration
%                         - c50--cytokine concentration for half-maximal APC maturation
%                         - DAMPS--concentration of cytokines released per dying cancer cell
%        TumorArray -- Array of suffixes for tumor comp. name
%        LNArray    -- Array of suffixes for LN comp. name
% Outputs: model -- SimBiology model object with new APC module


function model = APC_module(model,params,TumorArray,LNArray)

% number of tumors and LNs
ntumors= length(TumorArray);
uniqLNArray = unique(LNArray);
nLNs = length(uniqLNArray);

for array=1:ntumors
    % Add Species
    APC = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'APC',0,'InitialAmountUnits','cell'); 
        set(APC,'Notes','Number of naive antigen presenting cells in the tumour');
    mAPC = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'mAPC',0,'InitialAmountUnits','cell'); 
        set(mAPC,'Notes','Number of mature antigen presenting cells in the tumour');
    c = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'c',0,'InitialAmountUnits','molarity'); 
        set(c,'Notes','Concentration of maturation cytokines in the tumour');
    mAPC = addspecies(sbioselect(model, 'Name', ['V_LN' LNArray{array}]),['mAPC' TumorArray{array}],0,'InitialAmountUnits','cell'); 
    set(mAPC,'Notes','Number of mature antigen presenting cells in the lymph node');
end

for array=1:nLNs
    APC = addspecies(sbioselect(model, 'Name', ['V_LN' uniqLNArray{array}]),'APC',0,'InitialAmountUnits','cell'); 
    set(APC,'Notes','Number of naive antigen presenting cells in the lymph node');
end

% Add Parameters
k_APC_mat = addparameter(model,'k_APC_mat',params.k_APC_mat.Value,'ValueUnits',params.k_APC_mat.Units);
    set(k_APC_mat,'Notes',['Maximum rate of APC maturation ' params.k_APC_mat.Notes]);
k_APC_mig = addparameter(model,'k_APC_mig',params.k_APC_mig.Value,'ValueUnits',params.k_APC_mig.Units);
    set(k_APC_mig,'Notes',['Rate of APC migration ' params.k_APC_mig.Notes]);
k_APC_death = addparameter(model,'k_APC_death',params.k_APC_death.Value,'ValueUnits',params.k_APC_death.Units);
    set(k_APC_death,'Notes',['Rate of APC death ' params.k_APC_death.Notes]);
k_mAPC_death = addparameter(model,'k_mAPC_death',params.k_mAPC_death.Value,'ValueUnits',params.k_mAPC_death.Units);
    set(k_mAPC_death,'Notes',['Rate of mAPC death' params.k_mAPC_death.Notes]);
for array=1:ntumors
    APC0_T_temp = addparameter(model,['APC0_T' TumorArray{array}],params.(['APC0_T' TumorArray{array}]).Value,'ValueUnits',params.(['APC0_T' TumorArray{array}]).Units);
        set(APC0_T_temp,'Notes','Steady state density of CD103+ APCs in the tumor, calculated to fit basal APC equal to 26000/cm^3 from (Lavin 2017, PMID: 28475900)');
end
APC0_LN = addparameter(model,'APC0_LN',params.APC0_LN.Value,'ValueUnits',params.APC0_LN.Units);
    set(APC0_LN,'Notes','Steady state density of all APC subtypes in the LN (Catron 2004, PMID: 15357945)- divided by 10 to only take CD103+ based on (Lavin 2017, PMID: 28475900)');
k_c = addparameter(model,'k_c',params.k_c.Value,'ValueUnits',params.k_c.Units);
    set(k_c,'Notes',['Cytokine rate constant ' params.k_c.Notes]);
c0 = addparameter(model,'c0',params.c0.Value,'ValueUnits',params.c0.Units);
    set(c0,'Notes',['Baseline cytokine concentration ' params.c0.Notes]);
c50 = addparameter(model,'c50',params.c50.Value,'ValueUnits',params.c50.Units);
    set(c50,'Notes',['Cytokine concentration for half-maximal APC maturation ' params.c50.Notes]);
DAMPs = addparameter(model,'DAMPs',params.DAMPs.Value,'ValueUnits',params.DAMPs.Units);
    set(DAMPs,'Notes',['Concentration of cytokines released per dying cancer cell ' params.DAMPs.Notes]);
n_sites = addparameter(model,'n_sites_APC',params.n_sites_APC.Value,'ValueUnits',params.n_sites_APC.Units);
    set(n_sites,'Notes',['Maxium number of T Cells an APC can interact with ' params.n_sites_APC.Notes]);

for array=1:ntumors
    % Add Reactions
    % APC Reacruitment/Death in the Tumour
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.APC']);
        set(reaction,'ReactionRate',['k_APC_death*(APC0_T' TumorArray{array} '*V_T' TumorArray{array} '-V_T' TumorArray{array} '.APC)']);
        set(reaction,'Notes','APC recruitment/death in the tumour');
    % APC Maturation
    reaction = addreaction(model,['V_T' TumorArray{array} '.APC -> V_T' TumorArray{array} '.mAPC']);
        set(reaction,'ReactionRate',['k_APC_mat*V_T' TumorArray{array} '.c/(V_T' TumorArray{array} '.c+c50)*V_T' TumorArray{array} '.APC']);
        set(reaction,'Notes','APC maturation in the tumour');   
    % APC Migration to LN
    reaction = addreaction(model,['V_T' TumorArray{array} '.mAPC -> V_LN' LNArray{array} '.mAPC' TumorArray{array}]);
        set(reaction,'ReactionRate',['k_APC_mig*V_T' TumorArray{array} '.mAPC']);
        set(reaction,'Notes','APC migration to the lymph node');    
    % mAPC Death in T
    reaction = addreaction(model,['V_T' TumorArray{array} '.mAPC -> null']);
        set(reaction,'ReactionRate',['k_mAPC_death*V_T' TumorArray{array} '.mAPC']);
        set(reaction,'Notes','mAPC death in the tumour');
    % Baseline Cytokine Secretion/Degradation
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.c']);
        set(reaction,'ReactionRate',['k_c*(c0-V_T' TumorArray{array} '.c)']);
        set(reaction,'Notes','Baseline cytokine secretion/degradation');    
    % Cytokine Release in Response to the Tumour
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.c']);
        set(reaction,'ReactionRate',['R_Tcell' TumorArray{array} '*DAMPs']);
        set(reaction,'Notes','Cytokine release in response to the tumour');
end


for array=1:nLNs
    % thein: probably unused equation
    % APC Recruitment/Death in LN
    reaction = addreaction(model,['null -> V_LN' uniqLNArray{array} '.APC']);
        set(reaction,'ReactionRate',['k_APC_death*(APC0_LN*V_LN' uniqLNArray{array} '-V_LN' uniqLNArray{array} '.APC)']);
        set(reaction,'Notes','APC recruitment/death in LN');
end

for array=1:ntumors
    % mAPC Death in LN
    reaction = addreaction(model,['V_LN' LNArray{array} '.mAPC' TumorArray{array} ' -> null']);
        set(reaction,'ReactionRate',['k_mAPC_death*V_LN' LNArray{array} '.mAPC' TumorArray{array}]);
        set(reaction,'Notes','mAPC death in the lymph node');

    % Update Hill Functions for APCs
    addrule(model,['H_APC' TumorArray{array} ' = n_sites_APC*V_LN' LNArray{array} '.mAPC' TumorArray{array} '/(n_sites_APC*V_LN' LNArray{array} '.mAPC' TumorArray{array} '+V_LN' LNArray{array} '.nT0*n_T0_clones+cell)'],'repeatedAssignment');
    addrule(model,['H_mAPC' TumorArray{array} ' = n_sites_APC*V_LN' LNArray{array} '.mAPC' TumorArray{array} '/(n_sites_APC*V_LN' LNArray{array} '.mAPC' TumorArray{array} '+V_LN' LNArray{array} '.nT1*n_T1_clones+cell)'],'repeatedAssignment');
    addrule(model,['H_APCh' TumorArray{array} ' = n_sites_APC*V_LN' LNArray{array} '.mAPC' TumorArray{array} '/(n_sites_APC*V_LN' LNArray{array} '.mAPC' TumorArray{array} '+V_LN' LNArray{array} '.nT0*n_T1_clones+cell)'],'repeatedAssignment');
end

    

