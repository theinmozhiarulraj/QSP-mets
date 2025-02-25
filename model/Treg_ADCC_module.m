% Treg CTLA4 Module
%
% Models Treg depletion through antibody binding to CTLA4
% [Use after checkpoint and Treg module]
%
% Inputs: model        -- SimBiology model object with four compartments
%         params       -- object containing the default parameters
%         TumorArray   -- Array of suffixes for tumor comp. names
%         LNArray      -- Array of suffixes for LN comp. names
% Outputs: model -- SimBiology model object with new Treg-CTLA4 module

function model = Treg_ADCC_module(model,params,TumorArray,LNArray)

% number of tumors and LNs
ntumors = length(TumorArray);
uniqLNArray = unique(LNArray);
nLNs = length(uniqLNArray);

% Determine if some of the prameters that could have been defined in Checkpoint module are defined
first_call = true;
try % Add kon Values
    kon = addparameter(model,'kon_CTLA4_aCTLA4',params.kon_CTLA4_aCTLA4.Value,'ValueUnits',params.kon_CTLA4_aCTLA4.Units);
    set(kon,'Notes',['kon of CTLA4-aCTLA4 binding ' params.kon_CTLA4_aCTLA4.Notes]);
catch
    first_call = false;
end

if first_call
    % add aCTLA-4 PK if it does not exist yet
    params_aCTLA4  = pk_parameters('ipilimumab');
    model = pk_module(model,'aCTLA4' ,params_aCTLA4);
    % Add koff Values
    koff = addparameter(model,'koff_CTLA4_aCTLA4',params.koff_CTLA4_aCTLA4.Value,'ValueUnits',params.koff_CTLA4_aCTLA4.Units);
        set(koff,'Notes',['koff of CTLA4-aCTLA4 binding ' params.koff_CTLA4_aCTLA4.Notes]);
    % Bivalent anibody parameters
    p = addparameter(model,'Chi_CTLA4_aCTLA4' ,params.Chi_CTLA4_aCTLA4.Value ,'ValueUnits',params.Chi_CTLA4_aCTLA4.Units);
        set(p,'Notes',['Antibody cross-arm binding efficiency ' params.Chi_CTLA4_aCTLA4.Notes]);
end

% Add area of a T cell if not defined before in antigen module
try
    parameter = addparameter(model,'A_Tcell' ,params.A_Tcell.Value ,'ValueUnits',params.A_Tcell.Units);
    set(parameter,'Notes',['Surface area of the T cell ' params.A_Tcell.Notes]);
catch
end

% Total amount of CTLA4 on Treg
p = addparameter(model,'Treg_CTLA4_tot',params.Treg_CTLA4_tot.Value,'ValueUnits',params.Treg_CTLA4_tot.Units,'ConstantValue',false);
    set(p,'Notes',['Total number of CTLA4 on Treg cells ' params.Treg_CTLA4_tot.Notes]);
% CTLA4-related Hill parameters
p = addparameter(model,'Treg_CTLA4_50',params.Treg_CTLA4_50.Value,'ValueUnits',params.Treg_CTLA4_50.Units);
    set(p,'Notes',['CTLA4 occupancy for half-maximal Treg inactivation by macrophages ' params.Treg_CTLA4_50.Notes]);
p = addparameter(model,'n_Treg_CTLA4',params.n_Treg_CTLA4.Value,'ValueUnits',params.n_Treg_CTLA4.Units);
    set(p,'Notes',['CTLA4 occupancy Hill coefficient for Treg inactivation by macrophages ' params.n_Treg_CTLA4.Notes]);
% Treg ADCC rates
p = addparameter(model,'k_CTLA4_ADCC',params.k_CTLA4_ADCC.Value,'ValueUnits',params.k_CTLA4_ADCC.Units,'ConstantValue',false);
    set(p,'Notes',[params.k_CTLA4_ADCC.Notes]);
% Define Hill function for ADCC in tumor and peripheral compartments
p = addparameter(model,'H_Treg_P',0.0,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of Treg ADCC in peripheral compartment');

for array=1:ntumors

    p = addparameter(model,['H_Treg_T' TumorArray{array}],0.0,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of Treg ADCC in tumor');
    
    % Species for states of CTLA4 on Treg
    x = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'Treg_CTLA4',0,'InitialAmountUnits','molecule'); 
        set(x,'Notes','Number of free CTLA4 molecules on Treg in Tumour compartment');
    x = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'Treg_CTLA4_aCTLA4',0,'InitialAmountUnits','molecule'); 
        set(x,'Notes','Number of CTLA4-aCTLA4 complex on Treg in Tumour compartment');
    x = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'Treg_CTLA4_aCTLA4_CTLA4',0,'InitialAmountUnits','molecule'); 
        set(x,'Notes','Number of CTLA4-aCTLA4-CTLA4 complex on Treg in Tumour compartment');

end

x = addspecies(sbioselect(model, 'Name', 'V_P'),'Treg_CTLA4',0,'InitialAmountUnits','molecule'); 
    set(x,'Notes','Number of free CTLA4 molecules on Treg in Peripheral compartment');
x = addspecies(sbioselect(model, 'Name', 'V_P'),'Treg_CTLA4_aCTLA4',0,'InitialAmountUnits','molecule'); 
    set(x,'Notes','Number of CTLA4-aCTLA4 complex on Treg in Peripheral compartment');
x = addspecies(sbioselect(model, 'Name', 'V_P'),'Treg_CTLA4_aCTLA4_CTLA4',0,'InitialAmountUnits','molecule'); 
    set(x,'Notes','Number of CTLA4-aCTLA4-CTLA4 complex on Treg in Peripheral compartment');

% Initialize the total CTLA4 on the cells
addrule(model,'V_P.Treg_CTLA4 = Treg_CTLA4_tot' ,'initialAssignment');
% Add the Hill parameters and the rules for them
addrule(model,'H_Treg_P = ((V_P.Treg_CTLA4_aCTLA4+2*V_P.Treg_CTLA4_aCTLA4_CTLA4)/Treg_CTLA4_50)^n_Treg_CTLA4/(((V_P.Treg_CTLA4_aCTLA4+2*V_P.Treg_CTLA4_aCTLA4_CTLA4)/Treg_CTLA4_50)^n_Treg_CTLA4 + 1)','repeatedAssignment');

for array=1:ntumors

    addrule(model,['V_T' TumorArray{array} '.Treg_CTLA4 = Treg_CTLA4_tot'],'initialAssignment');
    addrule(model,['H_Treg_T' TumorArray{array} ' = ((V_T' TumorArray{array} '.Treg_CTLA4_aCTLA4+2*V_T' TumorArray{array} '.Treg_CTLA4_aCTLA4_CTLA4)/Treg_CTLA4_50)^n_Treg_CTLA4/(((V_T' TumorArray{array} '.Treg_CTLA4_aCTLA4+2*V_T' TumorArray{array} '.Treg_CTLA4_aCTLA4_CTLA4)/Treg_CTLA4_50)^n_Treg_CTLA4 + 1)'],'repeatedAssignment');

    % Binding and unbinding of aCTLA4 to CTLA4 on Treg on Tumour and Peripheral
    R = addreaction(model, ['V_T' TumorArray{array} '.Treg_CTLA4 <-> V_T' TumorArray{array} '.Treg_CTLA4_aCTLA4']);
        set (R, 'ReactionRate', ['kon_CTLA4_aCTLA4*(V_T' TumorArray{array} '.Treg_CTLA4 * V_T' TumorArray{array} '.aCTLA4/gamma_T_aCTLA4) -  koff_CTLA4_aCTLA4*V_T' TumorArray{array} '.Treg_CTLA4_aCTLA4']);
        set (R, 'Notes'       , 'binding and unbinding of CTLA4 to aCTLA4 on Treg surface in Tumor');
    R = addreaction(model, ['V_T' TumorArray{array} '.Treg_CTLA4_aCTLA4 + V_T' TumorArray{array} '.Treg_CTLA4 <-> V_T' TumorArray{array} '.Treg_CTLA4_aCTLA4_CTLA4']);
        set (R, 'ReactionRate', ['Chi_CTLA4_aCTLA4*kon_CTLA4_aCTLA4*(V_T' TumorArray{array} '.Treg_CTLA4 * V_T' TumorArray{array} '.Treg_CTLA4_aCTLA4)/A_Tcell -  koff_CTLA4_aCTLA4*V_T' TumorArray{array} '.Treg_CTLA4_aCTLA4_CTLA4']);
        set (R, 'Notes'       , 'binding and unbinding of CTLA4 to CTLA4-aCTLA4 on Treg surface in Tumor');
    
    % Treg Death through CTLA4 binding
    R = addreaction(model,['V_T' TumorArray{array} '.T0 -> null']);
        set(R,'ReactionRate',['k_CTLA4_ADCC*V_T' TumorArray{array} '.T0*H_Treg_T' TumorArray{array}]);
        set(R,'Notes','CTLA4 ADCC in the tumor compartment');

end

R = addreaction(model, 'V_P.Treg_CTLA4 <-> V_P.Treg_CTLA4_aCTLA4');
    set (R, 'ReactionRate', 'kon_CTLA4_aCTLA4*(V_P.Treg_CTLA4 * V_P.aCTLA4/gamma_P_aCTLA4) -  koff_CTLA4_aCTLA4*V_P.Treg_CTLA4_aCTLA4');
    set (R, 'Notes'       , 'binding and unbinding of CTLA4 to aCTLA4 on Treg surface in Peripheral compartment');
R = addreaction(model, 'V_P.Treg_CTLA4_aCTLA4 + V_P.Treg_CTLA4 <-> V_P.Treg_CTLA4_aCTLA4_CTLA4');
    set (R, 'ReactionRate', 'Chi_CTLA4_aCTLA4*kon_CTLA4_aCTLA4*(V_P.Treg_CTLA4 * V_P.Treg_CTLA4_aCTLA4)/A_Tcell -  koff_CTLA4_aCTLA4*V_P.Treg_CTLA4_aCTLA4_CTLA4');
    set (R, 'Notes'       , 'binding and unbinding of CTLA4 to CTLA4-aCTLA4 on Treg surface in Peripheral compartment');
R = addreaction(model,'V_P.T0 -> null');
    set(R,'ReactionRate','k_CTLA4_ADCC*V_P.T0*H_Treg_P');
    set(R,'Notes','CTLA4 ADCC in the peripheral compartment');
