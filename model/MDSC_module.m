% MDSC Module
%
% Inputs: model        -- simbio model object with four compartments
%         species      - MDSC--Myeloid-derived suppressor cell
%                      - ArgI--Arginase I
%                      - CCL2--Monocyte Chemoattractant Protein-1 (CCL2)
%                      - NO--Nitrite Oxide
%         params
%                      - k_MDSC_mig--MDSC recruitment rate by MCP-1
%                      - k_MDSC_death--MDSC death rate
%                      - k_CCL2_deg--Degradation rate of MCP-1
%                      - k_NO_deg--Degradation rate of NO
%                      - k_ArgI_deg--Degradation rate of Arg I
%                      - k_CCL2_sec--Secretion rate of MCP-1
%                      - k_NO_sec--Secretion rate of NO
%                      - k_ArgI_sec--Secretion rate of Arg I
%                      - IC50_inostat_C--Effective concentration of -inostat on inhibition of cancer proliferation
%                      - IC50_inostat_NO--Effective concentration of -inostat on inhibition of NO production
%                      - IC50_inostat_ArgI--Effective concentration of -inostat on inhibition of Arg I production
%                      - IC50_inostat_CCL2--Effective concentration of -inostat on inhibition of MCP-1 production
%                      - ArgI_50_Teff--Effective concentration of Arg I on inhibition of Teff activity
%                      - NO_50_Teff--Effective concentration of NO on inhibition of Teff activity
%                      - CCL2_50--Effective concentration of MCP-1 on recruitment of MDSC
%                      - ArgI_50_Treg--Effective concentration of Arg I on Th-to-Treg transdifferentiation
%        TumorArray    -- Array of suffixes for tumor comp. names
%        LNArray       -- Array of suffixes for LN comp. names
% Outputs: model -- SimBiology model object with new MDSC module


function model = MDSC_module(model,params,TumorArray,LNArray,varargin)

% number of tumors and LNs
ntumors = length(TumorArray);
uniqLNArray = unique(LNArray);
nLNs = length(uniqLNArray);

% Optional Inputs
in = inputParser;
addOptional(in,'cancer_types',{'C1'});
addParameter(in,'inostat',1);
addParameter(in,'drugName','entinostat');
% Parse Inputs
parse(in,varargin{:});
cancer_types = in.Results.cancer_types;
inostat = in.Results.inostat;
drugName = in.Results.drugName;

% Species Names
species_name = 'MDSC';

for array=1:ntumors
    % Add Species
    % MDSC in tumor
    MDSC = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'MDSC',0,'InitialAmountUnits','cell'); 
        set(MDSC,'Notes','Number of MDSCs in the tumour compartment');
    % add NO
    NO = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'NO',0,'InitialAmountUnits','nanomolarity'); 
        set(NO,'Notes','Concentration of NO in the tumor compartment');
    % add ArgI
    ArgI = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'ArgI',0,'InitialAmountUnits','mU'); 
        set(ArgI,'Notes','Concentration of Arg I in the tumor compartment');
    % CCL2 in tumor
    CCL2 = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'CCL2',0,'InitialAmountUnits','nanomolarity'); 
        set(CCL2,'Notes','Concentration of CCL2 in the tumor compartment');
end

first_call = true;
try
    % CCL2 Parameters
    p = addparameter(model,'k_CCL2_sec',params.k_CCL2_sec.Value,'ValueUnits',params.k_CCL2_sec.Units);
        set(p,'Notes',['rate of CCL2 secretion ' params.k_CCL2_sec.Notes]);
    p = addparameter(model,'k_CCL2_deg',params.k_CCL2_deg.Value,'ValueUnits',params.k_CCL2_deg.Units);
        set(p,'Notes',['rate of CCL2 degradation ' params.k_CCL2_deg.Notes]);
    p = addparameter(model,'CCL2_50',params.CCL2_50.Value,'ValueUnits',params.CCL2_50.Units);
        set(p,'Notes',['Half-maximal CCL2 level of MDSC recruitment ' params.CCL2_50.Notes]);
catch
    first_call = false;
end

% Add Parameters
for array=1:ntumors
    k_MDSC_mig = addparameter(model,['k_MDSC_mig' TumorArray{array}],params.(['k_MDSC_mig' TumorArray{array}]).Value,'ValueUnits',params.(['k_MDSC_mig' TumorArray{array}]).Units);
        set(k_MDSC_mig,'Notes',['Rate of MDSC recruitment into the tumor ' params.(['k_MDSC_mig' TumorArray{array}]).Notes]);
end

k_MDSC_death = addparameter(model,'k_MDSC_death',params.k_MDSC_death.Value,'ValueUnits',params.k_MDSC_death.Units);
    set(k_MDSC_death,'Notes',['Rate of MDSC death ' params.k_MDSC_death.Notes]);
p = addparameter(model,'k_NO_deg',params.k_NO_deg.Value,'ValueUnits',params.k_NO_deg.Units);
    set(p,'Notes',['rate of NO degradation ' params.k_NO_deg.Notes]);
p = addparameter(model,'k_ArgI_deg',params.k_ArgI_deg.Value,'ValueUnits',params.k_ArgI_deg.Units);
    set(p,'Notes',['rate of ArgI degradation ' params.k_ArgI_deg.Notes]);
p = addparameter(model,'k_NO_sec',params.k_NO_sec.Value,'ValueUnits',params.k_NO_sec.Units);
    set(p,'Notes',['rate of NO secretion from MDSCs ' params.k_NO_sec.Notes]);
p = addparameter(model,'k_ArgI_sec',params.k_ArgI_sec.Value,'ValueUnits',params.k_ArgI_sec.Units);
    set(p,'Notes',['rate of ArgI secretion ' params.k_ArgI_sec.Notes]);
p = addparameter(model,'ArgI_50_Teff',params.ArgI_50_Teff.Value,'ValueUnits',params.ArgI_50_Teff.Units);
    set(p,'Notes',['rate of ArgI-induced T cell death ' params.ArgI_50_Teff.Notes]);
p = addparameter(model,'NO_50_Teff',params.NO_50_Teff.Value,'ValueUnits',params.NO_50_Teff.Units);
    set(p,'Notes',['rate of NO-induced T cell death ' params.NO_50_Teff.Notes]);
p = addparameter(model,'ArgI_50_Treg',params.ArgI_50_Treg.Value,'ValueUnits',params.ArgI_50_Treg.Units);
    set(p,'Notes',['Half-maximal ArgI level of Treg expansion ' params.ArgI_50_Treg.Notes]);

for array=1:ntumors
    
    % Add Reactions
    % Recruitment of MDSC (Huang 2006, PMID: 17257744)
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.MDSC']);
        set(reaction,'ReactionRate',['k_MDSC_mig' TumorArray{array} '*V_T' TumorArray{array} '*(V_T' TumorArray{array} '.CCL2/(V_T' TumorArray{array} '.CCL2 + CCL2_50))']);
        set(reaction,'Notes','Recruitment of macrophage to tumor');
    % MDSC Death
    reaction = addreaction(model,['V_T' TumorArray{array} '.MDSC -> null']);
        set(reaction,'ReactionRate',['k_MDSC_death*V_T' TumorArray{array} '.MDSC']);
        set(reaction,'Notes','MDSC death in the tumor compartment');
    % MDSC death upon tumor eradication
    reaction = addreaction(model,['V_T' TumorArray{array} '.MDSC -> null']);
        set(reaction,'ReactionRate',['k_cell_clear*V_T' TumorArray{array} '.MDSC*(Kc_rec/(C_total' TumorArray{array} '^2 + Kc_rec))']);
        set(reaction,'Notes','MDSC death upon tumor eradication');
    % CCL2 secretion by MDSCs and Cancer cells
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.CCL2']);
        set(reaction,'ReactionRate',['k_CCL2_sec*C_total' TumorArray{array}]);
        set(reaction,'Notes','CCL2 secretion by MDSCs and cancer cells');
    % CCL2 Degradation
    reaction = addreaction(model,['V_T' TumorArray{array} '.CCL2 -> null']);
        set(reaction,'ReactionRate',['k_CCL2_deg*V_T' TumorArray{array} '.CCL2']);
        set(reaction,'Notes','CCL2 degradation');
    % NO Degradation
    reaction = addreaction(model,['V_T' TumorArray{array} '.NO -> null']);
        set(reaction,'ReactionRate',['k_NO_deg*V_T' TumorArray{array} '.NO']);
        set(reaction,'Notes','NO degradation');
    % ArgI Degradation
    reaction = addreaction(model,['V_T' TumorArray{array} '.ArgI -> null']);
        set(reaction,'ReactionRate',['k_ArgI_deg*V_T' TumorArray{array} '.ArgI']);
        set(reaction,'Notes','ArgI degradation');
    % NO Secretion by MDSC
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.NO']);
        set(reaction,'ReactionRate',['k_NO_sec*V_T' TumorArray{array} '.MDSC']);
        set(reaction,'Notes','NO secretion from MDSC');
    % ArgI Secretion by MDSC
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.ArgI']);
        set(reaction,'ReactionRate',['k_ArgI_sec*V_T' TumorArray{array} '.MDSC']);        
        set(reaction,'Notes','ArgI secretion from MDSC');
    
    % Set Default Hill Function for MDSC and -inostat
    H_NO = addparameter(model,['H_NO' TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
        set(H_NO,'Notes','Hill function for NO-mediated inhibition on Teff');
    H_ArgI_Teff = addparameter(model,['H_ArgI_Teff' TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
        set(H_ArgI_Teff,'Notes','Hill function for ariginase I-mediated inhibition on Teff');
    H_ArgI_Treg = addparameter(model,['H_ArgI_Treg' TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
        set(H_ArgI_Treg,'Notes','Hill function for ariginase I-mediated Th-to-Treg transdifferentiation');
    H_MDSC = addparameter(model,['H_MDSC' TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
        set(H_MDSC,'Notes','Hill function for overall effect of MDSC-mediated effect on Teff');
    
    % Add Hill Functions
    addrule(model,['H_NO' TumorArray{array} ' = V_T' TumorArray{array} '.NO/(NO_50_Teff+V_T' TumorArray{array} '.NO)'],'repeatedAssignment');
    addrule(model,['H_ArgI_Teff' TumorArray{array} ' = V_T' TumorArray{array} '.ArgI/(ArgI_50_Teff+V_T' TumorArray{array} '.ArgI)'],'repeatedAssignment');
    addrule(model,['H_ArgI_Treg' TumorArray{array} ' = V_T' TumorArray{array} '.ArgI/(ArgI_50_Treg+V_T' TumorArray{array} '.ArgI)'],'repeatedAssignment');
    addrule(model,['H_MDSC' TumorArray{array} ' = 1-(1-H_NO' TumorArray{array} ')*(1-H_ArgI_Teff' TumorArray{array} ')'],'repeatedAssignment');
end 
    
    model_rules = get(model,'Rules');

% Add -inostat drug mechanisms
if (inostat)
    p = addparameter(model,'IC50_inostat_C',params.IC50_inostat_C.Value,'ValueUnits',params.IC50_inostat_C.Units);
        set(p,'Notes',['-inostat concentration for half-maximal tumor cell death ' params.IC50_inostat_C.Notes]);
    p = addparameter(model,'IC50_inostat_NO',params.IC50_inostat_NO.Value,'ValueUnits',params.IC50_inostat_NO.Units);
        set(p,'Notes',['half-maximal -inostat concentration for NO inhibition ' params.IC50_inostat_NO.Notes]);
    p = addparameter(model,'IC50_inostat_CCL2',params.IC50_inostat_CCL2.Value,'ValueUnits',params.IC50_inostat_CCL2.Units);
        set(p,'Notes',['half-maximal -inostat concentration for CCL2 inhibition ' params.IC50_inostat_CCL2.Notes]);
    p = addparameter(model,'IC50_inostat_ArgI',params.IC50_inostat_ArgI.Value,'ValueUnits',params.IC50_inostat_ArgI.Units);
        set(p,'Notes',['half-maximal -inostat concentration for Arg I inhibition ' params.IC50_inostat_ArgI.Notes]);

    for array=1:ntumors

        H_inostat_C = addparameter(model,['H_inostat_C' TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
            set(H_inostat_C,'Notes','Hill function for -inostat-mediated inhibition on cancer cell growth');
        H_inostat_ArgI = addparameter(model,['H_inostat_ArgI' TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
            set(H_inostat_ArgI,'Notes','Hill function for -inostat-mediated inhibition on arginase I activity ');
    
        addrule(model,['H_inostat_C' TumorArray{array} ' = V_T' TumorArray{array} '.inostat/(V_T' TumorArray{array} '.inostat+IC50_inostat_C)'],'repeatedAssignment');
        addrule(model,['H_inostat_ArgI' TumorArray{array} ' = V_T' TumorArray{array} '.inostat/(V_T' TumorArray{array} '.inostat+IC50_inostat_ArgI)'],'repeatedAssignment');

        % Add -inostat drug-mediated effects
        for i = 1:length(model.reaction)
            if strcmp(model.reaction(i).reaction, ['null -> V_T' TumorArray{array} '.CCL2']) && isempty(strfind(model.reaction(i).reaction, 'IC50_inostat_CCL2'))
                model.reaction(i).ReactionRate = [model.reaction(i).ReactionRate, '*(1-V_T' TumorArray{array} '.inostat/(V_T' TumorArray{array} '.inostat+IC50_inostat_CCL2))'];
            elseif strcmp(model.reaction(i).reaction, ['null -> V_T' TumorArray{array} '.NO']) && isempty(strfind(model.reaction(i).reaction, 'IC50_inostat_NO'))
                model.reaction(i).ReactionRate = [model.reaction(i).ReactionRate, '*(1-V_T' TumorArray{array} '.inostat/(V_T' TumorArray{array} '.inostat+IC50_inostat_NO))'];
            elseif strcmp(model.reaction(i).reaction, ['null -> V_T' TumorArray{array} '.ArgI']) && isempty(strfind(model.reaction(i).reaction, 'IC50_inostat_ArgI'))
                model.reaction(i).ReactionRate = [model.reaction(i).ReactionRate, '*(1-V_T' TumorArray{array} '.inostat/(V_T' TumorArray{array} '.inostat+IC50_inostat_ArgI))'];
            end
            for j = 1:length(cancer_types)
                if strcmp(model.reaction(i).reaction, ['null -> V_T' TumorArray{array} '.' cancer_types{j}]) && isempty(strfind(model.reaction(i).reaction, ['*(1-H_inostat_C' TumorArray{array} ')']))
                    model.reaction(i).ReactionRate = [model.reaction(i).ReactionRate, '*(1-H_inostat_C' TumorArray{array} ')'];
                end
            end
        end
    
        % Update Hill function for MDSC-mediated inhibition
        for i = 1:length(model_rules)
            if ~isempty(strfind(model_rules(i).Rule, ['H_MDSC' TumorArray{array} ' = '])) && isempty(strfind(model.reaction(i).reaction, ['H_inostat_ArgI' TumorArray{array}]))
                model_rules(i).Rule = insertAfter(model_rules(i).Rule, ['H_ArgI_Teff' TumorArray{array}], ['*(1-H_inostat_ArgI)' TumorArray{array}]);
            end
        end

    end

    % Add Pharmacokinetics of -inostat
    params_inostat    = pk_parameters(drugName);
    model = pk_module(model,'inostat',TumorArray,LNArray,params_inostat,'o');
end

for array=1:ntumors
    %model_rules = get(model,'Rules'); %thein
    ruleobj=sbioselect(model,'Where','Rule', 'regexp', ['^R_Tcell' TumorArray{array} ' = ']); 
    for i = 1:length(model.reaction)
        for j = 1:length(cancer_types)

            % Update Teff-mediated killing rate
            if strcmp(model.reaction(i).reaction, ['V_T' TumorArray{array} '.' cancer_types{j} ' -> V_T' TumorArray{array} '.C_x']) ...
                    && isempty(strfind(model.reaction(i).reaction, ['*(1-H_MDSC' TumorArray{array} ')'])) && ~isempty(strfind(model.reaction(i).ReactionRate, 'k_C_T'))
                model.reaction(i).ReactionRate = [model.reaction(i).ReactionRate, '*(1-H_MDSC' TumorArray{array} ')'];
                ruleobj.Rule = insertAfter(ruleobj.Rule, ['*(1-H_PD1_' cancer_types{j} TumorArray{array} ')'], ['*(1-H_MDSC' TumorArray{array} ')']);
            % Update antigen release rate
            elseif ~isempty(strfind(model.reaction(i).reaction, ['null -> V_T' TumorArray{array} '.P']))
                model.reaction(i).ReactionRate = insertAfter(model.reaction(i).ReactionRate, ['*(1-H_PD1_' cancer_types{j} TumorArray{array} ')'], ['*(1-H_MDSC' TumorArray{array} ')']);
            end

        end

         % Update Th-to-Treg transdifferentiation
        if strcmp(model.reaction(i).reaction, ['V_T' TumorArray{array} '.Th -> V_T' TumorArray{array} '.T0'])
            model.reaction(i).ReactionRate = [model.reaction(i).ReactionRate, '*H_ArgI_Treg' TumorArray{array}];
        end

    end

end
