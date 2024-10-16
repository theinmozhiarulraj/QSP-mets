% Function to initialize simbiology object
%
% Inputs: model_name -- string containing model name
%         time       -- array of time points (in days)
%         solver     -- string containing solver (e.g. ode15s)
%         tol_abs    -- absolute tolerance of solver
%         tol_rel    -- relative tolerance of solver
%         params     -- object containing compartment parameters
%                       - V_C--central compartment volume
%                       - V_P--peripheral compartment volume
%                       - V_Tmin--cancer-free tumour compartment volume
%                       - V_LN--lymph node compartment volume
%                       - vol_cell--cancer cell volume
%                       - vol_Tcell--T cell volume
%                       - k_cell_clear--dead cell clearance rate
%       ncancer_clones -- Number of cancer clones
%       n_T_specs      -- number of T cell specificies
%   ag_fraction_clones -- fraction of different antigens in cancer clones
%       TumorArray     -- Array of suffixes for tumor comp. names
%       LNArray        -- Array of suffixes for LN comp. names
% Outputs: model -- simbiology model object containing user-specified
%                   configuration and four compartments: C,P,T,LN


function model = simbio_init(model_name,time,solver,tol_abs,tol_rel,params,ncancer_clones,n_T_specs,TumorArray,LNArray)

% % number of cancer clones and neoantigens
% [ncancer_clones nAgs]=size(ag_fraction_clones);

% number of tumors and LNs
ntumors = length(TumorArray);
uniqLNArray = unique(LNArray);
nLNs = length(uniqLNArray);

% Model Object
model = sbiomodel(model_name);
config = getconfigset(model);
options = get(config,'CompileOptions');
set(options,'UnitConversion',true);
set(config,'TimeUnits','day');
set(config,'SolverType',solver);
set(config.SolverOptions,'OutputTimes',time);
% set(config.SolverOptions,'MaxStep',0.01);
set(config,'StopTime',time(end));
set(config.SolverOptions,'AbsoluteTolerance',tol_abs);
set(config.SolverOptions,'RelativeTolerance',tol_rel);

% Setup Compartments
comp_C = addcompartment(model,'V_C',params.V_C.Value,'CapacityUnits',params.V_C.Units);
    set(comp_C,'Notes',['Central compartment (C) ' params.V_C.Notes]);
comp_P = addcompartment(model,'V_P',params.V_P.Value,'CapacityUnits',params.V_P.Units);
    set (comp_P,'Notes',['Peripheral compartment (P) ' params.V_P.Notes]);

for array=1:ntumors

    genvarname(['comp_T' TumorArray{array}]) = addcompartment(model,['V_T' TumorArray{array}],params.V_Tmin.Value,'CapacityUnits',params.V_Tmin.Units,'ConstantCapacity',false);
    set (genvarname(['comp_T' TumorArray{array}]),'Notes',['Tumor compartment (T) '  params.V_Tmin.Notes]);

end

for array=1:nLNs

    genvarname(['comp_LN' uniqLNArray{array}]) = addcompartment(model,['V_LN' uniqLNArray{array}],eval(['params.V_LN' uniqLNArray{array} '.Value']),'CapacityUnits',eval(['params.V_LN' uniqLNArray{array} '.Units']));
    set(genvarname(['comp_LN' uniqLNArray{array}]),'Notes',['Lymph node (LN) compartment volume ' eval(['params.V_LN' uniqLNArray{array} '.Notes']) ]);

end

% Define Cell and Time
p = addparameter(model,'cell',1.0,'ValueUnits','cell');
    set(p,'Notes','unit parameter for calculation');
p = addparameter(model,'day',1.0,'ValueUnits','day');
    set(p,'Notes','unit parameter for calculation');
p = addparameter(model,'V_Tmin',params.V_Tmin.Value,'ValueUnits',params.V_Tmin.Units);
    set(p,'Notes','minimum tumor compartment volume');
vol_cell = addparameter(model,'vol_cell',params.vol_cell.Value,'ValueUnits',params.vol_cell.Units);
    set(vol_cell,'Notes',['Average volume of cancer cell ' params.vol_cell.Notes]);
vol_Tcell = addparameter(model,'vol_Tcell',params.vol_Tcell.Value,'ValueUnits',params.vol_Tcell.Units);
    set(vol_Tcell,'Notes',['Average volume of T cells ' params.vol_Tcell.Notes]);

% Dead Cells
P = addparameter(model,'k_cell_clear',params.k_cell_clear.Value,'ValueUnits',params.k_cell_clear.Units);
    set(P,'Notes',['Rate of dead cell clearance from tumour compartment ' params.k_cell_clear.Notes]);
Ve_T = addparameter(model,'Ve_T',params.Ve_T.Value,'ValueUnits',params.Ve_T.Units,'ConstantValue',false);
    set(Ve_T,'Notes',['Void fraction of the tumor ' params.Ve_T.Notes]);

for nc=1:ncancer_clones

    %disp(nc)
    
    sum_agconc=0;
    
    % get proportions of all Ags 
    for nt=1:(n_T_specs+1)
        agconc = addparameter(model,['agconc_C' num2str(nc) '_' num2str(nt-1)],params.(['agconc_' 'C' num2str(nc) '_' num2str(nt-1)]).Value,'ValueUnits',params.(['agconc_' 'C' num2str(nc) '_' num2str(nt-1)]).Units,'ConstantValue',false);
            set(agconc,'Notes','Concentration of antigens');
    
        %rule = get(totalag_rule,'Rule');
        %set(totalag_rule,'Rule',[rule '+agconc_C' num2str(nc) '_' num2str(nt-1)]);

        if(nt>1) % don't add self ag conc
            sum_agconc=sum_agconc+agconc.Value;
        end
    end

    totalag = addparameter(model,['totalag_' 'C' num2str(nc)],sum_agconc,'ValueUnits','mole/cell','ConstantValue',false);
            set(totalag,'Notes','Total neoantigen concentration in a cancer cell clone');

end

for array=1:ntumors

    compname_tumor=['V_T' TumorArray{array}];

    S = addspecies(sbioselect(model, 'Name', compname_tumor),'C_x',0,'InitialAmountUnits','cell');
        set(S,'Notes','Dead cancer cells in the tumour compartment');
    S = addspecies(sbioselect(model, 'Name', compname_tumor),'T1_exh',0,'InitialAmountUnits','cell');
        set(S,'Notes','Exhausted CD8 effector T cells');
    S = addspecies(sbioselect(model, 'Name', compname_tumor),'Th_exh',0,'InitialAmountUnits','cell');
        set(S,'Notes','Exhausted CD4 effector T cells');
    R = addreaction(model,['V_T' TumorArray{array} '.C_x  -> null']);
        set(R,'ReactionRate',['k_cell_clear*V_T' TumorArray{array} '.C_x']);
        set(R,'Notes','Clearance of dead cancer cells from tumor');
    R = addreaction(model,['V_T' TumorArray{array} '.T1_exh -> null']);
        set(R,'ReactionRate',['k_cell_clear*V_T' TumorArray{array} '.T1_exh']);
        set(R,'Notes','Clearance of exhausted CD8 T cells from tumor');
    R = addreaction(model,['V_T' TumorArray{array} '.Th_exh -> null']);
        set(R,'ReactionRate',['k_cell_clear*V_T' TumorArray{array} '.Th_exh']);
        set(R,'Notes','Clearance of exhausted CD4 T cells from tumor');
    
    % Set Tumour Volume (Rule 1)
    % addrule(model,'V_T = (C_x+C_total)*vol_cell/Ve_T','repeatedAssignment');
    addrule(model,['V_T' TumorArray{array} '= V_Tmin+((V_T' TumorArray{array} '.C_x+C_total' TumorArray{array} ')*vol_cell+(V_T' TumorArray{array} '.T1_exh+V_T' TumorArray{array} '.Th_exh+T_total' TumorArray{array} ')*vol_Tcell)/Ve_T'],'repeatedAssignment'); % total cells in tumor / constant cellular volume fraction
    %addrule(model,['V_T' TumorArray{array} '= ((V_T' TumorArray{array} '.C_x+C_total' TumorArray{array} ')*vol_cell+(V_T' TumorArray{array} '.T1_exh+V_T' TumorArray{array} '.Th_exh+T_total' TumorArray{array} ')*vol_Tcell)/Ve_T'],'repeatedAssignment'); % total cells in tumor / constant cellular volume fraction

    % Set Total Number of Cancer Cells (Rule 2)
    p = addparameter(model,['C_total' TumorArray{array}],0,'ValueUnits','cell','ConstantValue',false);
        set(p,'Notes','Total number of cancer cells');
    addrule(model,['C_total' TumorArray{array} '= 0*cell'],'repeatedAssignment');
    
    % Set Total Number of T Cells in Tumour (Rule 3)
    p = addparameter(model,['T_total' TumorArray{array}],0,'ValueUnits','cell','ConstantValue',false);
        set(p,'Notes','Total number of activated T cells in tumor');
    addrule(model,['T_total' TumorArray{array} ' = 0*cell'],'repeatedAssignment');

    p = addparameter(model,['Tcyt_total' TumorArray{array}],0,'ValueUnits','cell','ConstantValue',false);
        set(p,'Notes','Total number of activated T cells in tumor');
    addrule(model,['Tcyt_total' TumorArray{array} ' = 0*cell'],'repeatedAssignment');

    for nc=1:ncancer_clones
        p = addparameter(model,['C' num2str(nc) '_Tcyt_total' TumorArray{array}],0,'ValueUnits','cell','ConstantValue',false);
            set(p,'Notes','Total number of T cells in tumor that can recognize cancer clone');
        cancer_spec_rule=addrule(model,['C' num2str(nc) '_Tcyt_total' TumorArray{array} ' = 0*cell'],'repeatedAssignment');

        total_conc=sbioselect(model, 'Name', ['totalag_C' num2str(nc)]).Value;

        % get proportions of all Ags 
        for nt=1:n_T_specs

            conc=sbioselect(model, 'Name', ['agconc_C' num2str(nc) '_' num2str(nt)]).Value;

            rule = get(cancer_spec_rule,'Rule');
            %disp(conc/total_conc)
            if((conc/total_conc)>=0.05)
                set(cancer_spec_rule,'Rule',[rule '+V_T' TumorArray{array} '.T' num2str(nt)]);
            end
        end

    end

    % Set Total Rate of Cancer Death by T Cells (Rule 5)
    p = addparameter(model,['R_Tcell' TumorArray{array}],'ValueUnits','cell/day','ConstantValue',false);
        set(p,'Notes','Rate of cancer cell death');
    addrule(model,['R_Tcell' TumorArray{array} ' = 0*cell/day'],'repeatedAssignment');

    % Set Default Number of Tregs
    p = addparameter(model,['Tregs_' TumorArray{array}],0,'ValueUnits','cell','ConstantValue',false);
        %set(p,'Notes','Total number of activated T cells in TDLNs');
        set(p,'Notes','Total number of Tregs in tumor');

    % Set Default Hill Function for PD1 Checkpoint
    p = addparameter(model,['H_PD1_C1' TumorArray{array}],0.90,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes','Hill function of PDL1 on tumor cells for T cell exhaustion');
    
    % Set Default Hill Function for CTLA4 Checkpoint
    p = addparameter(model,['H_CD28_C1' TumorArray{array}],0.1,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes','Hill function of CD28 on tumor cells');

    % Set Default Hill Function for APCs
    p = addparameter(model,['H_APC' TumorArray{array}],0.5,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes','Hill function of mAPC for self-antigen-specific Treg activation');
    
    % Set Default Hill Function for mAPCs
    p = addparameter(model,['H_mAPC' TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes','Hill function of mAPC for CD8 T cell activation');
    
    % Set Default Hill Function for mAPCs
    p = addparameter(model,['H_APCh' TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes','Hill function of mAPC for neoantigen-specific helper T cell activation');    

end

for array=1:nLNs

    % Set Total Number of T Cells in LN (Rule 4)
    p = addparameter(model,['T_total_LN' uniqLNArray{array}],0,'ValueUnits','cell','ConstantValue',false);
        set(p,'Notes','Total number of activated T cells in TDLNs');
    addrule(model,['T_total_LN' uniqLNArray{array} ' = 0*cell'],'repeatedAssignment');

    % addparameter(model,'H_PD1_C2',0.90,'ValueUnits','dimensionless','ConstantValue',false);
    p = addparameter(model,['H_PD1_APC' uniqLNArray{array}],0.90,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes','Hill function of PDL1 on APC for T cell exhaustion');
    
    % addparameter(model,'H_CD28_C2',0.1,'ValueUnits','dimensionless','ConstantValue',false);
    p = addparameter(model,['H_CD28_APC' uniqLNArray{array}],0.1,'ValueUnits','dimensionless','ConstantValue',false);
        set(p,'Notes','Hill function of CD28 on APCs for T cell activation');
    
end

% This is for metastatic tumor seeding

% delay can be -1 -> no delay
% delay can also be greater than 0
% delay = 0 -> might not work

addparameter(model,'seeding',1,'ValueUnits','dimensionless','ConstantValue',false);
% primary tumor, start should be always 1
addparameter(model,'start',1,'ValueUnits','dimensionless','ConstantValue',false);

for array=2:ntumors

    addparameter(model,['delay' TumorArray{array}],params.(['delay' TumorArray{array}]).Value,'ValueUnits','day');
    % if delay =-1, start should be 1, if delay is greater than 0 - set as 0
    if(params.(['delay' TumorArray{array}]).Value<=0)
        addparameter(model,['start' TumorArray{array}],1,'ValueUnits','dimensionless','ConstantValue',false);
    else
        addparameter(model,['start' TumorArray{array}],0,'ValueUnits','dimensionless','ConstantValue',false);
    end
    addevent(model,['seeding && time>=delay' TumorArray{array}],['start' TumorArray{array} '=1']);
    for nc=1:ncancer_clones
        addevent(model,['seeding && time>=delay' TumorArray{array}],['V_T' TumorArray{array} '.C' num2str(nc) '=' 'ncells_C' num2str(nc) TumorArray{array}]);
    end

end
   
%   addevent(model,['seeding && time>=delay' TumorArray{array}],['V_T' TumorArray{array} '.TGF=0*' 'V_T' TumorArray{array} '.TGF']);
%   addevent(model,['seeding && time>=delay' TumorArray{array}],['V_T' TumorArray{array} '.IL12=0*' 'V_T' TumorArray{array} '.IL12']);
%   addevent(model,['seeding && time>=delay' TumorArray{array}],['V_T' TumorArray{array} '.NO=0*' 'V_T' TumorArray{array} '.NO']);
%   addevent(model,['seeding && time>=delay' TumorArray{array}],['V_T' TumorArray{array} '.ArgI=0*' 'V_T' TumorArray{array} '.ArgI']);
%   addevent(model,['seeding && time>=delay' TumorArray{array}],['V_T' TumorArray{array} '.c=0*' 'V_T' TumorArray{array} '.c']);
%   addevent(model,['seeding && time>=delay' TumorArray{array}],['V_T' TumorArray{array} '.MDSC=0*' 'V_T' TumorArray{array} '.MDSC']);
%   addevent(model,['seeding && time>=delay' TumorArray{array}],['V_T' TumorArray{array} '.APC=0*' 'V_T' TumorArray{array} '.APC']);
%   addevent(model,['seeding && time>=delay' TumorArray{array}],['V_T' TumorArray{array} '.mAPC=0*' 'V_T' TumorArray{array} '.mAPC']);
% addevent(model,'seeding==0','start_Ln1=1');
% addrule(model,'k_C1_growth_Ln1= start*k_C1_growth','repeatedAssignment');
% addrule(model,'k_C2_growth_Ln1= start*k_C2_growth','repeatedAssignment');

% surgery - removal of primary tumor
% To switch off set do_surgery=0
addparameter(model,'do_surgery',0,'ValueUnits','dimensionless','ConstantValue',false);
addparameter(model,'surgery_time',params.surgery_time.Value,'ValueUnits','day');

% addevent(model,'do_surgery && time>=surgery_time','V_T.C1=0*cell');
% addevent(model,'do_surgery && time>=surgery_time','V_T.C2=0*cell');

for nc=1:ncancer_clones
    addevent(model,'do_surgery && time>=surgery_time',['V_T.C' num2str(nc) '=0*cell']);
end

for nt=1:n_T_specs
    addevent(model,'do_surgery && time>=surgery_time',['V_T.T' num2str(nt) '=0*cell']);
end

addevent(model,'do_surgery && time>=surgery_time','V_T.T0=0*cell');
addevent(model,'do_surgery && time>=surgery_time','V_T.Th=0*cell');
addevent(model,'do_surgery && time>=surgery_time','V_T.Mac_M1=0*cell');
addevent(model,'do_surgery && time>=surgery_time','V_T.Mac_M2=0*cell');
addevent(model,'do_surgery && time>=surgery_time','V_T.APC=0*cell');
addevent(model,'do_surgery && time>=surgery_time','V_T.mAPC=0*cell');
addevent(model,'do_surgery && time>=surgery_time','V_T.MDSC=0*cell');
addevent(model,'do_surgery && time>=surgery_time','V_T.T1_exh=0*cell');
addevent(model,'do_surgery && time>=surgery_time','V_T.Th_exh=0*cell');
addevent(model,'do_surgery && time>=surgery_time','V_T.C_x=0*cell');

% addevent(model,'do_surgery && time>=surgery_time','V_T.CCL2=0*V_T.CCL2');
% addevent(model,'do_surgery && time>=surgery_time','V_T.IL12=0*V_T.IL12');
% addevent(model,'do_surgery && time>=surgery_time','V_T.IL10=0*V_T.IL10');
% addevent(model,'do_surgery && time>=surgery_time','V_T.P1=0*V_T.P1');
% addevent(model,'do_surgery && time>=surgery_time','V_T.P0=0*V_T.P0');
% addevent(model,'do_surgery && time>=surgery_time','V_T.IFNg=0*V_T.IFNg');
% addevent(model,'do_surgery && time>=surgery_time','V_T.TGF=0*V_T.TGF');
% addevent(model,'do_surgery && time>=surgery_time','V_T.ArgI=0*V_T.ArgI');
% addevent(model,'do_surgery && time>=surgery_time','V_T.NO=0*V_T.NO');
% addevent(model,'do_surgery && time>=surgery_time','V_T.c_vas=0*V_T.c_vas');
% addevent(model,'do_surgery && time>=surgery_time','V_T.c=0*V_T.c');

end


