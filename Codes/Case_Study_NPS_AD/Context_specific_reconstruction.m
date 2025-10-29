%NPS_AD
setenv('GRB_LICENSE_FILE', '/depot/pbaloni/data/Lab_members/Boyu_Jiang/Software/gurobi_license/gurobi.lic');

%load("/home/jiang817/BaloniLab_Depot/data/Lab_members/Boyu_Jiang/Neuron_GEM/Refinement/Final_reconstruction/iNeuron_final_11092024.mat");
load("/home/jiang817/BaloniLab_Depot/data/Lab_members/Boyu_Jiang/Neuron_GEM/Refinement/Final_reconstruction/ineuron_metabolome_01062025.mat");
model=generateGrRules(Human_Neuron_GEM);


% Set media
exRxns = strmatch('EX_', model.rxns);
exRxns_names = model.rxns(exRxns);
model = changeRxnBounds(model,exRxns_names,0,'l');

%% CSF -- 100
model = changeRxnBounds(model,'EX_gln__L_e',-100,'l'); % not in the iNeuron
model = changeRxnBounds(model,'EX_gln__L_e',100,'u');
model = changeRxnBounds(model,'EX_acac_e',-100,'l');
model = changeRxnBounds(model,'EX_acac_e',-100,'u');
%model = changeRxnBounds(model,'r1391',-100,'l');
%model = changeRxnBounds(model,'r1391',0,'u');
%model = changeRxnBounds(model,'r1392',-100,'l');
%model = changeRxnBounds(model,'r1392',0,'u');
model = changeRxnBounds(model,'G3PD',-100,'l');
model = changeRxnBounds(model,'G3PD',0,'u');
model = changeRxnBounds(model,'EX_h2o2_e',-100,'l');
model = changeRxnBounds(model,'EX_h2o2_e',0,'u');
model = changeRxnBounds(model,'EX_hco3_e',-100,'l'); % not in the iNeuron
model = changeRxnBounds(model,'EX_hco3_e',100,'u');
model = changeRxnBounds(model,'EX_lac__L_e',-100,'l');
model = changeRxnBounds(model,'EX_lac__L_e',100,'u');
%model = changeRxnBounds(model,'EX_adrnl[e]',0,'l'); % not in the iNeuron
%model = changeRxnBounds(model,'EX_adrnl[e]',0,'u');
%model = changeRxnBounds(model,'EX_dopa[e]',0,'l');% not in the iNeuron
%model = changeRxnBounds(model,'EX_dopa[e]',0,'u');
%model = changeRxnBounds(model,'EX_nrpphr[e]',0,'l');
%model = changeRxnBounds(model,'EX_nrpphr[e]',0,'u');
model = changeRxnBounds(model,'EX_ile__L_e',-100,'l');
model = changeRxnBounds(model,'EX_ile__L_e',100,'u');
%model = changeRxnBounds(model,'EX_hxan[e]',-0.00045,'l');% not in the iNeuron
%model = changeRxnBounds(model,'EX_hxan[e]',0.00045,'u');
%model = changeRxnBounds(model,'EX_ins[e]',-0.00045,'l');% not in the iNeuron
%model = changeRxnBounds(model,'EX_ins[e]',0.00045,'u');
model = changeRxnBounds(model,'EX_ser__L_e',-100,'l');
model = changeRxnBounds(model,'EX_ser__L_e',100,'u');
model = changeRxnBounds(model,'EX_val__L_e',-100,'l');
model = changeRxnBounds(model,'EX_val__L_e',100,'u');
model = changeRxnBounds(model,'EX_asn__L_e',-100,'l');
model = changeRxnBounds(model,'EX_asn__L_e',100,'u');
model = changeRxnBounds(model,'EX_pro__L_e',-100,'l');
model = changeRxnBounds(model,'EX_pro__L_e',100,'u');
model = changeRxnBounds(model,'EX_pyr_e',-100,'l');
model = changeRxnBounds(model,'EX_pyr_e',100,'u');
model = changeRxnBounds(model,'EX_ala__L_e',-100,'l');
model = changeRxnBounds(model,'EX_ala__L_e',100,'u');
model = changeRxnBounds(model,'EX_gly_e',-100,'l');
model = changeRxnBounds(model,'EX_gly_e',100,'u');
model = changeRxnBounds(model,'EX_lys__L_e',-100,'l');
model = changeRxnBounds(model,'EX_lys__L_e',100,'u');
model = changeRxnBounds(model,'EX_glu__L_e',0,'l');
model = changeRxnBounds(model,'EX_glu__L_e',100,'u');
model = changeRxnBounds(model,'EX_glc__D_e',-100,'l');
model = changeRxnBounds(model,'EX_glc__D_e',-100,'u');
model = changeRxnBounds(model,'EX_phe__L_e',-100,'l');
model = changeRxnBounds(model,'EX_phe__L_e',100,'u');
model = changeRxnBounds(model,'EX_trp__L_e',-100,'l');
model = changeRxnBounds(model,'EX_trp__L_e',100,'u');
model = changeRxnBounds(model,'GLCt1',-100,'l');
model = changeRxnBounds(model,'GLCt1',-100,'u');
model = changeRxnBounds(model,'EX_thr__L_e',-100,'l');
model = changeRxnBounds(model,'EX_thr__L_e',100,'u');
model = changeRxnBounds(model,'EX_met__L_e',-100,'l');
model = changeRxnBounds(model,'EX_met__L_e',100,'u');
model = changeRxnBounds(model,'EX_tyr__L_e',-100,'l');
model = changeRxnBounds(model,'EX_tyr__L_e',100,'u');
model = changeRxnBounds(model,'EX_arg__L_e',-100,'l');
model = changeRxnBounds(model,'EX_arg__L_e',100,'u');
model = changeRxnBounds(model,'EX_his__L_e',-100,'l');
model = changeRxnBounds(model,'EX_his__L_e',100,'u');
model = changeRxnBounds(model,'EX_orn_e',-100,'l');
model = changeRxnBounds(model,'EX_orn_e',100,'u');
%model = changeRxnBounds(model,'EX_Lcystin[e]',0,'l'); % not in the iNeuron
%model = changeRxnBounds(model,'EX_Lcystin[e]',0.0045,'u');
model = changeRxnBounds(model,'EX_leu__L_e',-100,'l');
model = changeRxnBounds(model,'EX_leu__L_e',100,'u');
model = changeRxnBounds(model,'EX_o2_e',-100,'l');
model = changeRxnBounds(model,'EX_o2_e',-100,'u');
model = changeRxnBounds(model,'EX_co2_e',100,'l');
model = changeRxnBounds(model,'EX_co2_e',100,'u');
%model = changeRxnBounds(model,'EX_3aib[e]',-0.00023,'l');% not in the iNeuron
%model = changeRxnBounds(model,'EX_4abut_e',-100,'l');
%model = changeRxnBounds(model,'EX_4abut_e',0,'u');
model = changeRxnBounds(model,'EX_ac_e',-100,'l');
model = changeRxnBounds(model,'EX_ac_e',100,'u');
%model = changeRxnBounds(model,'EX_acald[e]',-0.0014,'l');% not in the iNeuron
%model = changeRxnBounds(model,'EX_acald[e]',1000,'u');
model = changeRxnBounds(model,'EX_ala_B_e',-100,'l');
model = changeRxnBounds(model,'EX_ala_B_e',100,'u');
model = changeRxnBounds(model,'EX_bhb_e',-100,'l');
model = changeRxnBounds(model,'EX_bhb_e',100,'u');
model = addExchangeRxn(model, {'btn_e', 'btn_c'});
model = changeRxnBounds(model,'EX_btn_e',-100,'l');
model = changeRxnBounds(model,'EX_btn_e',100,'u');
model = changeRxnBounds(model,'EX_crn_e',-100,'l');
model = changeRxnBounds(model,'EX_crn_e',100,'u');
model = changeRxnBounds(model,'EX_cys__L_e',-100,'l');
model = changeRxnBounds(model,'EX_cys__L_e',100,'u');
%model = changeRxnBounds(model,'EX_duri[e]',1000,'u');% not in the iNeuron
%model = changeRxnBounds(model,'EX_duri[e]',-1000,'l');
model = changeRxnBounds(model,'EX_etoh_e',-100,'l');
model = changeRxnBounds(model,'EX_etoh_e',100,'u');
model = changeRxnBounds(model,'EX_fe2[e]',-100,'l');% not in the iNeuron
model = changeRxnBounds(model,'EX_fe2[e]',100,'u');
model = changeRxnBounds(model,'EX_nh4_e',-100,'l');
model = changeRxnBounds(model,'EX_nh4_e',100,'u');
model = changeRxnBounds(model,'EX_pe_hs_e',-100,'l');
model = changeRxnBounds(model,'EX_pe_hs_e',100,'u');
%model = changeRxnBounds(model,'EX_sarcs[e]',-0.0053,'l'); % not in the iNeuron
%model = changeRxnBounds(model,'EX_sarcs[e]',0.0086,'u');
model = changeRxnBounds(model,'EX_h2o_e',-100,'l');
model = changeRxnBounds(model,'EX_h2o_e',100,'u');
model = changeRxnBounds(model,'EX_pi_e',-100,'l');
model = changeRxnBounds(model,'EX_h_e',-100,'l');
%% CSF -- ORI
model = changeRxnBounds(model,'EX_gln__L_e',-0.013,'l'); % not in the iNeuron
model = changeRxnBounds(model,'EX_gln__L_e',0.025,'u');
model = changeRxnBounds(model,'EX_acac_e',-0.012,'l');
model = changeRxnBounds(model,'EX_acac_e',-0.0015,'u');
model = changeRxnBounds(model,'r1391',-100,'l');
model = changeRxnBounds(model,'r1391',0,'u');
model = changeRxnBounds(model,'r1392',-100,'l');
model = changeRxnBounds(model,'r1392',0,'u');
%model = changeRxnBounds(model,'G3PD',-100,'l');
%model = changeRxnBounds(model,'G3PD',0,'u');
model = changeRxnBounds(model,'EX_h2o2_e',-100,'l');
model = changeRxnBounds(model,'EX_h2o2_e',0,'u');
model = changeRxnBounds(model,'EX_hco3_e',-100,'l'); % not in the iNeuron
model = changeRxnBounds(model,'EX_hco3_e',0,'u');
model = changeRxnBounds(model,'EX_lac__L_e',-0.0058,'l');
model = changeRxnBounds(model,'EX_lac__L_e',0.079,'u');
%model = changeRxnBounds(model,'EX_adrnl[e]',0,'l'); % not in the iNeuron
%model = changeRxnBounds(model,'EX_adrnl[e]',0,'u');
%model = changeRxnBounds(model,'EX_dopa[e]',0,'l');% not in the iNeuron
%model = changeRxnBounds(model,'EX_dopa[e]',0,'u');
%model = changeRxnBounds(model,'EX_nrpphr[e]',0,'l');
%model = changeRxnBounds(model,'EX_nrpphr[e]',0,'u');
model = changeRxnBounds(model,'EX_ile__L_e',-0.0041,'l');
model = changeRxnBounds(model,'EX_ile__L_e',0.0004,'u');
%model = changeRxnBounds(model,'EX_hxan[e]',-0.00045,'l');% not in the iNeuron
%model = changeRxnBounds(model,'EX_hxan[e]',0.00045,'u');
%model = changeRxnBounds(model,'EX_ins[e]',-0.00045,'l');% not in the iNeuron
%model = changeRxnBounds(model,'EX_ins[e]',0.00045,'u');
model = changeRxnBounds(model,'EX_ser__L_e',-0.011,'l');
model = changeRxnBounds(model,'EX_ser__L_e',0.0016,'u');
model = changeRxnBounds(model,'EX_val__L_e',-0.011,'l');
model = changeRxnBounds(model,'EX_val__L_e',0.005,'u');
model = changeRxnBounds(model,'EX_asn__L_e',-0.0009,'l');
model = changeRxnBounds(model,'EX_asn__L_e',0.0037,'u');
model = changeRxnBounds(model,'EX_pro__L_e',-0.0079,'l');
model = changeRxnBounds(model,'EX_pro__L_e',0.0066,'u');
model = changeRxnBounds(model,'EX_pyr_e',-0.0058,'l');
model = changeRxnBounds(model,'EX_pyr_e',0.007,'u');
model = changeRxnBounds(model,'EX_ala__L_e',-100,'l');
model = changeRxnBounds(model,'EX_ala__L_e',0.0079,'u');
model = changeRxnBounds(model,'EX_gly_e',-0.0053,'l');
model = changeRxnBounds(model,'EX_gly_e',0.0086,'u');
model = changeRxnBounds(model,'EX_lys__L_e',-0.0005,'l');
model = changeRxnBounds(model,'EX_lys__L_e',0.011,'u');
model = changeRxnBounds(model,'EX_glu__L_e',-0.0044,'l');
model = changeRxnBounds(model,'EX_glu__L_e',0.0047,'u');
model = changeRxnBounds(model,'EX_glc__D_e',-0.29,'l');
model = changeRxnBounds(model,'EX_glc__D_e',-0.196,'u');
model = changeRxnBounds(model,'EX_phe__L_e',0,'l');
model = changeRxnBounds(model,'EX_phe__L_e',100,'u');
model = changeRxnBounds(model,'EX_trp__L_e',0,'l');
model = changeRxnBounds(model,'EX_trp__L_e',100,'u');
model = changeRxnBounds(model,'GLCt1',-0.19,'l');
model = changeRxnBounds(model,'GLCt1',-0.16,'u');
model = changeRxnBounds(model,'EX_thr__L_e',0,'l');
model = changeRxnBounds(model,'EX_thr__L_e',0.0008,'u');
model = changeRxnBounds(model,'EX_met__L_e',0,'l');
model = changeRxnBounds(model,'EX_met__L_e',0.0017,'u');
model = changeRxnBounds(model,'EX_tyr__L_e',0,'l');
model = changeRxnBounds(model,'EX_tyr__L_e',0.0017,'u');
model = changeRxnBounds(model,'EX_arg__L_e',-0.004,'l');
model = changeRxnBounds(model,'EX_arg__L_e',0,'u');
model = changeRxnBounds(model,'EX_his__L_e',0,'l');
model = changeRxnBounds(model,'EX_his__L_e',0.0025,'u');
model = changeRxnBounds(model,'EX_orn_e',-0.0048,'l');
model = changeRxnBounds(model,'EX_orn_e',0.0041,'u');
%model = changeRxnBounds(model,'EX_Lcystin[e]',0,'l'); % not in the iNeuron
%model = changeRxnBounds(model,'EX_Lcystin[e]',0.0045,'u');
model = changeRxnBounds(model,'EX_leu__L_e',-0.0062,'l');
model = changeRxnBounds(model,'EX_leu__L_e',0.0011,'u');
model = changeRxnBounds(model,'EX_o2_e',-2.256,'l');
model = changeRxnBounds(model,'EX_o2_e',-1.351,'u');
model = changeRxnBounds(model,'EX_co2_e',0.515,'l');
model = changeRxnBounds(model,'EX_co2_e',0.530,'u');
%model = changeRxnBounds(model,'EX_3aib[e]',-0.00023,'l');% not in the iNeuron
model = changeRxnBounds(model,'EX_co2_e',1000,'u');
%model = changeRxnBounds(model,'EX_4abut[e]',-0.0018,'l');
%model = changeRxnBounds(model,'EX_4abut[e]',0,'u');
model = changeRxnBounds(model,'EX_ac_e',-0.0013,'l');
model = changeRxnBounds(model,'EX_ac_e',1000,'u');
%model = changeRxnBounds(model,'EX_acald[e]',-0.0014,'l');% not in the iNeuron
%model = changeRxnBounds(model,'EX_acald[e]',1000,'u');
model = changeRxnBounds(model,'EX_ala_B_e',-1000,'l');
model = changeRxnBounds(model,'EX_ala_B_e',1000,'u');
model = changeRxnBounds(model,'EX_bhb_e',-0.016,'l');
model = changeRxnBounds(model,'EX_bhb_e',-0.001,'u');
%model = changeRxnBounds(model,'EX_btn_e',-1000,'l');
%model = changeRxnBounds(model,'EX_btn_e',1000,'u');
model = changeRxnBounds(model,'EX_crn_e',-1000,'l');
model = changeRxnBounds(model,'EX_crn_e',1000,'u');
model = changeRxnBounds(model,'EX_cys__L_e',-0.0086,'l');
model = changeRxnBounds(model,'EX_cys__L_e',0.0033,'u');
%model = changeRxnBounds(model,'EX_duri[e]',1000,'u');% not in the iNeuron
%model = changeRxnBounds(model,'EX_duri[e]',-1000,'l');
model = changeRxnBounds(model,'EX_etoh_e',0,'l');
model = changeRxnBounds(model,'EX_etoh_e',0,'u');
%model = changeRxnBounds(model,'EX_fe2[e]',-1000,'l');% not in the iNeuron
%model = changeRxnBounds(model,'EX_fe2[e]',1000,'u');
model = changeRxnBounds(model,'EX_nh4_e',-1000,'l');
model = changeRxnBounds(model,'EX_nh4_e',0,'u');
model = changeRxnBounds(model,'EX_pe_hs_e',-0.00005,'l');
model = changeRxnBounds(model,'EX_pe_hs_e',1000,'u');
%model = changeRxnBounds(model,'EX_sarcs[e]',-0.0053,'l'); % not in the iNeuron
%model = changeRxnBounds(model,'EX_sarcs[e]',0.0086,'u');
model = changeRxnBounds(model,'EX_h2o_e',-1000,'l');
model = changeRxnBounds(model,'EX_h2o_e',1000,'u');

%% BrainPhys -- -100
%metList = {"na1_e", "ca2_e", "mg2_e", "zn2_e", "zn2_e", "fe3_e", "inost_e", "ncam_e", "pydxn_e", "thm_e"}
%model = addExchangeRxn(model, {"inost_e"}, -100, 100);

model = changeRxnBounds(model,'EX_na1_e',-100,'l');
model = changeRxnBounds(model,'EX_cl_e',-100,'l');
model = changeRxnBounds(model,'EX_k_e',-100,'l');
model = changeRxnBounds(model,'EX_ca2_e',-100,'l');
model = changeRxnBounds(model,'EX_mg2_e',-100,'l');
model = changeRxnBounds(model,'EX_so4_e',-100,'l');
model = changeRxnBounds(model,'EX_zn2_e',-100,'l');
model = changeRxnBounds(model,'EX_fe3_e',-100,'l');
model = changeRxnBounds(model,'EX_hco3_e',-100,'l');
model = changeRxnBounds(model,'EX_pi_e',-100,'l');
model = changeRxnBounds(model,'EX_h_e',-100,'l');
model = changeRxnBounds(model,'EX_ile__L_e',-100,'l');
model = changeRxnBounds(model,'EX_thr__L_e',-100,'l');
model = changeRxnBounds(model,'EX_leu__L_e',-100,'l');
model = changeRxnBounds(model,'EX_val__L_e',-100,'l');
model = changeRxnBounds(model,'EX_lys__L_e',-100,'l');
model = changeRxnBounds(model,'EX_gln__L_e',-100,'l');
model = changeRxnBounds(model,'EX_arg__L_e',-100,'l');
model = changeRxnBounds(model,'EX_tyr__L_e',-100,'l');
model = changeRxnBounds(model,'EX_phe__L_e',-100,'l');
model = changeRxnBounds(model,'EX_his__L_e',-100,'l');
model = changeRxnBounds(model,'EX_met__L_e',-100,'l');
model = changeRxnBounds(model,'EX_pro__L_e',-100,'l');
model = changeRxnBounds(model,'EX_cys__L_e',-100,'l');
model = changeRxnBounds(model,'EX_trp__L_e',-100,'l');
model = changeRxnBounds(model,'EX_asn__L_e',-100,'l');
model = changeRxnBounds(model,'EX_ala__L_e',-100,'l');
model = changeRxnBounds(model,'EX_gly_e',-100,'l');
model = changeRxnBounds(model,'EX_ser__L_e',-100,'l');
model = changeRxnBounds(model,'EX_chol_e',-100,'l');
model = changeRxnBounds(model,'EX_inost_e',-100,'l');
model = changeRxnBounds(model,'EX_ncam_e',-100,'l');
model = changeRxnBounds(model,'EX_pydxn_e',-100,'l');
model = changeRxnBounds(model,'EX_thm_e',-100,'l');
model = changeRxnBounds(model,'EX_pnto__R_e',-100,'l');
model = changeRxnBounds(model,'EX_fol_e',-100,'l');
model = changeRxnBounds(model,'EX_ribflv_e',-100,'l');
model = changeRxnBounds(model,'EX_glc__D_e',-100,'l');
model = changeRxnBounds(model,'EX_pyr_e',-100,'l');
model = changeRxnBounds(model,'EX_chsterol_e',-100,'l');
model = changeRxnBounds(model,'EX_o2_e',-100,'l');
model = changeRxnBounds(model,'EX_co2_e',-100,'l');
model = changeRxnBounds(model,'EX_h2o_e',-100,'l');
%% CSF --hmdb
[num, txt, raw] = xlsread('/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Refinement/Metabolomics/final_reconstruction_medium_constrains.xlsx')
for i = 2:1:282
    exchange = raw{i,1};
    disp(exchange)
    lower_bound = raw{i,6};
    upper_bound = raw{i,7};
    model = changeRxnBounds(model,exchange,lower_bound,'l');
    model = changeRxnBounds(model,exchange,upper_bound,'u');
end
model = changeRxnBounds(model,'BIOMASS_maintenance',1,'u');
model = changeRxnBounds(model,'BIOMASS_maintenance',1,'l');
%% set objective

GABA_obj = {'BIOMASS_maintenance', 'Synapse_GABAVESSEC_1'};
%GLU_obj = {'BIOMASS_maintenance', 'Synapse_GLUVESSEC'};

model = changeObjective (model, GABA_obj);
%model = changeObjective (model, GLU_obj);
checkObjective(model)

core = {'BIOMASS_maintenance', 'Synapse_GABAVESSEC_1', 'DM_4abut_syna'};

%folderPath = "/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/expressionSet_CPM/EN/"
folderPath = "/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/expressionSet_CPM/IN/"

BIOMASS_index = find(strcmp(model.rxns, "BIOMASS_maintenance"))
SYNAPSE_GABA_index = find(strcmp(model.rxns, "Synapse_GABAVESSEC_1"))
DM_4abut_syna_index = find(strcmp(model.rxns, "DM_4abut_syna"))

files = dir(folderPath);
files(1:2) = [];
fileNames = {files(~[files.isdir]).name};
for i = 1:1:1479
%for i = 1328:1:1479
    disp(i)
    index = i;
    index = num2str(index);
    disp(fileNames{i});
    name_of_file = fileNames{i};
    name_of_file = split(name_of_file, '.');
    name_of_file = name_of_file{1};
    name_of_file = split(name_of_file, '_CPM_ALRA');
    name_of_file = name_of_file{1};
    filepath = strcat(folderPath,fileNames{i});
    data = readtable(filepath, "FileType","text",'Delimiter', ',');
    exprData.gene=data.BIGG
    %exprData.value=data.EX_CPM
    exprData.value=data.IN_CPM
    rxn_expression=mapExpressionToReactions(model,exprData);

    % EN_75percentile = 44.97852
    IN_75percentile = 72.51399 
    th_ub = IN_75percentile;
    th_lb = 0;

  

    for k=1:length(rxn_expression)
        if isnan(rxn_expression(k))==1;
            rxn_expression(k)=-1;
        end
    end

    rxn_expression(BIOMASS_index,1) = prctile(rxn_expression,90);
    rxn_expression(SYNAPSE_GABA_index,1) = prctile(rxn_expression,90);
    rxn_expression(DM_4abut_syna_index,1) = prctile(rxn_expression,90);

    %for i = 2:1:146
    %    exchange = raw{i,1};
    %    exchange_index = find(strcmp(model.rxns, exchange));
    %    rxn_expression(exchange_index,1) = prctile(rxn_expression,90);
    %end
    tol = 1e-8;
    

   

    %iMAT_reconstruction = iMAT(model, rxn_expression,0,th_ub,tol, GLU_obj);
    iMAT_reconstruction = iMAT(model, rxn_expression, th_lb, th_ub, tol, core);
    

    SAVE_FILENAME = strrep(name_of_file,'_ExpressionSet','_GABA')

    %SAVE_FILENAME = strrep(name_of_file,'_ExpressionSet','_GLU')

    %options.core = GABA_obj;
    %options.core = GLU_obj;
    %options.solver ='iMAT';
    %options.expressionRxns=rxn_expression;
    %options.threshold_lb=th_lb;
    %options.threshold_ub=th_ub;
    %iMAT_reconstruction=createTissueSpecificModel(model, options);


    modesavePath = strcat('/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/iMAT_reconstruction_01162025/IN/', SAVE_FILENAME, ".mat")
    %modesavePath = strcat('/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/iMAT_reconstruction/EX/', SAVE_FILENAME, ".mat")

    %modesavePath = strcat('/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/PMID_38417436/iMAT_reconstructions_With_ALRA/GABA_CSF_100_th_ub_75/', SAVE_FILENAME, ".mat")
    %modesavePath = strcat('/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/PMID_38417436/iMAT_reconstructions_With_ALRA/GLUTAMATE_CSF_100_th_ub_75/', SAVE_FILENAME, ".mat")

    
    writeCbModel(iMAT_reconstruction, 'mat', modesavePath);
    
end

%% Check optimization of reconstruction
folderPath = "/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/iMAT_reconstruction_01132025/IN/"
files = dir(folderPath);
files(1:2) = [];
fileNames = {files(~[files.isdir]).name};
GABA_obj = {'BIOMASS_maintenance', 'Synapse_GABAVESSEC_1'};
%GLU_obj = {'BIOMASS_maintenance', 'Synapse_GLUVESSEC'};



for i = 1:1:1327
    disp(i)
    index = i;
    index = num2str(index);
    disp(fileNames{i});
    
    filepath = strcat(folderPath,fileNames{i});
    load(filepath);
    model=generateGrRules(Human_Neuron_GEM);
    model = changeObjective (model, GABA_obj);

    model = changeRxnBounds(model,'BIOMASS_maintenance',1,'u');
    model = changeRxnBounds(model,'BIOMASS_maintenance',1,'l');

    
    BIOMASS_index = find(strcmp(model.rxns, "BIOMASS_maintenance"))
    SYNAPSE_GABA_index = find(strcmp(model.rxns, "Synapse_GABAVESSEC_1"))

    sol = optimizeCbModel(model)

    
