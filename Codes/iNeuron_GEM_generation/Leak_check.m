%% LEAK TEAT
%load("/home/jiang817/BaloniLab_Depot/data/Lab_members/Boyu_Jiang/Neuron_GEM/Refinement/Final_reconstruction/iNeuron_final_09282024.mat");
load("/home/jiang817/BaloniLab_Depot/data/Lab_members/Boyu_Jiang/Neuron_GEM/Refinement/Final_reconstruction/iNeuron_final_10042024.mat");

model = generateGrRules (Human_Neuron_GEM)
% Set lower bounds of the biomass reactions to 0
model.lb(find(ismember(model.rxns, 'BIOMASS_maintenance'))) = 0
cnt = 1
tol = 1e-7
%% add DM_atp_c
if isempty(strmatch('DM_atp_c',model.rxns))
    [model, rxnIDexists] = addReaction(model,'DM_atp_c', 'reactionFormula', 'h2o_c + atp_c  -> adp_c + h_c + pi_c');
end
 
%% Close model
modelClosed = model;
modelexchanges1 = strmatch('EX_', modelClosed.rxns);
modelexchanges2 = strmatch('DM_', modelClosed.rxns);
modelexchanges3 = strmatch('SK_', modelClosed.rxns);
selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';
modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;selExc]);
modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
modelClosed.ub(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=1000;
modelClosedOri = modelClosed;
%% Start with tests.
% Perform leak test, i.e., whether the closed model can produce any exchanged 
% metabolite, as defined in the model, from nothing. 
modelClosed = modelClosedOri;
[LeakRxns,modelTested,LeakRxnsFluxVector] = fastLeakTest(modelClosed,modelClosed.rxns(selExc),'false');
TableChecks{cnt,1} = 'fastLeakTest 1';
if length(LeakRxns)>0
    warning('model leaks metabolites!')
    TableChecks{cnt,2} = 'Model leaks metabolites!';
else
    TableChecks{cnt,2} = 'Leak free!';
end
cnt = cnt + 1;
%% 
% Test if something leaks when demand reactions for each metabolite in the 
% model are added. Note that this step is time consuming.

modelClosed = modelClosedOri;
[LeakRxnsDM,modelTestedDM,LeakRxnsFluxVectorDM] = fastLeakTest(modelClosed,modelClosed.rxns(selExc),'true');

TableChecks{cnt,1} = 'fastLeakTest 2 - add demand reactions for each metabolite in the model';
if length(LeakRxnsDM)>0
    TableChecks{cnt,2} = 'Model leaks metabolites when demand reactions are added!';
else
    TableChecks{cnt,2} = 'Leak free when demand reactions are added!';
end
cnt = cnt + 1;
%% 
% Test if the model produces energy from water!

modelClosed = modelClosedOri;
% model does not have 'DM_atp_c_', so manually added it to the model
modelClosed = addDemandReaction(modelClosed, {'atp_c'})
modelClosedATP = changeObjective(modelClosed,'DM_atp_c');

modelClosedATP = changeRxnBounds(modelClosedATP,'DM_atp_c',0,'l');
modelClosedATP = changeRxnBounds(modelClosedATP,'EX_h2o_e',-1,'l');
FBA3=optimizeCbModel(modelClosedATP);
TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, except h2o';
if abs(FBA3.f) > 1e-6
    TableChecks{cnt,2} = 'model produces energy from water!';
else
    TableChecks{cnt,2} = 'model DOES NOT produce energy from water!';
end
cnt = cnt + 1;

%%
% Test if the model produces energy from water and oxygen!

modelClosed = modelClosedOri;
modelClosed = addDemandReaction(modelClosed, {'atp_c'})
modelClosedATP = changeObjective(modelClosed,'DM_atp_c');
modelClosedATP = changeRxnBounds(modelClosedATP,'DM_atp_c',0,'l');
modelClosedATP = changeRxnBounds(modelClosedATP,'EX_h2o_e',-1,'l');
modelClosedATP = changeRxnBounds(modelClosedATP,'EX_o2_e',-1,'l');

FBA6=optimizeCbModel(modelClosedATP);
TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, except h2o and o2';
if abs(FBA6.f) > 1e-6
    TableChecks{cnt,2} = 'model produces energy from water and oxygen!';
else
    TableChecks{cnt,2} = 'model DOES NOT produce energy from water and oxygen!';
end
cnt = cnt + 1;
%%
% Test if the model produces matter when atp demand is reversed!

modelClosed = modelClosedOri;
modelClosed = addDemandReaction(modelClosed, {'atp_c'})
modelClosed = changeObjective(modelClosed,'DM_atp_c');
modelClosed.lb(find(ismember(modelClosed.rxns,'DM_atp_c'))) = -1000;
FBA = optimizeCbModel(modelClosed);
TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, allow DM_atp_c to be reversible';
if abs(FBA.f) > 1e-6
    TableChecks{cnt,2} = 'model produces matter when atp demand is reversed!';
else
    TableChecks{cnt,2} = 'model DOES NOT produce matter when atp demand is reversed!';
end
cnt = cnt + 1;
%% 
% Test if the model has flux through h[m] demand !

modelClosed = modelClosedOri;
modelClosed = addDemandReaction(modelClosed,'h_m');
modelClosed = changeObjective(modelClosed,'DM_h_m');
modelClosed.ub(find(ismember(modelClosed.rxns,'DM_h_m'))) = 1000;
FBA = optimizeCbModel(modelClosed,'max');
TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, test flux through DM_h_m (max)';
if abs(FBA.f) > 1e-6
    TableChecks{cnt,2} = 'model has flux through h_m demand (max)!';
else
    TableChecks{cnt,2} = 'model has NO flux through h_m demand (max)!';
end
cnt = cnt + 1;
%% 
% Test if the  model has flux through h_c demand !

modelClosed = modelClosedOri;
modelClosed = addDemandReaction(modelClosed,'h_c');
modelClosed = changeObjective(modelClosed,'DM_h_c');
modelClosed.ub(find(ismember(modelClosed.rxns,'DM_h_c'))) = 1000;
FBA = optimizeCbModel(modelClosed,'max');
TableChecks{cnt,1} = 'Exchanges, sinks, and demands have  lb = 0, test flux through DM_h[c] (max)';
if abs(FBA.f) > 1e-6
    TableChecks{cnt,2} = 'model has flux through h_c demand (max)!';
else
    TableChecks{cnt,2} = 'model has NO flux through h_c demand (max)!';
end
cnt = cnt +1 
%% 
% Test if the  model produces too much atp demand from glucose under aerobic 
% condition. Also consider using the tutorial testModelATPYield to test if the 
% correct ATP yield from different carbon sources can be realized by the model.

modelClosed = modelClosedOri;
modelClosed = addDemandReaction(modelClosed, {'atp_c'});
modelClosed = changeObjective(modelClosed,'DM_atp_c');
modelClosed.lb(find(ismember(modelClosed.rxns,'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns,'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns,'EX_h2o_e'))) = 1000;
modelClosed.ub(find(ismember(modelClosed.rxns,'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns,strcat('EX_glc__D_e')))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns,strcat('EX_glc__D_e')))) = -1;
modelClosed.ub(selExc)=1000;

FBAOri_colon = optimizeCbModel(modelClosed,'max');
FBAOri = optimizeCbModel(modelClosed,'max');
TableChecks{cnt,1} = 'ATP yield ';
if abs(FBAOri.f) > 31 % this is the theoretical value
     TableChecks{cnt,2} = 'model produces too much atp demand from glc!';
else
   TableChecks{cnt,2} ='model DOES NOT produce too much atp demand from glc!';
end
cnt = cnt + 1;
%% 
% Test metabolic objective functions with open sinks. Note this step 
% is time consuming and may only work reliably on Recon 3D derived models due 
% to different usage of abbreviations.
model_test = model
TableChecks{cnt,1} = 'Test metabolic objective functions with open sinks';
if 1 % perform test function
    [TestSolution,TestSolutionNameOpenSinks, TestedRxnsSinks, PercSinks] = test4HumanFctExt(model,'all');
    TableChecks{cnt,2} = strcat('Done. See variable TestSolutionNameOpenSinks for results. The model passes ', num2str(length(find(abs(TestSolution)>tol))),' out of ', num2str(length(TestSolution)), 'tests');
else
    TableChecks{cnt,2} = 'Not performed.';
end
cnt =  cnt + 1;
%% 
% Test metabolic objective functions with closed sinks (lb). Note this step 
% is time consuming and may only work reliably on Recon 3D derived models due 
% to different usage of abbreviations.

TableChecks{cnt,1} = 'Test metabolic objective functions with closed sinks (lb)';
if 1 % perform test functions
    [TestSolution,TestSolutionNameClosedSinks, TestedRxnsClosedSinks, PercClosedSinks] = test4HumanFctExt(model,'all',0);
    TableChecks{cnt,2} = strcat('Done. See variable TestSolutionNameClosedSinks for results. The model passes ', num2str(length(find(abs(TestSolution)>tol))),' out of ', num2str(length(TestSolution)), 'tests');
else
    TableChecks{cnt,2} = 'Not performed.';
end
cnt =  cnt + 1;
%% 
% Compute ATP yield. This test is identical to the material covered in the 
% tutorial testModelATPYield.

TableChecks{cnt,1} = 'Compute ATP yield';
if 1 % test ATP yield
    [Table_csources, TestedRxns, Perc] = testATPYieldFromCsources(model);
    TableChecks{cnt,2} = 'Done. See variable Table_csources for results.';
else
    TableChecks{cnt,2} = 'Not performed.';
end
cnt = cnt + 1;
%% 
% Check for duplicated reactions in the model.

TableChecks{cnt,1} = 'Check duplicated reactions';
method='FR';
removeFlag=0;
[modelOut,removedRxnInd, keptRxnInd] = checkDuplicateRxn(model,method,removeFlag,0);
if isempty(removedRxnInd)
    TableChecks{cnt,2} = 'No duplicated reactions in model.';
else
    TableChecks{cnt,2} = 'Duplicated reactions in model.';
end
cnt = cnt + 1;
%% 
% Check empty columns in 'model.rxnGeneMat'.

TableChecks{cnt,1} = 'Check empty columns in rxnGeneMat';
E = find(sum(model.rxnGeneMat)==0);
if isempty(E)
    TableChecks{cnt,2} = 'No empty columns in rxnGeneMat.';
else
    TableChecks{cnt,2} = 'Empty columns in rxnGeneMat.';
end
cnt = cnt + 1;
%% 
% Check that demand reactions have a lb >= 0.

TableChecks{cnt,1} = 'Check that demand reactions have a lb >= 0';
DMlb = find(model.lb(strmatch('DM_',model.rxns))<0);
if isempty(DMlb)
    TableChecks{cnt,2} = 'No demand reaction can have flux in backward direction.';
else
    TableChecks{cnt,2} = 'Demand reaction can have flux in backward direction.';
end
cnt = cnt + 1;
%% 
% Check for flux consistency.

TableChecks{cnt,1} = 'Check for flux consistency';
param.epsilon=1e-4;
param.modeFlag=0;
%param.method='null_fastcc';
param.method='fastcc';
printLevel = 1;
[fluxConsistentMetBool,fluxConsistentRxnBool,fluxInConsistentMetBool,fluxInConsistentRxnBool,model] = findFluxConsistentSubset(modelOut,param,printLevel);
if isempty(find(fluxInConsistentRxnBool))
    TableChecks{cnt,2} = 'Model is flux consistent.';
else
    TableChecks{cnt,2} = 'Model is NOT flux consistent';
end
cnt = cnt + 1;
%% 
% Save all results.

resultsFileName = 'LeaktestResult_ineuron_10032024';
save(strcat(resultsFileName,'.mat'));


