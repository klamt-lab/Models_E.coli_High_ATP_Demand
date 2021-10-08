function [Fluxes, Conc, N, N_red, L, S, R, ParametersID,ParameterValue,MO,REQ,ntwkreactions,intracellular_index,Cnc] = ReadSBML()

% import SBML model. Simbiology toolbox required. 

% There are external and internal species 
MO = sbmlimport('ModelVersoni2_MonteCarlo_WT_CM.xml');
%%
% Analysis performed only on intracellular species. 
% Extracellular species and realtive entries in the stoichiometric matrix
% will be erased
% the compartment is indicated as BoundaryCondition
% External Metabolites: Boundary condition = 1
% Internal Boundary condition = 0

% symbiology function to identify conserved moieties
[SI, SD, L0, NR, ND] = sbioconsmoiety(MO,'link');


% Erase unwanted strings in metabolite names
% Find species with fixed values
vctDependent = [0];
for i = 1:1:numel(SD)
    current = SD(i);
    if contains(current,'Cytoplasm.')
       SD(i) = erase(current,'Cytoplasm.');
    else  
       SD(i) = erase(current,'Out.');
    end 
    
    for j=1:1:length(MO.Species)
        tofind = SD(i);
        ttlespec = MO.Species(j).Name;
        if strcmp(tofind,ttlespec)
           vctDependent(i) = j;
        end 
    end       
end 


index_m_in = 1;  
index_m_ex = 1; 
% Define intracellular and extracellular species. Get their index position
for i=1:1:length(MO.Species)
    cn = MO.Species(i).Parent.Name; 
    if contains(cn,'Out') 
        extracellular_index(index_m_ex) = i;
        index_m_ex = index_m_ex +1;
    else 
        intracellular_index(index_m_in) = i;
        index_m_in = index_m_in+1;
    end 
end 


% merge extracellular species indexes and dependent species 
merge_excell_dependent = [extracellular_index, vctDependent];
% erase reduntant values in the vector 
med = unique(merge_excell_dependent);

% define vector of species for the reduced stoichiometry 
indepSpecies = setdiff(intracellular_index,med);

%%
% Define number of reactions
nr = length(MO.Reactions);
ntwkreactions = linspace(1,nr,nr);

%%       
% retrieve stoichiometric matrix and reduce it 
% Stoichiometric matrix m rows x n columns 
% m = metabolites 
% n = reactions 

[N,objSpecies,objReactions] = getstoichmatrix(MO);
Nf = full(N);
% Delete external metabolites 
% retrieve and delete from the stoichtiometrix matrix all not necessary
% reactions (i.e. dilution reactions, transport-uptake-txcretion reactions)

% N_red = N(indepSpecies,ntwkreactions);
N = Nf(intracellular_index,ntwkreactions);
N_red = Nf(indepSpecies,ntwkreactions);


% Calculate link matrix.  A*x = B  ||  L*N_red = N  
L = N/N_red;

%% Import from .xlsx flies (generated from Copasi) steady state concentrations and fluxes
Cnc = readtable('SteadyStateConcentrations.xlsx');
Flx = readtable('SteadyStateFluxes.xlsx');

%% Define reaction and concentrations names 
S = {MO.Species.Name};           % Names of the metabolites
S = S(intracellular_index);
R = {MO.Reactions.Name};         % Names of the reactions
R = R(ntwkreactions);

%% Assign steady state concentration to MO structure 

% saved under MO.Species(i).Value
for i = 1:1:numel(MO.Species)
    in = MO.Species(i).Name;
    for j = 1:1:numel(Cnc.Name)
        ic = Cnc.Name(j);
        if contains(in,ic)
            MO.Species(i).Value = Cnc.Concentration(j);
        end 
    end 
end 


for k= 1:1:numel(intracellular_index)
    vl = intracellular_index(k);
    ssc(k) = MO.Species(vl).Value;
end 

% Select only intracellular concentrations 
Conc = ssc;

%% Assign steady state fluxes to MO structure. 
% The flux will be stored in MO.Reactions(i).UserData

for i = 1:1:numel(MO.Reactions)
    in = MO.Reactions(i).Name;
    for j = 1:1:numel(Flx.Name)
        ic = Flx.Name(j);
        if contains(in,ic)
           MO.Reactions(i).UserData = Flx.Flux(j);
        end 
    end 
end 

for k= 1:1:numel(ntwkreactions)
    vl = ntwkreactions(k);
    ssf(k) = MO.Reactions(vl).UserData;
end 

Fluxes = ssf;

%% Retrieve kinetic law from SBML file 
% All kinetics are written with Cytoplasm. or Out. as structure. 
% Eliminate structure from kinetic laws. 

for i = 1:1:numel(MO.Reactions)
    kl = MO.Reactions(i).ReactionRate;
        
    e_kl = erase(kl,'Cytoplasm.');
    e2_kl = erase(e_kl,'Out.');
    e3_kl = erase(e2_kl,'Growth_');
    
    MO.Reactions(i).ReactionRate = e3_kl;
       
end 

for i = 1:1:numel(ntwkreactions)
        idx = ntwkreactions(i);
        REQ{i} = MO.Reactions(idx).ReactionRate;
end 

%% Retrieve parameters
% in MO.Parameters are reported the Global quantity parameters. 
% MO.Reactions(11, 1).KineticLaw.Parameters
index = 1;

% For some reason growth rate parameters are not reported in the SBML model
% and therefore to be added manually to the parameter values matrix

for i=1:1:numel(ntwkreactions)
    intrac_reac = ntwkreactions(i);
    if ~isempty((MO.Reactions(intrac_reac, 1).KineticLaw))
    for j = 1:1:numel(MO.Reactions(intrac_reac, 1).KineticLaw.Parameters)
    ParametersID(index)     = {getfield(MO.Reactions(intrac_reac, 1).KineticLaw.Parameters(j),'Name')};
    ParameterValue(index)   = getfield(MO.Reactions(intrac_reac, 1).KineticLaw.Parameters(j),'Value');
    index = index + 1;    
    end 
    end 
end 


% Growth rate parameter names and values
PID_growth = {'mu_max'	'ks_GAP'	'ks_PEP'	'ks_PYR' 'n_gr'...
     	'ks_ATP'	'Ks_ACoA'	'Ks_AKG'	'Ks_OAA'	'Ks_NADH'...
        'ks_F6P'	'ks_G6P'	'ks_D3PG'	'ks_FORin'}; %'ks_FORin'};
 
 
 PVL_growth = [0.0003323	0.00215	0.0075	0.007	1.7  2.814	0.0043	0.00143	...
              0.009	0.00072	0.015	0.075	0.0018  0.0018];

ParametersID = [PID_growth  ParametersID(4:end )];

 ParameterValue = [PVL_growth ParameterValue(4:end)];

MO = MO;

end 




