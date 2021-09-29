%% ------------------------------------------------------------------------
% AUTHOR: GIULIA SLAVIERO 
% SEPTEMBER 2021

% THIS CODE IS PART OF THE FOLLOWING PUBBLICATION 
% Boecker S., Slaviero G., Schramm T., Steuer R., Szymanski W., Link H.,
% Klamt S., (2021), "Deciphering the Physiological Response 
% of Escherichia coli Under High ATP Demand", Submitted     

% MORE INFORMATION ON THE MODEL CAN BE FOUND AT 
% https://github.com/klamt-lab/Models_E.coli_High_ATP_Demand

% PART OF THE CODE IS TAKEN FROM 
% Murabito E, Verma M, Bekker M, Bellomo D, Westerhoff HV, et al. (2014) 
% Monte-Carlo Modeling of the Central Carbon Metabolism of Lactococcus lactis: 
% Insights into Metabolic Regulation. PLOS ONE 9(9): e106453. 
% https://doi.org/10.1371/journal.pone.0106453
% -------------------------------------------------------------------------
%%
close all
clear all 

% import SBML data 
[Fluxes, Conc, N, N_red, L, S, R, ParametersID,ParameterValue,MO,REQ,ntwkreactions,intracellular_index,Cnc] = ReadSBML()


%% Define values of dependent model quantities 
       
  totAxP  = 2.70058;
  totNAD  = 2;
  totCoA  = 3.67;
  CO2  = 2.478;
  Cytoplasm = 1;
  Out = 1;  
  
%% Declaring variables 


iterations = 10000;                     % No. of sampling iterations

nOfSS = 0;                              % Initialize number of stable steady states

F1 = 0.1;                               % Multiplicative factor defining the lower bound of the sampling intervals
F2 = 10;                                % Multiplicative factor defining the upper bound of the sampling intervals


MaxRealEigens = zeros(iterations,1);    % At iteration i, the maximal real
                                        % part of the eigenvalues is stored in
                                        % MaxRealEigens(i,1)
                                     
% Rename variables
FullSto = N;                            % Full stoichtiometric Matrix
RedSto  = N_red;                        % Reduced Stoichiometric matrix
Link    = L;                            % Link matrix 


[m,n] = size(FullSto);                  % Size of the full stoichiometric matrix
[m1,n] = size(RedSto);                  % Size of the reduced stoichiometric matrix

CS = zeros(m,n);                        % Initialize Concentation Control Coefficients
CJ = zeros(n,n);                        % Initialize Flux Control Coefficients
                 
CJ_rec = zeros(n,n,iterations);         % Flux control coefficients of stable systems
CS_rec = zeros(m,n,iterations);         % Concentration control coefficients of stable systems
E_rec = zeros(n,m,iterations);          % Elasticity coefficients of stable systems

dfodc = zeros(n,m);                     % Matrix containing the numeric value of the derivative
                                        % of the fluxes with respect to the metabolite
                                        % concentrations.

RedJac1 = zeros(m1,m1);                 % Reduced Jacobian Matrix

ID = eye(n);                            % Unitary matrix with n. of columns/rows = n. of reactions

Parameters = zeros(iterations,length(ParametersID)); % At iteration i, all the
                                              % parameters (sampled and 
                                              % not) are stored in row i
                                              % of this matrix.

%% Compute notrmalization matrices 
Dj = diag(Fluxes);                      % Matrix with steady state fluxes on the diagonal
Dj_inv = inv(Dj);                       
Ds = diag(Conc);                        % Matrix with steady state concnetrations on the diagonal
Ds_inv = inv(Ds);
nrmlz3 = Conc' * ((1./Fluxes));

%% -- RETRIEVE THE FUNCTIONAL FORM OF THE RATE EQUATIONS --%
tic
DER = compSymDeriv(R,S,ParametersID,REQ,MO); % Compute simbolically the symbolic partial derivatives of
                                             % the fluxes with respect to the concentrations       
toc

%-- COMPUTE THE SAMPLING INTERVAL OF EACH PARAMETER --%
%-- ACCORDING TO THE PARAMETER ALREADY PRESENT IN THE MODEL  --%

[ParValMin,ParValMax,vmax_indeces] = compParInt(S,F1,F2,ParameterValue,ParametersID)

fprintf('%s\n', 'Done');

eval('default = 1;');

tic

if ~exist('seed','var'); seed = 1; 
end
rng(seed);


tic
for c = 1:iterations
    
    %-- SAMPLE/EVALUATE PARAMETERS --%
    for p=1:length(ParametersID)        
        ParAux = log(ParValMin(p)) + rand(1)*(log(ParValMax(p)) - log(ParValMin(p)));
        ParAux = exp(ParAux);
        eval([ParametersID{p},' = ParAux;']);
        Parameters(c,p) = ParAux;
    end
    
    %-------------------------------------------------------------
    %-- DETERMINING THE PARTIAL DERIVATIVES FOR STATE 1 --%
    
    %--  Evaluating concentration values  --%
    for i = 1:length(S)
        eval([S{i},' = Cnc(i,:).Number;']);
    end   
    
    
     AMP  = totAxP - ATP - ADP;
     NAD  = totNAD - NADH;
     CoA  = totCoA - ACoA;
    
    %-- Computing Vmaxs --%
    for j = 1:length(vmax_indeces)
        eval([ParametersID{vmax_indeces(j)},' = 1;']);
    end 
    for j=1:numel(vmax_indeces)
        pippo = eval([REQ{j},';']);
        Fluxes(j);
        Vmax = Fluxes(j)/pippo;
        eval([ParametersID{vmax_indeces(j)},' = Vmax;']);
        Parameters(c,vmax_indeces(j)) = Vmax;
    end
    
    
    %-- Evaluating numerical value of partial derivatives --%
    for i=1:size(DER,1)
        for j=1:size(DER,2)
            dfodc(j,i) = eval([DER{i,j}]);
            if (~isreal(dfodc(j,i)))
              dfodc(j,i) = 0;
            elseif (isnan(dfodc(j,i)))
                dfodc(j,i) = 0;
            end 
        end
    end
      
    % ---- COMPUTING THE REDUCED JACOBIANS ----%
    RedJac1(:,:) = RedSto(:,:)*dfodc(:,:)*Link(:,:);
      
    %---- COMPUTING THE EIGENVALUES ----%
    eigenvalues = eig(RedJac1);
    MaxRealEigen = max(real(eigenvalues));
    MaxRealEigens(c,1) = MaxRealEigen;
    %---- CHECKING FOR STABILITY ----%
    ok = 1;
    if(MaxRealEigen >= 0)
        ok = 0;
    end

    if(ok == 1)
        nOfSS = nOfSS + 1;
        
    %---- COMPUTING AND RECORDING ELASTICITY COEFFICIENTS ----%
    E_rec(:,:,nOfSS) = dfodc .* (nrmlz3(:,:))';
    %---- COMPUTING CONCENTRATION and FLUX CONTROL COEFFICIENTS ----%
    CS_f(:,:) = -( Link*(inv(RedJac1(:,:)))*RedSto);
    CJ_f(:,:) = ( dfodc(:,:)*CS_f(:,:) + ID);
    CS(:,:) = Ds_inv * CS_f * Dj;
    CJ(:,:) = Dj_inv * CJ_f * Dj;
        
    CJ_rec(:,:,nOfSS) = CJ;
    CS_rec(:,:,nOfSS) = CS;
    end
    
    ratio = (nOfSS/c)*100;
    if(~mod(c,10) || c == iterations)
    	fprintf('%d%s%d%s%6.2f%s\n',nOfSS,' -> ',c,' (',ratio,'%)');
    end

end 
toc

Pcontrol = [ParametersID; string(ParameterValue)];  % generate matrix to control sampled values 

CJ_rec(:,:,(nOfSS+1):end) = [];
CS_rec(:,:,(nOfSS+1):end) = [];
E_rec(:,:,(nOfSS+1):end) = [];
numSim = nOfSS;

% Save result
 save('HC_MonteCarlo_result', 'CJ_rec','CS_rec','E_rec','iterations','MaxRealEigens','Parameters')


%%
%==========================================================================
function [ParValMin,ParValMax,vmax_indeces] = compParInt(S,F1,F2,ParameterValue,ParametersID)

    fprintf('%s', 'Evaluating parameter sampling intervals... ');

    ParValMin = ParameterValue;
    ParValMax = ParameterValue;

    count=1;
    for p = 1:length(ParametersID)

        par_aux = ParametersID{p};
        % exclude some parameters from sampling procedure. Hill
        % coefficients are not sampled
        if ((~contains(par_aux,'max')) ...
                && (~contains(par_aux,'_eq'))  ...
                && (~contains(par_aux,'_EQ')) ...    
                && (~contains(par_aux,'n_'))...
                && (~contains(par_aux,'GLC_ex'))...
                && (~contains(par_aux,'ks_'))...
                && (~contains(par_aux,'Ks_'))...
                && (~contains(par_aux,'Ph'))...
                && (~contains(par_aux,'mu'))... 
                && (~contains(par_aux,'mu_max'))) ...
                && (~contains(par_aux,'CO2'))...
    
              ParValMin(p) = ParameterValue(p)*F1;
              ParValMax(p) = ParameterValue(p)*F2;
        end

        % Set Vmax =1 and register their position
        if contains(par_aux,'max') 
            ParValMin(p) = 1;
            ParValMax(p) = 1;    
            vmax_indeces(count) = p;
            count = count+1;      
        end
    end
end