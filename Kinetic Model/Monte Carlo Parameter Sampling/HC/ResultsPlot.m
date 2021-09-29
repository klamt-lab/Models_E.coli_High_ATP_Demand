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
%% Plot instograms of GCC
% First of all normalize all FCC between -1 1. 

close all
clear all
load('HC_MonteCarlo_result.mat')


% load Model structure 
[Fluxes, Conc, N, N_red, L, S, R, ParametersID,ParameterValue,MO,REQ,ntwkreactions,intracellular_index,Cnc] = ReadSBML()

%%
% Check that FCC and CCC are wright. 
% Summation theorem
% sum(FCC) = 1;
% sum(CCC) = 0;
for i=1:1:numel(CS_rec(1,1,:))
    sm(i) = sum(sum(CS_rec(:,:,i)));
end 

for i=1:1:numel(CJ_rec(1,1,:))
    sj(i) = sum(sum(CJ_rec(:,:,i)));
end 

%%
CJ_norm = CJ_rec;
index = 1;

%%
% Selected vector index to avoid ploting transport reactions for LAC ETH ACE SUC and FOR 
index_r_ex = 1;
ing = 1;
for i=1:1:length(R)
    cnr = R(i);
    if  contains(cnr,'trsp') || contains(cnr,'SUCtr')
        g(ing) = i;
        ing = ing +1;
    else 
        ntwkreactionsNoTrsp(index_r_ex) = ntwkreactions(i);
        tpntrsp(index_r_ex) = i;
        index_r_ex = index_r_ex +1;
       
    end 
end 

nsubplots = numel(ntwkreactionsNoTrsp);

%%
%% FCC values calculated for WT or HC with the deterministic model 
% To plot values of FCC calculates with copasi for the model described in
% the documentation
P=readtable(fullfile('MCA_FCC_HC'));

fn = [fieldnames(P)];
fn = regexprep(fn,'^x','','emptymatch');
fn = regexprep(fn,'^_','','emptymatch');
fn = regexprep(fn,'_$','','emptymatch');

% Identify the rows and column in the MCA to be plot   
for i = 1:1:nsubplots    
    ii = ntwkreactionsNoTrsp(i);
    nametofind = {MO.Reactions(ii).Name};
    for z = 1:1:numel(fn)
        fntc = fn(z);
        if  strcmp(fntc,nametofind)
            mca_indx(i) = z;
        end 
    end     
end

%% Plot instogram

F = figure(1);
for i = 1:1:nsubplots % effect 
    
    ii = tpntrsp(i);
    iii = mca_indx(i) -1;
        
    for j=1:1:nsubplots  % cause 
       
        jj = tpntrsp(j);
        jjj = mca_indx(j);
        
        X = CJ_norm(ii,jj,:);           % Plot all columns and then the rows
        
        % generate matrix with deterministic FCC of the reactions to plot
        FCCfx(i,j) = table2array(P(iii,jjj));
        FCCtoplot = (FCCfx(i,j));
        
        % Compute standard deviation
        SDsp(i,j) = std(X);
        MVsp(i,j) = mean(X);
        Mednsp(i,j) = median(X);
           
        cdx = std(X);   
        
        subplot(nsubplots,nsubplots,index)
        
        IS = histogram(X,'FaceColor','k');
%        % draw red lines for deterministic FCC
%         hold on
%         xline(FCCtoplot, 'r', 'LineWidth',1.5)
        IS.BinWidth = 0.01;
        
        xlim([-1 1])

        index= index+1;
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        
        if jj == 1
           y_effect = ntwkreactionsNoTrsp(i);
           yl = ylabel(MO.Reactions(y_effect).Name,'fontweight','bold');
           set(yl,'rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
           
        end
        if i == nsubplots
           x_cause  = ntwkreactionsNoTrsp(j);
           xlabel(MO.Reactions(x_cause).Name,'fontweight','bold')
        end            
        
%         

           
    end   
    
    rraction_names(i) = {MO.Reactions(ii).Name};
    
end

%% calculate average absolute standard deviation
smv = numel(MVsp(1,:));
MV1 = abs(SDsp);
MV2 = sum(MV1);
MV3 = sum(MV2);
MV4 = MV3 / (smv * smv);

%% save plots in A4 format
jFrame = get(handle(F), 'JavaFrame');

jFrame.setMaximized(1);
set(F, 'PaperSize',[42 23]);

title_plot = strcat('MV2_MV_HC','.svg');
saveas(F, title_plot,'svg'); 

title_plot = strcat('MV2_MV_HC','.pdf');
saveas(F, title_plot,'pdf'); 

%% 



% Put in a single array all metrics Standard Deviation, Median and Average
% Set the values as rows. 
indx   = 1;
ncol   = numel(rraction_names);
totcol = ncol +2;

% stup = {{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'}}
stup = {'-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'};
  % ,'-','-','-','-','-'

for j = 1:1:numel(SDsp(:,1)) % rows
    
    jj = tpntrsp(j); 
    rn = {MO.Reactions(jj).Name};
        % median
        % standard deviation
        % average
        % try to proceede for rows
        medianRow  = Mednsp(j,:);
        SD_Row     = SDsp(j,:);
        AverageRow = MVsp(j,:);

         T1(indx,1) = cell2table(rn);
         T1(indx,2) = cell2table({'Median'});
         T1(indx,3:totcol) = array2table(medianRow);

         indx = indx+1;
         T1(indx,1) = cell2table(rraction_names(j));
         T1(indx,2) = cell2table({'StdDev'});
         T1(indx,3:totcol) = array2table(SD_Row);
      
         indx = indx+1;
         T1(indx,1) = cell2table(rraction_names(j));
         T1(indx,2) = cell2table({'Average'});
         T1(indx,3:totcol) = array2table(AverageRow);
      
         indx = indx+1;
         
         % Insert a blank line [NaN] between each reactions
         T1(indx,1) = cell2table({'-'});
         T1(indx,2) = cell2table({'-'});

         indx = indx+1;

end 
   
r1 = {'Reaction'};
r2 = {'Metrics'};
r3 = rraction_names(:);
T1.Properties.VariableNames(:) = [r1;r2;r3];
filename = 'Metrics_n_half_10000s.xlsx';
writetable(T1,filename,'Sheet',1,'WriteVariableNames',true);      
        
