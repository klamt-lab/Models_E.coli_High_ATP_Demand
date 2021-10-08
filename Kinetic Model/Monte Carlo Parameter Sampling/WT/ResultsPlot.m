%% ------------------------------------------------------------------------
% AUTHOR: GIULIA SLAVIERO 
% SEPTEMBER 2021

% THIS CODE IS PART OF THE FOLLOWING PUBBLICATION 
% Boecker S., Slaviero G., Schramm T., Szymanski W., Steuer R., Link H.,
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
load('WT_MonteCarlo_result.mat')

%%
% load excel table with FCC computed with the parameter set currently
% present in the model 
P=readtable(fullfile('MCA_FCC_WT'));

%% 
% load Model structure 
[Fluxes, Conc, N, N_red, L, S, R, ParametersID,ParameterValue,MO,REQ,ntwkreactions,intracellular_index,Cnc] = ReadSBML();

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
% select reactions to plot. Discard dilution and transport reactions
index_r_in = 1;  
index_r_ex = 1;
for i=1:1:length(MO.Reactions)
    cnr = MO.Reactions(i).Name;
    if  contains(cnr,'exc')|| contains(cnr,'D_') || contains(cnr,'trsp') || contains(cnr,'SUCtr')
        suppreactions(index_r_in) = i;
        index_r_in = index_r_in+1;
    else 
        ntwkreactionsN(index_r_ex) = i;
        index_r_ex = index_r_ex +1;
    end 
end 
%%
CJ_norm = CJ_rec;
index = 1; 

tpntrsp = ntwkreactionsN;
nsubplots = numel(tpntrsp);
%%
fn = [fieldnames(P)];
fn = regexprep(fn,'^x','','emptymatch');
fn = regexprep(fn,'^_','','emptymatch');
fn = regexprep(fn,'_$','','emptymatch');

% Identify the rows and column in the MCA to be plot   
for i = 1:1:nsubplots    
    nametofind = {MO.Reactions(i).Name};
    for z = 1:1:numel(fn)
        fntc = fn(z);
        if  strcmp(fntc,nametofind)
            mca_indx(i) = z;
        end 
    end    
end 

Parr = table2array(P(:,2:end-1));
%% Plot instogram

F = figure(1);
for i = 1:1:nsubplots % effect 
    
    ii = tpntrsp(i);
        
    for j=1:1:nsubplots  % cause 
       
        jj = tpntrsp(j);
        
        X = CJ_norm(ii,jj,:);           % Plot all columns and then the rows
        FCCtoplot = (Parr(ii,jj));
        
        % Compute standard deviation
        SDsp(i,j) = std(X);
        MVsp(i,j) = mean(X);
        Mednsp(i,j) = median(X);
           
        cdx = std(X);   
        
        subplot(nsubplots,nsubplots,index)
        IS = histogram(X,'FaceColor','k');
       % draw red lines for deterministic FCC
        hold on
        xline(FCCtoplot, 'r', 'LineWidth',1.5)
        
        IS.BinWidth = 0.01;
        
        if std(X) > 340
            xlim([-10 10]) 
            set(gca,'color',[222, 226, 230]/255)
        else   
            xlim([-1 1])
        end 
        

        index= index+1;
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        
        if jj == 1
           y_effect = ntwkreactionsN(i);
           yl = ylabel(MO.Reactions(y_effect).Name,'fontweight','bold');
           set(yl,'rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
           
        end
        if i == nsubplots
           x_cause  = ntwkreactionsN(j);
           xlabel(MO.Reactions(x_cause).Name,'fontweight','bold')
        end              
    end   
    
    rraction_names(i) = {MO.Reactions(ii).Name};
    
end


%% save plots in A4 format
jFrame = get(handle(F), 'JavaFrame');

jFrame.setMaximized(1);
set(F, 'PaperSize',[42 23]);

title_plot = strcat('MV2_MV_WT','.svg');
saveas(F, title_plot,'svg'); 

title_plot = strcat('MV2_MV_WT','.pdf');
saveas(F, title_plot,'pdf'); 

