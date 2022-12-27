% % Modified by Mingxi Dang (18813059828@163.com) at 20221227
% Based on stepwise connectivity analysis: Sepulcre et al, Journal of Neuroscience 2012; Sepulcre, Cerebral Cortex 2015; Ortiz-Teran et al, PNAS 2017; Ibai & Sepulcre, Nat Commun 2018; Bueichekú et al, PNAS 2020.
clear;
clc;

q_fdr=0.00001; % q for the fdr 0.00001
num_steps=7; % number of steps for the SCA

%Read the PET matri
[num,txt] = xlsread('AD_Abeta_SUVR.xlsx');
PET_matrix = num(:,:);
seed_indx = 85;
seed_indx2 = 8;
%%Compute connectivity matrix
[fc,pvals]=corr(PET_matrix);
fc(1:(1+90):end)=0; %diagonal to zero
fc= 0.5 * log( ( 1 + fc)./( 1 - fc) ); %fisher (optional)

%%Take positive links and filter using fdr (all optional but recommended)
fc(find(fc<0))=0;
indx=find(fc>0);
[n_signif,indx_signif]=fdr(pvals(indx),q_fdr,'original','mean');
fdr_matrix=zeros(size(pvals));
fdr_matrix(indx(indx_signif))=1;
fdr_matrix=uint8(fdr_matrix);
adj=fc.*double(fdr_matrix);

%%Stepwise connectivity analysis 
adj=(adj - min(adj(:)))./ (max(adj(:)) - min(adj(:))); %normalize (optional)
[m n]=size(adj);
adj(1:(m+1):end)=0; %diagonal to zero
step1=adj;
degree_centrality=sum(step1,2);
adjp=step1;
step1_seed=adjp(seed_indx,:);
step1_seed2=adjp(seed_indx2,:);
for s=2:num_steps
  adj=adjp*adjp; % for exponential increase, multiply by output; for arithmetic increase multiply always by step1
  adj=(adj - min(adj(:)))./ (max(adj(:)) - min(adj(:))); %normalize (optional)
  adj(1:(m+1):end)=0; %diagonal to zero
  adjp=adj;
  eval(sprintf('step%s=adjp;',num2str(s)));
  eval(sprintf('step%s_seed=adjp(seed_indx,:);',num2str(s)));
%   eval(sprintf('step%s_seed2=adjp(seed_indx2,:);',num2str(s)));
end


%%Output result
%Transfer the label of the regions in the template to the value specified in excel(optional)
templatePath_1 = 'AAL_Contract_90_2MM.nii';
templateVol_1 = spm_read_vols(spm_vol(templatePath_1));
temp1 = spm_vol(templatePath_1);
tmpBinTemplate_1 = templateVol_1;
for i = 1:90
    if i == seed_indx
        tmpBinTemplate_1(tmpBinTemplate_1==i) = 1.3;
    else
    tmpBinTemplate_1(tmpBinTemplate_1==i) = step1_seed(i);
    end
end
tmpBinTemplate_1(tmpBinTemplate_1 ~= 0);
temp1.fname = fullfile('result','step1_seed.nii');
temp1.dt = [64, 0];
spm_write_vol(temp1,tmpBinTemplate_1);


templatePath_2 = 'AAL_Contract_90_2MM.nii';
templateVol_2 = spm_read_vols(spm_vol(templatePath_2));
temp2 = spm_vol(templatePath_2);
tmpBinTemplate_2 = templateVol_2;
for i = 1:90
    if i == seed_indx
        tmpBinTemplate_2(tmpBinTemplate_2==i) = 1.3;
    else
    tmpBinTemplate_2(tmpBinTemplate_2==i) = step2_seed(i);
    end
end
tmpBinTemplate_2(tmpBinTemplate_2 ~= 0);
temp2.fname = fullfile('result','step2_seed.nii');
temp2.dt = [64, 0];
spm_write_vol(temp2,tmpBinTemplate_2);

templatePath_3 = 'AAL_Contract_90_2MM.nii';
templateVol_3 = spm_read_vols(spm_vol(templatePath_3));
temp3 = spm_vol(templatePath_3);
tmpBinTemplate_3 = templateVol_3;
for i = 1:90
    if i == seed_indx
        tmpBinTemplate_3(tmpBinTemplate_3==i) = 1.3;
    else
    tmpBinTemplate_3(tmpBinTemplate_3==i) = step3_seed(i);
    end
end
tmpBinTemplate_3(tmpBinTemplate_3 ~= 0);
temp3.fname = fullfile('result','step3_seed.nii');
temp3.dt = [64, 0];
spm_write_vol(temp3,tmpBinTemplate_3);

templatePath_4 = 'AAL_Contract_90_2MM.nii';
templateVol_4 = spm_read_vols(spm_vol(templatePath_4));
temp4 = spm_vol(templatePath_4);
tmpBinTemplate_4 = templateVol_4;
for i = 1:90
    if i == seed_indx
        tmpBinTemplate_4(tmpBinTemplate_4==i) = 1.3;
    else
    tmpBinTemplate_4(tmpBinTemplate_4==i) = step4_seed(i);
    end
end
tmpBinTemplate_4(tmpBinTemplate_4 ~= 0);
temp4.fname = fullfile('result','step4_seed.nii');
temp4.dt = [64, 0];
spm_write_vol(temp4,tmpBinTemplate_4);

templatePath_5 = 'AAL_Contract_90_2MM.nii';
templateVol_5 = spm_read_vols(spm_vol(templatePath_5));
temp5 = spm_vol(templatePath_5);
tmpBinTemplate_5 = templateVol_5;
for i = 1:90
    if i == seed_indx
        tmpBinTemplate_5(tmpBinTemplate_5==i) = 1.3;
    else
    tmpBinTemplate_5(tmpBinTemplate_5==i) = step5_seed(i);
    end
end
tmpBinTemplate_5(tmpBinTemplate_5 ~= 0);
temp5.fname = fullfile('result','step5_seed.nii');
temp5.dt = [64, 0];
spm_write_vol(temp5,tmpBinTemplate_5);

templatePath_6 = 'AAL_Contract_90_2MM.nii';
templateVol_6 = spm_read_vols(spm_vol(templatePath_6));
temp6 = spm_vol(templatePath_6);
tmpBinTemplate_6 = templateVol_6;
for i = 1:90
    if i == seed_indx
        tmpBinTemplate_6(tmpBinTemplate_6==i) = 1.3;
    else
    tmpBinTemplate_6(tmpBinTemplate_6==i) = step6_seed(i);
    end
end
tmpBinTemplate_6(tmpBinTemplate_6 ~= 0);
temp6.fname = fullfile('result','step6_seed.nii');
temp6.dt = [64, 0];
spm_write_vol(temp6,tmpBinTemplate_6);

templatePath_7 = 'AAL_Contract_90_2MM.nii';
templateVol_7 = spm_read_vols(spm_vol(templatePath_7));
temp7 = spm_vol(templatePath_7);
tmpBinTemplate_7 = templateVol_7;
for i = 1:90
    if i == seed_indx
        tmpBinTemplate_7(tmpBinTemplate_7==i) = 1.3;
    else
    tmpBinTemplate_7(tmpBinTemplate_7==i) = step7_seed(i);
    end
end
tmpBinTemplate_7(tmpBinTemplate_7 ~= 0);
temp7.fname = fullfile('result','step7_seed.nii');
temp7.dt = [64, 0];
spm_write_vol(temp7,tmpBinTemplate_7);

