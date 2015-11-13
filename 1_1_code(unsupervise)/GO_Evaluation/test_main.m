clc;
clear all;
load('../../5_6_result_2015/CMNMF/CMNMF_LF/predicted_pathway_LF_20151029174257.mat');
ratio = get_Ratio(pathway_gene_mgi_id{1});
disp('..........');