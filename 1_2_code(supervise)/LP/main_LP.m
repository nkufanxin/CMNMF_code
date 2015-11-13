
path(path,'../../3_usefuldata');
load('g_p_netwrok_2015_3_4&5.mat');
load('similarity.mat');

%max_pathway_gene_num为筛选出G0的阈值（pathway的长度最大值）
%T为筛选出G0_NolessThan_T的阈值（每个pathway基因与mgi_id相交后pathway的最小长度）
max_pathway_gene_num=200;
T=5;
iter = 300;
build_G0(max_pathway_gene_num,T);
load('G0_data.mat');
load('mgi_id');

%t_circle为初始化次数
t_circle=5;
          
 %roc_n_set为计算roc时的传参值n
 roc_n_set = 168;
 %roc_gene为存放结果的向量，每一个元素代表去掉第i个基因训练后的auc结果值
 %roc_gene=zeros(length(roc_n_set),length(SelectedGene));
 roc_gene=zeros(length(SelectedGene),t_circle);
 %每一个roc_n_set对应一个roc cell，每个cell为一个矩阵，存储不同参数下，不同初始结果的auc值
 roc = zeros(1,t_circle);
 for t=1:t_circle                  
      for i=1:length(SelectedGene)
            %target为需要去掉的基因所对应的mgi_id位置
            target=find(mgi_id==SelectedGene(i));
            
            [W_out,rmse,f_u] = lpa(simi,G0_NoLessThan_T_norm,target,iter);

            ranked_pathways_for_specific_gene=rankgene(W_out,G0_NoLessThan_T,target);


            roc_gene(i,t)=ROC(ranked_pathways_for_specific_gene,roc_n_set);  

      end
       roc(t)=mean(roc_gene(:,t));  
 end

datetime=fix(clock);
      s='';
      for i=1:6
      s=[s num2str(datetime(i))];
      end
      fn2=[ '../../5_8_result_2015/LP/ROC_' s '.mat'];
      save(fn2,'roc');

