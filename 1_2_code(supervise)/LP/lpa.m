%lpa
function [pathway_matrix, rmse,f_u] = lpa(g_g_sim, pathway_matrix, target, iter)

for i = 1:size(g_g_sim)
 g_g_sim(i,:) = g_g_sim(i,:)/sum(g_g_sim(i,:));
end

f_u = rand(1,size(pathway_matrix,2));
f_u0 = pathway_matrix(target, :);
rmse = zeros(1, iter);
y_l = pathway_matrix;
y_l(target, :) = [];
p_ul = g_g_sim(target, :);
p_ul(target) = [];
%p_lu = g_g_sim(1:-2, -1);
p_uu = g_g_sim(target, target);
for i = 1:iter
    temp=f_u;
    f_u = p_uu*f_u + p_ul*y_l;
    rmse(i) = sum((f_u - temp).^2)/size(f_u, 2); 
    disp(rmse(i));
end
%I = eyes(1);
%f_u = (I - p_uu)^-1*p_ul*y_l;
[~,index]=sort(f_u,'descend');

f_u = zeros(1,size(f_u,2));
f_u(index(1:25)) = 1;
pathway_matrix(target, :) = f_u;
end