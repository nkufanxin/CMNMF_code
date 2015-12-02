clc;
clear;
fid = fopen( 'symbols_1274.txt', 'r' );
C = textscan(fid, '%s'); % 一行有两列数据
%data = cell2mat(C); % 转化成普通矩阵plot(data);hold on;