% 时频熵计算
% 计算HHT谱 然后对二维平面进行分块，
% 接着对分块的能量进行归一化
% 计算信息谱的时频熵

% 二维平面 hs
cell = mat2cell(hs,[10,10,10,10,10,10,10,10,10,11],[4096,4096,4096,4096,4096,4096,4096,4096,4096,4097]);
energy_of_each = cellfun(@(x) sum(full(x(:))), cell);
% 接着对分块的能量进行归一化
energy_of_all = sum(energy_of_each(:));
en_q = energy_of_each/energy_of_all;
% 时频熵
res_en = 0;
size_en = size(energy_of_each);

for i=1:size_en(1)
    for j=1:size_en(2)
        single_en = -1*en_q(i,j)* log(en_q(i,j));
        res_en = res_en+single_en;
    end
end
