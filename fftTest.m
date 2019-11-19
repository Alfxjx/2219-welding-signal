% 加载数据
data = importdata('1.txt',',',0);
% 电流电压

% 取一小段测试一下 需要对应
current = data(23000:36000,2);
volt = data(23000:36000,1);

% 预处理数据，将电压/电流信号转化成非负信号

minVIndex = find(volt==min(volt));

out = volt(minVIndex);

for i=1:len(volt)
    volt(i) = volt(i)+out;
end
for j=1: