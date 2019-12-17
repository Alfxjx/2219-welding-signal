% 入口文件，整个流程依托 matlab 的工作区完成，在本文件夹下存放按逗号分隔的数据。
% 数据类型是第一列电压，第二列电流
% 首先导入数据 分别保存为 dadi_ 和 gm_ 开头的变量。
importData

% 对数据进行一个简单的观察
showbrief

% 设置取样位置和取样的长度
start1 = input('input start dadi: ');
step1 = input('input step dadi: ');
% start2 = input('input start gm: ');
% step2 = input('input step gm: ');
% 显示一下待分析的小波
dadi_volt = dadi(start1:(start1+step1),1);
dadi_curr = dadi(start1:(start1+step1),2);

% gm_volt = gm(start2:(start2+step2),1);
% gm_curr = gm(start2:(start2+step2),2);

showall

% 对数据进行 emd 分析
% 对所得的IMF研究其频谱
volt = dadi_volt;
current = dadi_curr;
emdTest
fftTest
% 计算对应的IMF的 HHT 变换
hhtTest

% 切换一下注释

% volt = gm_volt;
% current = gm_curr;
% emdTest
% % 计算对应的IMF的 HHT 变换
% hhtTest

% 计算时频熵
tfEN


