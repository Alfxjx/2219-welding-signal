% emd 分解法使用


% 不需要 参数就是过零的变极性焊接
% min_curr = min(current);
% min_volt = min(volt);
% 
% for k=1:length(current)
%     current(k) = current(k)-min_curr;    
% end
% for m=1:length(volt)
%     volt(m) = volt(m)-min_volt;
% end

[imf_volt,residual_volt] = emd(volt);
[imf_curr,residual_curr] = emd(current);

volt_len = size(imf_volt,2);
curr_len = size(imf_curr,2);

% plot(imf_volt(:,4));

% 需要注意的是，下面的循环次数和产生的imf分量的个数相同，可能需要手动改一下

% 绘制电压的imf图
figure(1);
% for i=1:volt_len
%     subplot(volt_len,1,i);
%     plot(imf_volt(:,i));
%     title(['volt',num2str(i)]);
%     % subplot(10,4,(4*i-1));
%     % plot(fftshift(imf_volt(:,i)));
% end
% % 绘制对应imf的傅里叶频谱
% figure(2);
% for i=1:volt_len
%     subplot(volt_len,1,i);
%     plot(fftshift(imf_volt(:,i)));
%     title(['volt-fft',num2str(i)]);
%     % subplot(10,4,(4*i-1));
%     % plot(fftshift(imf_volt(:,i)));
% end
% % 电流的imf & fft
% figure(3);
for i=1:curr_len
    subplot(curr_len,1,i);
    plot(imf_curr(:,i));
    title(['curr',num2str(i)]);
end    
% figure(4);
% for i=1:curr_len
%     subplot(curr_len,1,i);
%     plot(fftshift(imf_curr(:,i)));
%     title(['curr-fft',num2str(i)]);
% end  

% figure(5);
% plot(volt);
% title('volt');
% figure(6);
% plot(current);
% title('current');

% hht transform
% figure(7);
% for i=1:volt_len
%     subplot(volt_len,1,i);
%     plot(hht(imf_volt(:,i)));
%     title(['hht-volt', num2str(i)]);
% end
% figure(8);
% for i=1:curr_len
%     subplot(curr_len,1,i);
%     plot(hht(imf_curr(:,i)));
%     title(['hht-curr', num2str(i)]);
% end


