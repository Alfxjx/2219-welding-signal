% emd �ֽⷨʹ��


% ����Ҫ �������ǹ���ı伫�Ժ���
% min_curr = min(current);
% min_volt = min(volt);
% 
% for k=1:length(current)
%     current(k) = current(k)-min_curr;    
% end
% for m=1:length(volt)
%     volt(m) = volt(m)-min_volt;
% end
% 
% [imf_volt,residual_volt] = emd(volt);
% [imf_curr,residual_curr] = emd(current);

% imf_volt = eemd(volt,0.2,20);
imf_curr = eemd(current,0,1);

% volt_len = size(imf_volt,2);
curr_len = size(imf_curr,2);

% plot(imf_volt(:,4));

% ��Ҫע����ǣ������ѭ�������Ͳ�����imf�����ĸ�����ͬ��������Ҫ�ֶ���һ��

% ���Ƶ�ѹ��imfͼ
figure(1);
% for i=1:volt_len
%     subplot(volt_len,1,i);
%     plot(imf_volt(:,i));
%     title(['volt',num2str(i)]);
%     % subplot(10,4,(4*i-1));
%     % plot(fftshift(imf_volt(:,i)));
% end
% % ���ƶ�Ӧimf�ĸ���ҶƵ��
% figure(2);
% for i=1:volt_len
%     subplot(volt_len,1,i);
%     plot(fftshift(imf_volt(:,i)));
%     title(['volt-fft',num2str(i)]);
%     % subplot(10,4,(4*i-1));
%     % plot(fftshift(imf_volt(:,i)));
% end
% % ������imf & fft
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


