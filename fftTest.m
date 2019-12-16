for i=1:size(imf_curr,2)
    n = size(imf_curr,1);
    y = fft(imf_curr(:,i));
    yshift = fftshift(y);
    power = abs(yshift).^2/n;
%     f = (0:length(y)-1)*100000/length(y);
    f = (-n/2:n/2-1)*(100000/n);
    figure(i);
    plot(f,power);
    xlim([-3e4,3e4]);
%     plot(f,abs(y));
end    