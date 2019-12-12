function psnr = getPSNR(noisy, clean, maxVal)
diff = noisy - clean;
rmse = sqrt(mean(diff(:).^2));
psnr = 20*log10(maxVal/rmse);
end