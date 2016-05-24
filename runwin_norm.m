function dataout = runwin_norm(datain)
	N = 5;
	smamp = smooth(abs(datain),N);
	dataout = datain(:)./smamp(:);
return
