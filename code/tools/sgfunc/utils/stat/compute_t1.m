function t = compute_t1(x)
dim = 1;
samplesize = size(x,dim);
xmean = mean(x,dim);
sdpop = std(x,[],dim);
sqrtn = sqrt(samplesize);
ser = sdpop ./ sqrtn;
t = xmean ./ (std(x,[],1) ./ sqrt(n);
end
