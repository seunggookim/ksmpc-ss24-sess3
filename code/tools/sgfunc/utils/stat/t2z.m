function zstats = t2z(tstats, df)
% zstats = t2z(tstats, df)
zstats = norminv(tcdf(tstats, df), 0, 1);

end
