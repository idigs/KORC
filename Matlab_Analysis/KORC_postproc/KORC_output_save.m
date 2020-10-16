%% Reading KORC output files

tsF=20000;

readbegin=cputime;
ST=diagnoseKORC('../OUT_00/',false,[0,tsF]);
readone=cputime;

FOread=readone-readbegin;

fprintf('Data read in %6.5e s \n',readone-readbegin)

save('ST','ST','-v7.3')
