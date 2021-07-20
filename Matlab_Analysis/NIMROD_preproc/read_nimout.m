format longE

%% reading in data

filein='/rmp_m3/nimrod.out';
fileID = fopen(filein);
strings = textscan(fileID,'%s');
fclose(fileID);

strings=strings{1};

%% finding all time steps and recording timestep number, time, dt

times=[];
jj=2;
for ii=1:size(strings,1)
    if strcmp(strings{ii},'Time')
       times(1,jj)=str2num(strings{ii-1});
       times(2,jj)=str2num(strings{ii+2});
       times(3,jj)=str2num(strings{ii+5});
       
       jj=jj+1;
    end
end

times(1,1)=times(1,2)-1;
times(2,1)=times(2,2)-times(3,2);

%% finding cycles and times corresponding to dump files

dumpit=50;

times1={};
jj=2;
times1{1,1}=strcat('dump.',num2str(times(1,1)));
times1{2,1}=times(2,1);
for ii=2:size(times,2)
   if mod((times(1,ii)-times(1,1)),dumpit)==0
       times1{1,jj}=strcat('dump.',num2str(times(1,ii)));
       times1{2,jj}=times(2,ii);       
       jj=jj+1;
   end
end

%% writing to file

ifwrite=1;
if ifwrite==1
    
    fileout='dumpfiles_and_times_rmpm3.txt';
    fileID = fopen(fileout,'w');
    fprintf(fileID,'%i \n',size(times1,2));
    for ii=1:size(times1,2)
        fprintf(fileID,'%15s %11.10f\n',times1{1,ii},times1{2,ii});
    end
    fclose(fileID);
        
end
