% Output wave-xxx.dat files and this program must be in the same directory.

clear
clc

dir = ls;
[a,~] = size(dir);
strdir = string(zeros(a,1));
for n = (1:a)
    strdir(n,:) = strtrim(convertCharsToStrings(dir(n,:)));
end
clear a dir
dat = strdir(contains(strdir,".dat"));
clear strdir
if(~isempty(dat))
    noFile = length(dat);
else
    noFile = 0;
end
clear dat

opt = fixedWidthImportOptions("NumVariables",2,"VariableWidths",[15,15],...
"VariableTypes",{'double','double'});

data = zeros(noFile,401,2);

for n =[1:noFile]
    no = sprintf('%03d',n-1);
    filename = append("wave-",no,".dat");
    data(n,:,:) = readmatrix(filename,opt);
end

for n=[1:noFile]
   no = sprintf('%03d',n-1);
   filename = append("wave-",no,".dat");
   plot(data(n,:,1),data(n,:,2));
   title(filename);
   xlim([-20,20]);
   ylim([-3,3]);
   pause(0.05);
end