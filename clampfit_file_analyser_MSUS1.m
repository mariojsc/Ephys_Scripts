%cd the Rec folder with animals
%the code is A# for animal, s# for slice, c# for cell, 'holding current'.
%txt for file

r=0;
r2=0;
fld1=dir('A*');

for i=1:length(fld1);
%go to animal
oldfolder=cd(fld1(i).name);

fld11=dir('s*'); 
    for j=1:length(fld11);
    %go to slice
    oldfolderr=cd(fld11(j).name);    
    
    fld111=dir('c*');
        for x=1:length(fld111);
        %go to cell
        oldfolderrr=cd(fld111(x).name);
        
        % Get all files in the current cell
        
        filesc= dir('*.atf');
        if length(filesc)>=1
            for id = 1:length(filesc)
                    % Get the file name (minus the extension)
                    [~, f] = fileparts(filesc(id).name);
                    % rename
                    movefile(filesc(id).name, sprintf('%s.txt',f));
            end
        end
        
        files = dir('-*.txt');
        if (length(files))~=0
            for id = 1:length(files);
                        % creates tables with values from the converted files
                        
                        B=readtable(files(id).name,'HeaderLines',2,'Delimiter','tab','ReadVariableNames',false);
                        %stores in struct
                        rec(id+r).type=files(id).name;
                        rec(id+r).values=B;
                        rec(id+r).animal=fld1(i).name;
                        rec(id+r).slice=fld11(j).name;
                        rec(id+r).cell=fld111(x).name; 
            end
            r=r+id;
        end
        
        files = dir('A*.txt');
        if (length(files))~=0
            for id = 1:length(files);
                        % creates tables with values from the converted files
                        
                        B=readtable(files(id).name,'HeaderLines',2,'Delimiter','tab','ReadVariableNames',false);
                        %stores in struct
                        recA(id+r2).type=files(id).name;
                        recA(id+r2).values=B;
                        recA(id+r2).animal=fld1(i).name;
                        recA(id+r2).slice=fld11(j).name;
                        recA(id+r2).cell=fld111(x).name; 
            end
            r2=r2+id;
        end
        
        cd(oldfolderrr)
        end
    cd(oldfolderr) 
    end
cd(oldfolder)
end

%%
    %take the size of the table from biggest size peak tables
    for i=1:length(rec)
        a(i)=(height(rec(i).values));
    end 
   
%% mini Analysis

types={'-70.txt','-55_E.txt','-55_I.txt','-70_ttx.txt'}; %defines columns of o (add here for more :))
variables={'Inst. Freq. (Hz)','Peak Amp (pA)','Interevent Interval (ms)'}; %defines rows of o

for j=1:length(types)
    c1=char(types(j));    
    for i=1:length(variables)
        o{i,j}=ephyst_extract(rec,c1,variables(i),a); %turns negative values to absolute
    end

end

%statistics will come here

%exports files and gives names to the csv
[n,m]=size(o);
for j=1:m
   for i=1:n
       T=cell2table(o{i,j});
       writetable(T,strcat(char(types(j)),'--',char(variables(i)),'.csv'));
   end
end

clear o
%% APs Analysis

types={'APs.txt'}; %defines columns of o (add here for more :))
variables={'Trace','Interevent Interval (ms)'};%,'Interevent Interval (ms)'}; %defines rows of o

%for the special case of trace
for j=1:length(types)
    c1=char(types(j));    
    for i=1:length(variables)
        o{i,j}=ephyst_extract(recA,c1,variables(i),a); %turn into cell to export later
    end
end

[n,m]=size(o);
% ot=zeros(n,m);

% for i=1:m
%     j=1;
%     jt=1;
%     while j<=n && jt<=n
%         if (o{j,i})>0
%         ot(jt,i)=ot(jt,i)+1;
%         j=j+1;
%         if j<=n
%             while o{j,i}==o{j-1,i} 
%             ot(jt,i)=ot(jt,i)+1;
%             j=j+1;
%             end
%         end
%         jt=jt+1;
%         end
%         j=j+1;
%     end
% end

% statistics will come here

% exports files and gives names to the csv

for j=1:m
   for i=1:n
       T=cell2table(o{i,j});    
       writetable(T,strcat(char(types(j)),'--',char(variables(i)),'.csv'));
   end
end