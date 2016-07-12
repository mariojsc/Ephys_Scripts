%cd the Rec folder with animals
%the code is A# for animal, s# for slice, c# for cell, 'holding current'.
%txt for file

%there are Nan values on top of values in value cells of o and
%such. This is the result of ephyst_extract starting applying str2double
%from the very top of the files, so variable name gets Nan, and for instant
%frequency, also the first N/A

cd

%% import recording data sets

r=0;
r1=0;
r2=0;
r3=0;
r4=0;
fld1=dir('A*');

for i=1:length(fld1);
%go to animal
oldfolder=cd(fld1(i).name);

fld11=dir('s*'); 
    
    for j=1:length(fld11);
        %% 
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
        
        filesc= dir('-55_bl_results_E.txt'); %change name of this type of file
        if length(filesc)>=1
            for id = 1:length(filesc)
                    % Get the file name (minus the extension)
                    [~, f] = fileparts(filesc(id).name);
                    % rename
                    movefile(filesc(id).name,'-55_E.txt');
            end
        end
        
        files = dir('-55_E.txt');
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
        
        files = dir('*.png');
        if (length(files))~=0
            for id = 1:length(files);
                
                        % stores in struct
                        recD(id+r4).type=files(id).name;
                        recD(id+r4).values=files(id).date;
                        recD(id+r4).animal=fld1(i).name;
                        recD(id+r4).slice=fld11(j).name;
                        recD(id+r4).cell=fld111(x).name;
            end
            r4=r4+id;
        end
        
        files = dir('CC1_curve.txt');
        if (length(files))~=0
            for id = 1:length(files);
                        % creates tables with values from the converted files
                        
                        B=readtable('CC1_curve.txt','HeaderLines',2,'Delimiter','tab','ReadVariableNames',false);
                        %stores in struct
                        recA(id+r1).type=files(id).name;
                        recA(id+r1).values=B;
                        recA(id+r1).animal=fld1(i).name;
                        recA(id+r1).slice=fld11(j).name;
                        recA(id+r1).cell=fld111(x).name; 
            end
            r1=r1+id;
        end
        
        files = dir('CC.txt');
        if (length(files))~=0
            for id = 1:length(files);
                        % creates tables with values from the converted files
                        
                        B=readtable('CC.txt','HeaderLines',0,'Delimiter','tab','ReadVariableNames',false);
                        %stores in struct
                        CC(id+r2).type=files(id).name;
                        CC(id+r2).values=B;
                        CC(id+r2).animal=fld1(i).name;
                        CC(id+r2).slice=fld11(j).name;
                        CC(id+r2).cell=fld111(x).name;
                        
            end
            r2=r2+id;
        end
        
        files = dir('CC1_spiket.txt');
            if (length(files))~=0
                for id = 1:length(files);
                            % creates tables with values from the converted files

                            B=readtable(files(id).name,'HeaderLines',2,'Delimiter','tab','ReadVariableNames',false);
                            %stores in struct
                            recB(id+r3).type=files(id).name;
                            recB(id+r3).values=B;
                            recB(id+r3).animal=fld1(i).name;
                            recB(id+r3).slice=fld11(j).name;
                            recB(id+r3).cell=fld111(x).name; 
                end
                r3=r3+id;
            end
        cd(oldfolderrr)
        end
    cd(oldfolderr) 
    end
cd(oldfolder)
end

% take rec structure and generate a blinding.txt with fields for each
% parameter. This text file has to be filed in order to separate groups and
% perform statistics

clear cache
cache=dir;

v=0;
for i=1:length(cache)
    ccache=cache(i).name;
        if strcmp(ccache,'~blinding.txt')
        v=v+1;
        end
end

if v==0

    
fID=fopen('~blinding.txt','wt');
fprintf(fID,'type animal slice cell pharmacology experimental_group\n');

for i=1:length(rec)
    type=rec(i).type;
    animal=rec(i).animal;
    slice=rec(i).slice;
    cell=rec(i).cell;
    fprintf(fID,'%s %s %s %s\n',type,animal,slice,cell);
end
    fclose(fID);
    
end


