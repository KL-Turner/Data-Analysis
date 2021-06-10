%% subfunctions
function [the_code]=GetFunctionCode_2P(filename)
%this function opens the file and returns it as a string
a = fopen([ filename '.m']);
the_code_numbers=fread(a);
the_code = char(the_code_numbers);
end
function [s1,s2]=MatchStructs(s1,s2)
fields1=fieldnames(s1);
fields2=fieldnames(s2);
mpP=s2;
for n=1:length(fields1)
    if (strcmp(fields1{n},fields2)==0)
        for k=1:length(s2)
            s2(k).(fields1{n})=[];%=setfield(s2(k),fields1{n},[]);%put an empty field in
        end
    end
end

for n=1:length(fields2)
    if (strcmp(fields2{n},fields1)==0)
        for k=1:length(s1)
            s1(k).(fields2{n})=[]%s1(k)=setfield(s1(k),fields2{n},[]);%put an empty field in
        end
    end
end
end