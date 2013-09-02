function [pos,n]=evFindGroups(v,maxdiff,minsize)
%function [pos,n]=group(v,maxdiff,minsize)
%determine the number of groups of an increasing vector with as criterium that it's group 
%all the elements which are distant of striclty less than 'maxdiff' elements. It does not consider
%group that have less than 'minsize' elements.
%pos gives the positions related to the vector v of the begins and ends of each group
%example group([1 2 4 6 10 11 15 16 18],2,2) gives n=2 and pos = [1 7;      (modified by ben)
%                                                                 4 9]

minsize = minsize - 1; % Bug in GROUP: error of one element. Fix from ben, 6 march 2007

[r l]=size(v);
if r~=1
    v=v';
end

temp=v(2:end)-v(1:end-1);
if min(temp)<1|min([r l])>1 %not increasing vector | a single element
    pos=[];
    n=0;
else
    s=find(temp>maxdiff);
    pos(1)=1;
    t=2;
    for i=1:length(s)
        if (s(i)-pos(t-1))>=minsize
            pos(t)=s(i);
            t=t+1;
            pos(t)=s(i)+1;
            t=t+1;
        else
            pos=setdiff(pos,pos(t-1));
            pos(t-1)=s(i)+1;
        end
    end
    if (length(v)-pos(t-1))>=minsize
        pos(t)=length(v);
    else
        pos=setdiff(pos,pos(t-1));
        t=t-1;
    end   
    n=length(pos)/2;
end

pos = [pos(1:2:end); pos(2:2:end)]; %ben 
