function out = plotdata(alldata,group)

data = alldata(alldata(:,4)==group,1:3);
mouseid = unique(data(:,3))';

%figure
%for ii = min(mouseid):max(mouseid)
for ii = mouseid
    index = find(data(:,3)==ii);
    if(~isempty(index))
        plot(data(index,1),data(index,2),'ko-')
        hold on
    end
end
hold off