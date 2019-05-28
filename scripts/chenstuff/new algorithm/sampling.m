function y=sampling(x0,y0,x,control)
%sampling, the controal command can be 
%'zero': insert zeros if data is not valid (default)
%'ones': insert ones.....
%'near': insert the nearest number....
% and the length of x0 and y0 should be match

l=length(x);
y=zeros(1,l);

for i=1:l
    if x(i)>max(x0)
        id=x0==max(x0);
        switch lower(control)
            case 'ones'
                y(i)=1;
            case 'near'
                y(i)=y0(id);
            otherwise
                y(i)=0;
        end
    elseif x(i)<min(x0)
        id=x0==min(x0);
        switch lower(control)
            case 'ones'
                y(i)=1;
            case 'near'
                y(i)=y0(id);
            otherwise
                y(i)=0;
        end
    elseif x(i)==max(x0)
        id=x0==max(x0);
        y(i)=y0(id);
    elseif x(i)==min(x0)
        id=x0==min(x0);
        y(i)=y0(id);
    elseif x(i)==x0(end)
        y(i)=y(end);        
    elseif x0(1)<x0(end)
        id=sum(x(i)>=x0);
        y(i)=x(i)*(y0(id+1)-y0(id))/(x0(id+1)-x0(id))+(x0(id+1)*y0(id)-x0(id)*y0(id+1))/(x0(id+1)-x0(id));
    else
        id=sum(x(i)<=x0);
        y(i)=x(i)*(y0(id+1)-y0(id))/(x0(id+1)-x0(id))+(x0(id+1)*y0(id)-x0(id)*y0(id+1))/(x0(id+1)-x0(id));
    end
end


