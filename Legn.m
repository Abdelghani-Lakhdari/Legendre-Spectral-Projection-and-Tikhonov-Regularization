function [y] = Legn(k,x)
if(k==0)
    y =1;
else if(k==1)
        y=x;
    else
        y =Legn(k-1,x)*x*(2*k-1)/(k)-((k-1)/k)*Legn(k-2,x);
    end
end

end