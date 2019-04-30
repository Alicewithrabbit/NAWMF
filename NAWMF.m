function DenoisedImg = NAWMF(x1, ds, Ds)
% Copyright (c) 2019, Houwang Zhang & Chong Wu, All rights reserved.
% x1 is the image to be processed, ds should be set as 20, Ds should be set as 2.


    Y_max = max(max(x1));  Y_min = min(min(x1));
    [row,col]=size(x1);
    x1=double(x1);
    x2=x1;
    f=zeros(size(x1));
    wmax=39;
    w0=1;
    h=1;
    for i=1:row
        for j=1:col
            
        if (x1(i, j) > Y_min)&&(x1(i,j) < Y_max)
            continue
        end
        
            w = w0;
            while w <= wmax
                if w==w0
                    wsize=w;
                    c=x1(max(i-wsize,1):min(i+wsize,row),max(j-wsize,1):min(j+wsize,col));
                    W1=sort(c(:));
                    currmin = min(W1(:));
                    currmax = max(W1(:));
                else
                    W1=W2;
                    currmin=nextmin;
                    currmax=nextmax;
                end
                wsize=w+h;
                c=x1(max(i-wsize,1):min(i+wsize,row),max(j-wsize,1):min(j+wsize,col));
                W2 = c(:);
                nextmin = min(W2(:));
                nextmax = max(W2(:));
                if ~(currmin==nextmin && currmax==nextmax)
                    w=w+h;
                else
                    W11 = W1.*(W1>currmin).*(W1<currmax);
                    temp = W11(W11(:)~=0);
                    if ~isempty(temp)
                        currmean = mean(temp);
                        break;
                    else
                        w=w+h;
                    end
                end
            end
            if ~((currmin<x1(i,j)) && (x1(i,j)<currmax))
                x2(i,j) = currmean;
                f(i,j)=1;
            end
        end
    end
    dst=double((x2));
    
    noise = sum(sum(f))/(row*col);
    hs = 4.5595 + 6.0314*noise + 2.2186*(noise^2);

    I =dst;
    [m,n]=size(I);
    PaddedImg = padarray(I,[Ds+ds+1,Ds+ds+1],'symmetric','both');
    PaddedV = padarray(I,[Ds,Ds],'symmetric','both');
    average=zeros(m,n);
    wmax=average;
    sweight=average;
    h2=hs*hs;
    d=(2*ds+1)^2;
    for t1=-Ds:Ds
        for t2=-Ds:Ds
            if(t1==0&&t2==0)
                continue;
            end
            Sd=integralImgSqDiff(PaddedImg,Ds,t1,t2);
            SqDist2=Sd(2*ds+2:end-1,2*ds+2:end-1)+Sd(1:end-2*ds-2,1:end-2*ds-2)...
                -Sd(2*ds+2:end-1,1:end-2*ds-2)-Sd(1:end-2*ds-2,2*ds+2:end-1);
            SqDist2=SqDist2/d;
            w=exp(-SqDist2/h2);
            v = PaddedV(1+Ds+t1:end-Ds+t1,1+Ds+t2:end-Ds+t2);
            average=average+w.*v;
            wmax=max(wmax,w);
            sweight=sweight+w;
        end
    end
    average=average+0.*I;
    average=average./(0+sweight);
    DenoisedImg = average;

    for i = 1:m
        for j = 1:n
            if f(i, j) == 0
                DenoisedImg(i, j) = dst(i, j);
            end
        end
    end
    
end

function Sd = integralImgSqDiff(PaddedImg,Ds,t1,t2)
Dist2=(PaddedImg(1+Ds:end-Ds,1+Ds:end-Ds)-PaddedImg(1+Ds+t1:end-Ds+t1,1+Ds+t2:end-Ds+t2)).^2;
Sd = cumsum(Dist2,1);
Sd = cumsum(Sd,2);
end
