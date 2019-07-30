clear;
[I,map]=imread('/Users/fuweinan/Downloads/sky-and-birds.gif');
G=ind2gray(I,map);
figure;
imshow(G);
n=4;
ca=mat2cell(G,n*ones(1,size(G,1)/n),n*ones(1,size(G,2)/n));
DC=[];
Term=[];
AC=[];
RecTerm=[];
E=0;
I=0;
OriTerm=[];
Test=[];
for c=1:size(ca,1)
    for r=1:size(ca,2)
        temp = double(ca{c,r}); %// Change - cast to double for precision
        OriTerm=[OriTerm;temp];
        temp=temp-128;
        temp=dct2(temp);
        ca{c,r} = temp;
        DC=[DC,temp(1,1)]; 
        Term=[Term;ZigZagscan(temp)];
        AC=Term(:,2:n^2);
     end
end

LDC = min(DC);
HDC = max(DC);
partition = (LDC:(HDC-LDC)/16:HDC);
codebook = (LDC-(HDC-LDC)/32:(HDC-LDC)/16:HDC+(HDC-LDC)/32);
[index,quantsDC] = quantiz(DC,partition,codebook);

LAC = min(min(AC));
HAC = max(max(AC))+0.000001;
partitionAC1 = (LAC:(HAC-LAC)/8:HAC);
codebookAC1 = (LAC-(HAC-LAC)/16:(HAC-LAC)/8:HAC+(HAC-LAC)/16);

partitionAC2 = (LAC:(HAC-LAC)/4:HAC);
codebookAC2 = (LAC-(HAC-LAC)/8:(HAC-LAC)/4:HAC+(HAC-LAC)/8);
floor = floor((n^2-1)/10);
rem(1:n^2-2*floor-1) = 0;

for c=1:size(ca,1)
    for r=1:size(ca,2)
        if n == 2
            temp(1,1)=quantsDC(size(ca,2)*(c-1)+r);
            temp(1,2) = 0;
            temp(2,2) = 0;
            temp(2,1) = 0;
        else
            AC1 = Term(size(ca,2)*(c-1)+r,2:floor+1);
            AC2 = Term(size(ca,2)*(c-1)+r,floor+2:2*floor+1);
            [index2,quantsAC1] = quantiz(AC1,partitionAC1,codebookAC1);
            [index3,quantsAC2] = quantiz(AC2,partitionAC2,codebookAC2);
            array =[quantsDC(size(ca,2)*(c-1)+r),quantsAC1,quantsAC2,rem];
            array = double(array);
            temp = izigzag(array, n, n);                  
        end
        ca{c,r} = temp; 
     end
end

for c=1:size(ca,1)
    for r=1:size(ca,2)
        temp = ca{c,r}; %// Change - cast to double for precision
        temp=idct2(temp);
        temp=temp+128;
        for X=1:size(temp,1)
            for Y= 1:size(temp,2)
                if temp(X,Y)>255
                    temp(X,Y)= 255;
                elseif temp(X,Y)<0
                    temp(X,Y)= 1;
                end
            end
        end
        Test = [Test;ZigZagscan(temp)];
        ca{c,r} = temp;
        RecTerm=[RecTerm;temp];
     end
end

for x=1:size(OriTerm,1)
    for y=1:size(OriTerm,2)
        E = E+(OriTerm(x,y)-RecTerm(x,y))^2;
        I = I+OriTerm(x,y)^2;
    end
end
comrat = 8*(n^2)/(4+floor*5);
SNR = 10*(log10(I/E));
IG=cell2mat(ca);
imshow(IG/255);  
imshow(IG,[]);
imshow(uint8(IG));
Snr1 = [27.7103, 21.0493, 19.7494, 15.9944,12.7933,14.7067];
Snr2 = [24.5574, 20.6917, 16.8772, 11.6817, 9.8302,6.4777];
N = [2,4,8,16,32,64];
plot(N,Snr1,'-',N,Snr2,'--');
ylim([0 30]);
legend('Snr1','Snr2') 
xlabel('n');
ylabel('Snr');


