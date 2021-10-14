function newInd=DynamicTol2(Y1,num,changes,Tol)

for i=1:size(Y1,1)
    if num(i)>0
        [~,I]=find(changes(i,:)>0);
        Yd=max(Y1(i,1),Y1(i,I(1)+1))/min(Y1(i,1),Y1(i,I(1)+1));
        if length(I)==1
            Yd=[Yd max(Y1(i,end),Y1(i,I(1)+1))/min(Y1(i,end),Y1(i,I(1)+1))];
        else
            for j=1:length(I)-1   
                Yd=[Yd max(Y1(i,I(j)+1),Y1(i,I(j+1)+1))/min(Y1(i,I(j)+1),Y1(i,I(j+1)+1))];
            end
            Yd=[Yd max(Y1(i,end),Y1(i,I(end)+1))/min(Y1(i,end),Y1(i,I(end)+1))];
        end
        numTol=sum(Yd>Tol);
        if numTol>1
            newInd(i,1)=numTol; %The number of significant sign changes is added to the array that counts the total 3 of sign changes
            aa=sort(Yd,'descend');
            newInd(i,2)=min(aa(1:2));
        else
            newInd(i,1:2)=0;
        end
    else
        newInd(i,1:2)=0;
    end
end

end