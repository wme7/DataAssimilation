function eval = Ntildeyy(lambda_red)

% global Bred Fred k1 
global k1 TMP

eval = zeros(k1,k1);
for ii=1:k1
    for jj=1:ii
        % eval(ii,jj)=lambda_red'*Bred*(Fred(:,ii).*Fred(:,jj));      
        eval(ii,jj) = lambda_red'*TMP(:,ii,jj);
    end
end

eval = eval + eval' - diag(diag(eval)); % Hessian is symmetric




end



