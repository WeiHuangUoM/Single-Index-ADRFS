function [vii,vij] = BasisFUNv3(T,X,N,p1,p2)

    X = 2.*(X - min(X,[],1))./(max(X,[],1) - min(X,[],1)) - 1;
    T = 2.*(T - min(T))./(max(T) - min(T)) - 1;

    Tp1 = T.^(0:p1);

    d = size(X,2);

    
    X1 = (X(:,1)).^(0:p2(1));
    for j = 2:d
        lX = size(X1,2);
        Xj = zeros(N,(p2(j)+1)*lX);
        for i = 1:p2(j)+1
            Xj(:,(i-1)*lX+1:i*lX) = X1.*(X(:,j)).^(i-1);
        end
        X1 = Xj;
    end

    lX = size(X1,2);
    vii = zeros(N,(p1+1)*lX);
     for i = 1:p1+1
        vii(:,(i-1)*lX+1:i*lX) = X1.*Tp1(:,i);
    end

    if nargout > 1
        [indT,indX] = ndgrid((1:N)', (1:N)');
        indT = indT(:);
        indX = indX(:);
        Tij = T(indT);
        Xij = X(indX,:);
    
        Tijp1 = Tij.^(0:p1);

        X1 = (Xij(:,1)).^(0:p2(1));
        for j = 2:d
            lX = size(X1,2);
            Xj = zeros(N^2,(p2(j)+1)*lX);
            for i = 1:p2(j)+1
                Xj(:,(i-1)*lX+1:i*lX) = X1.*(Xij(:,j)).^(i-1);
            end
            X1 = Xj;
        end

        lX = size(X1,2);
        vij = zeros(N^2,(p1+1)*lX);
        for i = 1:p1+1
            vij(:,(i-1)*lX+1:i*lX) = X1.*Tijp1(:,i);
        end
    end
end