function [vii,vij,partial2_vii,partial2_vij,...
    partial22_vii,partial22_vij] = BasisFUNv4(T,X,N,p1,p2)
%poly basis
    X = 2.*(X - min(X,[],1))./(max(X,[],1) - min(X,[],1)) - 1;
    T = 2.*(T - min(T))./(max(T) - min(T)) - 1;

    Tp1 = T.^(0:p1);

    d = size(X,2);

    Xp2 = zeros(N, d*p2);
    for j = 1:d
        Xp2(:,(j-1)*p2+1:j*p2) = (X(:,j)).^(1:p2);
    end
    Xp2 = [ones(N,1),Xp2];

    lX = size(Xp2,2);
    vii = zeros(N,(p1+1)*lX);
     for i = 1:p1+1
        vii(:,(i-1)*lX+1:i*lX) = Xp2.*Tp1(:,i);
    end

    if nargout > 1
        [indT,indX] = ndgrid((1:N)', (1:N)');
        indT = indT(:);
        indX = indX(:);
        Tij = T(indT);
        Xij = X(indX,:);
        
        Tijp1 = Tij.^(0:p1);
        %% vij
        Xijp2 = zeros(N^2, d*p2);
        for j = 1:d
            Xijp2(:, (j-1)*p2+1:j*p2) = (Xij(:,j)).^(1:p2);
        end
        Xijp2 = [ones(N^2,1),Xijp2];

        vij = zeros(N^2,(p1+1)*lX);
        for i = 1:p1+1
            vij(:,(i-1)*lX+1:i*lX) = Xijp2.*Tijp1(:,i);
        end
        
        %% partial2_vii
        d_Xp2 = zeros(N, d*p2);
        for j = 1:d
            d_Xp2(:,(j-1)*p2+1:j*p2) = (X(:,j)).^(0:p2-1).*repmat((1:p2),[N,1]);
        end
        d_Xp2 = [zeros(N,1),d_Xp2];

        partial2_vii = zeros(N,(p1+1)*lX);
        for i = 1:p1+1
            partial2_vii(:,(i-1)*lX+1:i*lX) = d_Xp2.*Tp1(:,i);
        end

        %% partial2_vij
        d_Xijp2 = zeros(N^2, d*p2);
        for j = 1:d
            d_Xijp2(:, (j-1)*p2+1:j*p2) = (Xij(:,j)).^(0:p2-1).*repmat((1:p2),[N^2,1]);
        end
        d_Xijp2 = [zeros(N^2,1),d_Xijp2];

        partial2_vij = zeros(N^2,(p1+1)*lX);
        for i = 1:p1+1
            partial2_vij(:,(i-1)*lX+1:i*lX) = d_Xijp2.*Tijp1(:,i);
        end

        %% partial22_vii
        dd_Xp2 = zeros(N, d*p2);
        for j = 1:d
            dd_Xp2(:,(j-1)*p2+1:j*p2) = (X(:,j)).^([0,0:p2-2]).*repmat((1:p2),[N,1]).*repmat((0:p2-1),[N,1]);
        end
        dd_Xp2 = [zeros(N,1),dd_Xp2];

        partial22_vii = zeros(N,(p1+1)*lX);
        for i = 1:p1+1
            partial22_vii(:,(i-1)*lX+1:i*lX) = dd_Xp2.*Tp1(:,i);
        end

         %% partial22_vij
        dd_Xijp2 = zeros(N^2, d*p2);
        for j = 1:d
            dd_Xijp2(:, (j-1)*p2+1:j*p2) = (Xij(:,j)).^([0,0:p2-2]).*repmat((1:p2),[N^2,1]).*repmat((0:p2-1),[N^2,1]);
        end
        dd_Xijp2 = [zeros(N^2,1),dd_Xijp2];

        partial2_vij = zeros(N^2,(p1+1)*lX);
        for i = 1:p1+1
            partial22_vij(:,(i-1)*lX+1:i*lX) = dd_Xijp2.*Tijp1(:,i);
        end
    end
    
end