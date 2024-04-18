function problem2
N = 20;
Nt = N^2;
k = 1/Nt;
h = 1/(N-1);
grid = 0:1/(N-1):1;
U = zeros(N,N);
e = ones(N,1);
D2 = spdiags([e -2*e e],-1:1,N,N);
D2(1,2) = 2;
D2(N,N-1) = 2;
for nt = 1:Nt
    U = oneADIstep(U,(nt-1)*k);
    if nt == Nt*0.2
        Upt2 = U;
    elseif nt == Nt*0.4
        Upt4 = U;
    elseif nt == Nt*0.6
        Upt6 = U;
    elseif nt == Nt*0.8
        Upt8 = U;
    end
end

figure;
contourf(zeros(N,N));
figure;
contourf(Upt2);
figure;
contourf(Upt4);
figure;
contourf(Upt6);
figure;
contourf(Upt8);
figure;
contourf(U);


    function Unew = oneADIstep(Uold,t)
        Ustar = zeros(N,N);
        Unew = zeros(N,N);
        I = speye(N);
        b = zeros(N,1);

        for i = 1:N
            b(i) = b(i)+f(grid(i),grid(1),t);
        end
        b = b+(k/2)*(1/h^2)*(2*Uold(:,2)-2*Uold(:,1));
        Ustar(:,1) = (I-(k/2)*(1/h^2)*D2)\b;
        for j = 2:N-1
            for i = 1:N
                b(i) = b(i)+f(grid(i),grid(j),t);
            end
            b = b+(k/2)*(1/h^2)*(Uold(:,j+1)-2*Uold(:,j)+Uold(:,j-1));
            Ustar(:,j) = (I-(k/2)*(1/h^2)*D2)\b;
        end
        for i = 1:N
            b(i) = b(i)+f(grid(i),grid(N),t);
        end
        b = b+(k/2)*(1/h^2)*(2*Uold(:,N-1)-2*Uold(:,N));
        Ustar(:,N) = (I-(k/2)*(1/h^2)*D2)\b;
        
        b = zeros(N,1);
        for j = 1:N
            b(j) = b(j)+f(grid(1),grid(j),t);
        end
        b = b+(k/2)*(1/h^2)*(2*Ustar(2,:)-2*Ustar(1,:))';
        Unew(1,:) = ((I-(k/2)*(1/h^2)*D2)\b)';
        for i = 2:N-1
            for j = 1:N
                b(j) = b(j)+f(grid(i),grid(j),t);
            end
            b = b+(k/2)*(1/h^2)*(Ustar(i+1,:)-2*Ustar(i,:)+Ustar(i-1,:))';
            Unew(i,:) = ((I-(k/2)*(1/h^2)*D2)\b)';
        end
        for j = 1:N
            b(j) = b(j)+f(grid(N),grid(j),t);
        end
        b = b+(k/2)*(1/h^2)*(2*Ustar(N-1,:)-2*Ustar(N,:))';
        Unew(:,N) = ((I-(k/2)*(1/h^2)*D2)\b)';
    end

    function y = f(x,y,t)
        y = exp(-10*((x-0.6*cos(2*pi*t))^2+(y-sin(2*pi*t))^2));
    end
end