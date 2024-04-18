function problem3
clear all; close all;
tol = 10^-5;
fd = @(p) max(-dcircle(p,0,0,1),dcircle(p,0,0,2));
%fh = @(p) min(min(0.03+0.3*dcircle(p,0,0,1),0.03+abs(0.3*dcircle(p,0,0,2))),2);
[pts,tri] = distmesh2d(fd,@huniform,0.1,[-2,-2;2,2],[0,2;2,0;0,-2;0,2;0,1;1,0;0,-1;0,1]);
dt = 0.01;

D0 = find(sqrt(pts(:,1).^2+pts(:,2).^2)<1+tol);
D1 = find(sqrt(pts(:,1).^2+pts(:,2).^2)>2-tol);
Npts = size(pts,1);
Ntri = size(tri,1);
FreeNodes = setdiff(1:Npts,union(D0,D1));
A = sparse(Npts,Npts);
b = sparse(Npts,1);
B = sparse(Npts,Npts);
u = sparse(Npts,1);

for j = 1:length(FreeNodes)
    k = FreeNodes(j);
    u(k) = norm(pts(k,:))+cos(atan(pts(k,2)/pts(k,1)));    
end

for j = 1:Ntri       
        A(tri(j,:),tri(j,:)) = A(tri(j,:),tri(j,:)) + stima3(pts(tri(j,:),:));
        B(tri(j,:),tri(j,:)) = B(tri(j,:),tri(j,:)) + det([1,1,1;pts(tri(j,:),:)'])*[2,1,1;1,2,1;1,1,2]/24;
        b(tri(j,:)) = b(tri(j,:))+det([1,1,1;pts(tri(j,:),:)'])*(1/6);%f equiv 1
end

figure;
trisurf(tri,pts(:,1),pts(:,2),full(u)','facecolor','interp')
hold on
axis ij
view(2)

for nt = 1:100
    RHS = (B-0.5*dt*A)*u+dt*b;
    LHS = 0.5*dt*A+B;
    u(FreeNodes) = LHS(FreeNodes,FreeNodes)\RHS(FreeNodes);
    if nt == 10
        upt1 = u;
    end
end

figure;
trisurf(tri,pts(:,1),pts(:,2),full(upt1)','facecolor','interp')
hold on
axis ij
view(2)

figure;
trisurf(tri,pts(:,1),pts(:,2),full(u)','facecolor','interp')
hold on
axis ij
view(2)

figure;
hold on
r = sqrt(pts(:,1).^2 + pts(:,2).^2);
[rsort,isort] = sort(r,'ascend');
usort = u(isort);
uexact = zeros(Npts,1);

for i = 1:length(rsort)
    r = rsort(i);
    uexact(i) = (1-r^2)/4+(3*log(r))/(4*log(2));
end

%disp(rsort)

a1 = plot(rsort,usort,'Linewidth',2); M1 = 'u';
a2 = plot(rsort,uexact,'Linewidth',2); M2 = 'uexact';
legend(M1,M2);

function M = stima3(vertices)
    d = size(vertices,2);
    G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
    M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);    
end

end
