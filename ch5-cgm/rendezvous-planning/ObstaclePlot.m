gridmum = 200;
obc = length(Roc);
[cyl_x,cyl_y,cyl_z]=cylinder(1,gridmum-1);
for j = 1:obc
    for jj = 1:2
        Cyl_x = Roc(j)*Lnorm*cyl_x(jj,:);
        Cyl_y = Roc(j)*Lnorm*cyl_y(jj,:);
        Cyl_z = cyl_z(jj,:)*Lnorm*0.3;
        Cyl = Aoc{j,1}*[Cyl_x;Cyl_y];
        Cyl_Xtemp(jj,:) = Cyl(1,:);
        Cyl_Ytemp(jj,:) = Cyl(2,:);
        Cyl_Ztemp(jj,:) = Cyl_z;
    end
    Cyl_X{j,1} = Cyl_Xtemp+Xoc(j)*Lnorm;
    Cyl_Y{j,1} = Cyl_Ytemp+Yoc(j)*Lnorm;
    Cyl_Z{j,1} = Cyl_Ztemp;
end
for i =1:obc
    sh=surf(Cyl_X{i,1},Cyl_Y{i,1},Cyl_Z{i,1});%
    shading flat;
    set(sh,'linestyle','none');  %Òþ²ØÍø¸ñ
    alpha(0.65);
%     colormap(red);
%     colormap default;
    hold on;
end
%%
Cx=sin(0:0.1:2.1*pi);Cy=cos(0:0.1:2.1*pi);
Cxx = zeros(1,length(Cx));
Cyy = zeros(1,length(Cx));
for pli = 1:length(Roc)
    for i = 1:length(Cx)
        Cxy = Aoc{pli,1}*[Cx(i)*Roc(pli)*Lnorm;Cy(i)*Roc(pli)*Lnorm];
        Cxx(i) = Cxy(1);
        Cyy(i) = Cxy(2);
    end
    h = fill(Cxx+Xoc(pli)*Lnorm,Cyy+Yoc(pli)*Lnorm,'r');
    set(h,'EdgeColor','b','FaceColor','b');
    axis equal;
    hold on;
end
% title([method ' ' 'Method']);
text(300,100,['Iterations = ' num2str(s)],'FontName','Times New Roman','FontSize',15);