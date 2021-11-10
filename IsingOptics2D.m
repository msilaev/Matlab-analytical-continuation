clear all

G=0.01;

b=7;
%h=5.0;
Temp=0.01;

file=['Op'] ;  
load (file,'OOp', 'hh');

namebh=['DiffCondbEq7hEq5'];
Titlebh=['$\beta=7;\; h=5$'];

JJa=[];    
Oom=[];

StepE=G/2;
E=-20-0.1*pi:StepE:20+0.1*pi;

ER= E +j*G;
EA= E - j*G;

hhh=[];
ooom=[];
JJJa=[];

for indh=1:length(hh)
    
    D=OOp(indh);
    h=hh(indh)

for dOm=1:1:300

    OmRF=0.05*dOm;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
E1= E+OmRF ;
E2= E ;

Om1= j*E1; 
Om2=j*E2;
        
E1R= E1 +j*G;
E2R= E2 +j*G;

E1A= E1 - j*G;
E2A= E2 - j*G;

Om1R= j*E1R; 
Om2R=j*E2R;

Om1A= j*E1A; 
Om2A=j*E2A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JaR1R2 = Cond(Om1R,Om2R,D,h,b);
JaA1R2 = Cond(Om1A,Om2R,D,h,b);   
JaA1A2 = Cond(Om1A,Om2A,D,h,b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Ja= tanh(E1/(2*Temp)).*( JaR1R2 - JaA1R2 ) + tanh(E2/(2*Temp)).*( JaA1R2 - JaA1A2  );        
 
    %     JJa=[JJa StepE*sum(Ja)];
    
%Oom=[Oom OmRF];

hhh=[hhh h];
ooom=[ooom OmRF];
JJJa=[JJJa StepE*sum(Ja)];
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

end

tT = min(min(ooom)):0.01:max(max(ooom)); 
tL = min(min(hhh)):0.01:max(max(hhh)); 
[XI,YI] = meshgrid(tT,tL);

Zre = griddata(ooom,hhh, real(JJJa)./(16*ooom), XI,YI);

Zim = griddata(ooom,hhh, imag(JJJa)./(16*ooom), XI,YI);

Zre1=Zre;

Zre1=Zre;
Zim1=Zim;

% for ind=1:length(Zre(1,:))
%     
%     Zre1(:,ind)=Zre(:,ind)-Zre(:,length(Zim(:,1)));
%     
%         Zim1(:,ind)=Zim(:,ind)-0*Zim(:,1);
%     
% end 

for ind=1:length(Zre(:,1))
    
    Zre1(ind,:)=Zre(ind,:)-0*Zre(1,:);
    Zim1(ind,:)=Zim(ind,:)-0*Zim(1,:);
    
    end 
%%%%%%%%%%%%%%%%%%%%%%%%
figure
[h,h]=contourf(XI,YI,Zre,50)
set(h,'edgecolor','none');
colormap(jet)
colorbar

figure
[h,h]=contourf(XI,YI,Zim1,50)
set(h,'edgecolor','none');
colormap(jet)
colorbar

figure
[h,h]=contourf(XI,YI,Zim,50)

file=['Cond2D'] ;  
save (file,'ooom','hhh', 'JJJa');

% figure
%  plot(Oom, real(JJa)./(16*Oom),'linewidth',2)
%  hold on
%  plot(Oom, -imag(JJa)./(16*Oom),'linewidth',2)
%  
%  legend({'$Re (\sigma)$','$ Im (\sigma)$'},'interpreter','latex','Orientation','vertical','Location','east')
%  grid on
%  
%  set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',30)
% ylabel('$\sigma/\sigma_n$','interpreter','latex','FontSize',35)
% xlabel('$\omega/\Delta$','interpreter','latex','FontSize',35)
% 
% ylim([-0.25 1.5])
% title(Titlebh,'interpreter','latex')
% 
%   fname1=[namebh '.png']
%  print(gcf,fname1,'-dpng','-r300')
