 %% Script to Run M-SSA and Monte Carlo significance test using Allen & Robertson and/or Procrustres
% Gaston Manta, Eviatar Bach, Stefanie Talento, Marcelo Barreiro, Sabrina Speich and Michael Ghil
% 2023 

%% Import and cut the data, identify dependences,etc
% %%%[fList,pList] = matlab.codetools.requiredFilesAndProducts('mssa_nov_2023.m');
% 
% set(0,'defaultfigurecolor','w')
% 
% % SET UP THE DATA
% %era5_sst_slp_winds_world_1959_2022_winds_world_1959_2022.nc 
% ncdisp('era5_sst_slp_winds_world_1959_2022.nc')
% 
% %ue=ncread('era5_monthly_south_atlantic_winds_sst_slp_fluxes.nc','u10');
% sste=ncread('era5_sst_slp_winds_world_1959_2022.nc','sst');sste=sste-273.15;
% lone=ncread('era5_sst_slp_winds_world_1959_2022.nc','longitude');
% late=ncread('era5_sst_slp_winds_world_1959_2022.nc','latitude');
% timee=ncread('era5_sst_slp_winds_world_1959_2022.nc','time');timee=double(timee/24 +datenum('1900-01-01'));
% msl=ncread('era5_sst_slp_winds_world_1959_2022.nc','msl');msl=msl/100;
% 

% sste=squeeze(sste(:,:,1,:));
% msl=squeeze(msl(:,:,1,:));
% 
% 
% %ue=ue.*top';
% % msl=msl.*top;
% 
% top=sste(:,:,1)./sste(:,:,1);
% 
% 
% X=lone;
% X(X>180)=X(X>180)-360;
% [X,id_sort]=sort(X);
% sste=sste(id_sort,:,:);
% msl=msl(id_sort,:,:);
% lone=X;
% 
% 
% % Nmax=768;
% % j=1;
% % for i=1:12:Nmax
% %     sst_anual(:,:,j)=nanmean(sste(:,:,i:i+11),3);
% %     msl_anual(:,:,j)=nanmean(msl(:,:,i:i+11),3);
% %     j=j+1;
% % end
% % 
% 
% 
% 
% time=timee;
% timee=datevec(timee);
% 
% for i=1:((max(timee(:,1)))-(min(timee(:,1)))+1);
% 
% j=min(timee(:,1))+i-1;
% ind=find(timee==j);
% 
% sst_yearly(:,:,i)=nmean(sste(:,:,ind),3);
% %uee(:,:,i)=nmean(ue(:,:,ind),3);
% msl_yearly(:,:,i)=nmean(msl(:,:,ind),3);
% 
% end
% %yy=timee(1,1):timee(end,1);
% 
% %save regss_world sst_yearly msl_yearly yy top

load regss_world


%load regs_cera
% REGS ARE THE DATASETS PREPARED
%reg1=sst
%reg2=slp




lon=lone;
lat=late;

[lonn,latt]=meshgrid(lon,lat);

% % reg1=regg1;
% % reg2=regg2;


% choose domain
WL=-63;EL=20;SL=-50;NL=-10;


lat1=find(late==NL);lat2=find(late==SL);
lon1=find(lone==WL);lon2=find(lone==EL);

latew=late(lat1:lat2);
lonew=lone(lon1:lon2);
sstew=sst_yearly(lon1:lon2,lat1:lat2,:);
%ue=ue(lon1:lon2,lat1:lat2,:);
mslw=msl_yearly(lon1:lon2,lat1:lat2,:);


% % % figure;pcolor(lone,late,nmean(sste,3)');shading interp
% % % hold on
% % % contour(lone,late,nmean(ue,3)',[0 0],'k')
% % % 
% % % [llat,llon]=meshgrid(late,lone);llon=llon';llat=llat';
% % 

% get the data into spatially averaged ''spacestep''
spacestep=4;% now 1 degree

lat=latew(1:spacestep:end);
lon=lonew(1:spacestep:end);

reg1=sstew;reg2=mslw;
reg1=movmean(reg1,spacestep,1);reg1=movmean(reg1,spacestep,2);reg1=reg1(1:spacestep:end,1:spacestep:end,:);
reg2=movmean(reg2,spacestep,1);reg2=movmean(reg2,spacestep,2);reg2=reg2(1:spacestep:end,1:spacestep:end,:);

%builds a mask to remove data over land
top=reg1(:,:,1)./reg1(:,:,1);
reg2=reg2.*top;
time=yy;


lonf=lon;
latf=lat;




%load enso_yearly; yy=yearss;
%yearss=1950:2021;
%[C,ia,ib] = intersect(yy,yearss);
%enso_years=enso_yearly(ia,:);enso_years=enso_years-nmean(enso_years);

% mean
meanreg1=mean(reg1,3);meanreg2=mean(reg2,3);

%regs are now anomalies
reg1=reg1-meanreg1;reg2=reg2-meanreg2;
% reg3=reg3-meanreg3;

%the original values
regg1=reg1+meanreg1;
regg2=reg2+meanreg2;

% normalized anomalies
reg1=normalize(reg1,3);
reg2=normalize(reg2,3);

%% MSSA and Monte carlo parameters 

nn = 1000;        %number of times I want to repeat for the Monte Carlo
M = 14;           %window legth for M_SSA
N = size(reg1,3); %longitude of ts (#observations)
D = N;            %number of channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EOFs 
% Run PCA on original data so that the number of real channels gets
% reduced to L=N

%
[eof_maps1,pc1,expvar1]= eof(reg1,N);

if isnan(pc1(end,:))
pc1(end,:)=pc1(end-1,:);
eof_maps1(:,:,end)=eof_maps1(:,:,end-1);
end


if exist('/Users/user/Desktop/mssa/')
cd /Users/user/Desktop/mssa/
end

if not(isfolder('mssa_coupled'))
    mkdir('mssa_coupled')
end

if exist('/Users/user/Desktop/mssa/mssa_coupled/')
cd mssa_coupled
end

% %get pc from -1 to 1
% for k = 1:size(pc1,1)
% 
%    % Find the the maximum value in the time series of each principal
%   [maxval,ind] = max(abs(pc1(k,:)));
% 
%    % Divide the time series by its maximum value:
%    pc1(k,:) = pc1(k,:)/maxval;
% 
%    % Multiply the corresponding EOF map: eof_maps1(:,:,k) =
%    eof_maps1(:,:,k)*maxval;
% end


[eof_maps2,pc2,expvar2]= eof(reg2,N);

if isnan(pc2(end,:))
pc2(end,:)=pc2(end-1,:);
eof_maps2(:,:,end)=eof_maps2(:,:,end-1);
end

% I case you want to check the frequency
% figure
% for mm=1:10
%   subplot(10,1,mm);
% %   plot(PC1(mm,:),'k-');
%   ylabel(sprintf('PC %d',m));
% plotpsd(pc2(mm,:),1,'logx','lambda','r'); hold on
% 
% xlabel 'Periodicity (years)'
% 
% set(gca,'xtick',[0.1 .2 .3 .5 1:15 17 20])
% grid on; box on; axis tight
% 
% xlim([2   25])
% 
% % legend('SST')
%   title(sprintf('Spectral density of PC %d',mm));
% 
% set(gcf,'color','w');
% 
%   %ylim([-10 10]);
% end
% 
% 
% 
% %s = [-1 1 -1 1 -1 1]; % (sign multiplier to match Messie and Chavez 2011)
% 
% figure('pos',[100 100 500 700])
% for k = 1:6
%    subplot(3,2,k)
%  %  imagescn(lon,lat,eof_maps(:,:,k)*s(k));
%  %  imagescn(lon,lat,eof_maps(:,:,k)*s(k));
% pcolor(lonf,latf,eof_maps1(:,:,k)');shading interp
% cmocean('balance'); colorbar
%    title(['Mode ',num2str(k),' (',num2str(expvar1(k),'%0.1f'),'%)'])
%   % caxis([-2 2])
% end
% 
% 
% 
% 
% figure('pos',[100 100 500 700])
% for k = 1:6
%    subplot(3,2,k)
%  %  imagescn(lon,lat,eof_maps(:,:,k)*s(k));
%  %  imagescn(lon,lat,eof_maps(:,:,k)*s(k));
% pcolor(lonf,latf,eof_maps2(:,:,k)');shading interp
% cmocean('balance'); colorbar
%    title(['Mode ',num2str(k),' (',num2str(expvar2(k),'%0.1f'),'%)'])
%   % caxis([-2 2])
% end
% 
% figure('pos',[100 100 500 700])
% for k = 1:6
%    subplot(3,2,k)
%  %  imagescn(lon,lat,eof_maps(:,:,k)*s(k));
%  %  imagescn(lon,lat,eof_maps(:,:,k)*s(k));
% pcolor(lonf,latf,eof_maps5(:,:,k)');shading interp
% cmocean('balance'); colorbar
%    title(['Mode ',num2str(k),' (',num2str(expvar2(k),'%0.1f'),'%)'])
%   % caxis([-2 2])
% end


% From now on I will work over the PC matrix PC has size (number of pc,time)
X2 = pc2'; % X has N rows (#observations) & D==N columns (#channels) 
% % Make anomalies X at each channel
% for i=1:N
%     X2(:,i) = ( X2(:,i)-mean(X2(:,i)) );
% end

% From now on I will work over the PC matrix PC has size (number of pc,time)
X1 = pc1'; % X has N rows (#observations) & D==N columns (#channels) 
% % Make anomalies X at each channel
% for i=1:N
%     X1(:,i) = ( X1(:,i)-mean(X1(:,i)) );
% end

X=cat(2,X1,X2);
%X=cat(2,X2,X1);

% IN CASE YOU WANT TO RUN THE MSSA ONLY WITH SST
%X=X1;

% IN CASE YOU WANT TO RUN THE MSSA ONLY WITH SLP
%X=X2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run the M-SSA applied to X (it is a rank-defficient problem!):
    NN=size(X,2);

    %%Calculate covariance_surrogate matrix (trajectory approach)


    Y=zeros(N-M+1,M,NN);

    for m=1:M      % create time-delayed embedding of X
    for d=1:NN
        Y(:,m,d) = X((1:N-M+1)+m-1,d);
    end
    end

    Y=Y(:,:);
    Cemb=Y*Y' / (NN*M); %this is the reduced covariance matrix
    C=Cemb;

    
    %%%% Calculate eigenvalues LAMBDA and eigenvectors RHO
    [RHO_original,LAMBDA_original] = eig(C);
    [eigvalues,ind] = sort(diag(LAMBDA_original),'descend');
    LAMBDA_original = LAMBDA_original(ind,ind);
    RHO_original = RHO_original(:,ind);    % and eigenvectors: first M columns 
                                  % correspond to the first grid point
    %%% Calculate sigma to scale the eigenvectors later
    sigma_original = LAMBDA_original.^(1/2); %I need this scale matrix for later
    



explained_variance=eigvalues/(nsum(eigvalues))*100;





% 
% figure
% for mm=1:10
%   subplot(10,1,mm);
% %   plot(PC1(mm,:),'k-');
%  ylabel(sprintf('PC %d',m));
% plotpsd(RHO_original(:,mm),1,'logx','lambda','k'); hold on
% 
% xlabel 'Periodicity (years)'
% 
% set(gca,'xtick',[0.1 .2 .3 .5 1:15 17 20])
% grid on; box on; axis tight
% 
% xlim([2   25])
% 
% % legend('SST')
%   title(sprintf('Spectral density of PC %d',mm));
% 
% set(gcf,'color','w');
% 
%   %ylim([-10 10]);
% end
% 

% figure
% set(gcf,'name','Principal components PCs')
% clf;
% for m=1:8
%   subplot(8,1,m);
% % plot(pc2(m,:),'r-'); hold on
% % yyaxis right
% %  yyaxis right
%  plot(RHO_original(:,m),'k-'); hold on
% 
% %corpear=corr(pc2(m,:)',pc1(m,:)')
% 
%   ylabel(sprintf('PC %d',m));
% 
%   %title(sprintf(' corr %d',corpear));
% 
% 
%   %ylim([-10 10]);
% end
% 

%% Montecarlo 

% Surrogate data

%Initialize matrices, the last dimension is for each repetition of the
%Monte Carlo
LAMBDA_montecarlo = zeros(N-M+1,N-M+1,nn);
RHO_montecarlo = zeros(N-M+1,N-M+1,nn);
sigma_montecarlo= zeros(N-M+1,N-M+1,nn);
LAMBDA_montecarlo_AR= zeros(N-M+1,N-M+1,nn); %Allen & Robertson
LAMBDA_montecarlo_PRO= zeros(N-M+1,N-M+1,nn); %Procrustres

for i=1:nn;
    [x]=generate_surrogate_ts(X);
    %x has N rows (#observations) & D==N columns (#channels) 
    %At each grid point, x has mean 0 --> we work with anomalies 
    %x is fitted to have the same AR(1) & same std, at each point, as X 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Run the M-SSA applied to the surrogate data x:
    
    %%Calculate covariance_surrogate matrix (trajectory approach)
    Ym=zeros(N-M+1,M,N);
    for m=1:M      % create time-delayed embedding of X
    for d=1:N
        Ym(:,m,d) = x((1:N-M+1)+m-1,d);
    end
    end
    Ym=Ym(:,:);
    Cemb=Ym*Ym' / (N*M); %this is the reduced covariance matrix
    C=Cemb;
    
    %%%% Calculate eigenvalues LAMBDA and eigenvectors RHO
    [RHO,LAMBDA] = eig(C);
    [eigvalues,ind] = sort(diag(LAMBDA),'descend');
    LAMBDA = LAMBDA(ind,ind);
    RHO = RHO(:,ind);    % and eigenvectors: first M columns 
                         % correspond to the first grid point
    %%% Calculate sigma to scale the eigenvectors later
    sigma = LAMBDA.^(1/2); %I need this scale matrix for later
    
    LAMBDA_montecarlo(:,:,i) = LAMBDA;
    RHO_montecarlo(:,:,i) = RHO;
    sigma_montecarlo(:,:,i) = sigma;
    
    %%% Allen & Robertson
    LAMBDA_montecarlo_AR(:,:,i) = RHO_original' * C * RHO_original;
    
    %%% Procrustres target rotation
    %[distance,Z,transform] = procrustes(RHO_original*sigma_original,RHO_montecarlo(:,:,i)*sigma_montecarlo(:,:,i),'scaling',false);
    %T_P = transform.b.*transform.T;
    [result,transform,distance] = procrust(RHO_original*sigma_original,RHO_montecarlo(:,:,i)*sigma_montecarlo(:,:,i));
    T_P=transform;

    LAMBDA_montecarlo_PRO(:,:,i) = T_P' * LAMBDA_montecarlo(:,:,i) * T_P;

end
clear LAMBDA RHO




%calculate dominant frequency of each PC

%% Fast Fourier transform to find dominant frequency of all PCs
Fs=1; %Sampling frequency: 1 year
%find next power of 2
i=0;
while 2^i < N-M+1;
i = i+1;
end
nfft = 2^i;
for j=1:N-M+1    
%y = fft(PC(j,:),nfft); % Fast Fourier Transform
y = fft(RHO_original(:,j),nfft); % Fast Fourier Transform
y = abs(y.^2); % raw power spectrum density
y = y(1:1+nfft/2); % half-spectrum
[v,k] = max(y); % find maximum
f_scale = (0:nfft/2)* Fs/nfft; % frequency scale
%plot(f_scale, y),axis('tight'),grid('on'),title('Dominant Frequency')
freq_est(j) = f_scale(k); % dominant frequency estimate
period_est(j) = 1./freq_est(j);
end

%% Allen & Robertson

for j=1:N-M+1
    p1AR(j)=prctile(LAMBDA_montecarlo_AR(j,j,:),1);
    p99AR(j)=prctile(LAMBDA_montecarlo_AR(j,j,:),99); 
end

%% Procrustres
for j=1:N-M+1
    p1(j)=prctile(LAMBDA_montecarlo_PRO(j,j,:),1);
    p99(j)=prctile(LAMBDA_montecarlo_PRO(j,j,:),99); 
end

for j=1:N-M+1
    p5(j)=prctile(LAMBDA_montecarlo_PRO(j,j,:),5);
    p95(j)=prctile(LAMBDA_montecarlo_PRO(j,j,:),95);
end

% 
% figure 
% subplot(2,1,1)
% plot([period_est(1),period_est(1)],[p1AR(1),p99AR(1)],'-k')
% hold on; grid
% for j=2:N-M+1
% plot([period_est(j),period_est(j)],[p1AR(j),p99AR(j)],'-k');
% end
% plot(period_est(1:N-M+1),diag(LAMBDA_original),'*b') %ind has the order of the PCs
% title('Allen and Robertson')
% axis([1 N/4 min(diag(LAMBDA_original))*0.9 max(diag(LAMBDA_original))*1.1]) %only plot periods between sampling (1yr) % 1/4 of the total longitude f time series
% xlabel('Period of oscillation (yr)')
% ylabel('Power')
% 
% 
% 
% subplot(2,1,2)
% plot([period_est(1),period_est(1)],[p1(1),p99(1)],'-k')
% plot([period_est(1),period_est(1)],[p5(1),p95(1)],'-r')
% hold on; grid
% for j=2:N-M+1
% plot([period_est(j),period_est(j)],[p1(j),p99(j)],'-k');
% plot([period_est(j),period_est(j)],[p5(j),p95(j)],'-r');
% end
% plot(period_est(1:N-M+1),diag(LAMBDA_original),'*b') %ind has the order of the PCs
% title('Allen and Robertson')
% axis([1 N/4 min(diag(LAMBDA_original))*0.9 max(diag(LAMBDA_original))*1.1]) %only plot periods between sampling (1yr) % 1/4 of the total longitude f time series
% xlabel('Period of oscillation (yr)')
% ylabel('Power')
% 
% 
% 
% 
% figure 
% subplot(2,1,1)
% plot([freq_est(1),freq_est(1)],[p1AR(1),p99AR(1)],'-k')
% hold on; grid
% for j=2:N-M+1
% plot([freq_est(j),freq_est(j)],[p1AR(j),p99AR(j)],'-k');
% end
% 
% plot(freq_est(1:N-M+1),diag(LAMBDA_original),'*b') %ind has the order of the PCs
% 
% title('Allen and Robertson')
% 
% %axis([1 N/4 min(diag(LAMBDA_original))*0.9 max(diag(LAMBDA_original))*1.1]) %only plot periods between sampling (1yr) % 1/4 of the total longitude f time series
% xlabel('Frequency (cycles/year)')
% ylabel('Power')
% 
% xlim([-0.01 0.51])
% 
% ylim([0 51])
% 
% %aa=diag(LAMBDA_original)
% %explained_variance
% %(nsum(eigvalues))*100





% %% get it to 1 decimal
explained_variance=round(explained_variance,1);
period_est=round(period_est,1);


figure('Renderer', 'painters', 'Position', [50 50 250 450])

hw = subplot(5,1,1:2);
plot([freq_est(1),freq_est(1)],[p5(1),p95(1)],'-k')
hold on
%plot([freq_est(1),freq_est(1)],[p5(1),p99(1)],'-r')

hold on; grid

for j=2:N-M+1
plot([freq_est(j),freq_est(j)],[p5(j),p95(j)],'-k');
%plot([freq_est(j),freq_est(j)],[p5(j),p95(j)],'-r');

end
plot(freq_est(1:N-M+1),diag(LAMBDA_original),'db') %ind has the order of the PCs
%title('Procrustres')
%axis([1 N/4 min(diag(LAMBDA_original))*0.9 max(diag(LAMBDA_original))*1.1]) %only plot periods between sampling (1yr) % 1/4 of the total longitude f time series
xlabel('Frequency (cycles/year)')
ylabel('Power')

xlim([-0.01 0.51])

 ylim([0 250])

set(gcf,'color','w');

 pos = get(hw, 'Position') ;

 posnew = pos; posnew(2) = posnew(2) + 0.08; posnew(4) = posnew(4) - 0.06;

set(hw, 'Position', posnew)
posnew = pos;  set(hw, 'Position', posnew)

title('SST + SLP','fontweight','normal')

subplot(5,1,3);

plot(RHO_original(:,1),'k-');
title(['PC 1 | ',num2str(explained_variance(1)),'% exp. var. | ',num2str(period_est(1)), ' yr'],'fontweight','normal')

axis tight
subplot(5,1,4);

plot(RHO_original(:,2)+RHO_original(:,3),'k-');
title(['PCs 2-3 | ',num2str(explained_variance(2)+explained_variance(3)),'% exp. var. | ',num2str(period_est(3)), ' yr'],'fontweight','normal')

axis tight
subplot(5,1,5);

plot(RHO_original(:,4)+RHO_original(:,5),'k-');
title(['PCs 4-5 | ',num2str(explained_variance(4)+explained_variance(5)),'% exp. var. | ',num2str(period_est(5)), ' yr'],'fontweight','normal')
xlabel('Years')
% plot(RHO_original(:,5)+RHO_original(:,6),'k-');
% title(['PCs 5-6 | ',num2str(explained_variance(5)+explained_variance(6)),'% exp. var. | ',num2str(period_est(5)), ' yr'],'fontweight','normal')
% xlabel('Years')

axis tight
% 
% print(gcf,'-dpng','-r300','figure1_coupled')
% savefig('figure1_coupled')
% print(gcf,'-painters','-depsc2','-r600','figure1_coupled')
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%[rho,pval]=corr(enso_years(end-49:end,:),rho45)


%figure;plot(enso_years(1:50,end))

%% Calculate reconstructed components RC
% In order to determine the reconstructed components RC, we have to invert
% the projecting PC = Y*RHO; i.e. RC = Y*RHO*RHO'=PC*RHO'. Averaging along
% anti-diagonals gives the RCs for the original input X.

% Hay tantos RCs como grid points. Cada grid point tiene un conjunto de
% RCs.

var_eofs = eof_maps1;  % Change depending on variable
mask = ~isnan(var_eofs(:, :, 1)) & ~isinf(var_eofs(:, :, 1));

filtered = zeros(N, sum(sum(mask)));
for i=1:N
	slice = var_eofs(:, :, i);
	filtered(i, :) = slice(mask);
end

%var_eofs = permute(var_eofs, [3, 1, 2]);
%var_eofs = var_eofs(:, :);

RC_reduced=zeros(N-M+1,N,D);


Y1=Y(:,1:N*M);

RHO_full = Y1'*RHO_original;  % eigenvectors in D*M dimension

% Normalize eigenvectors
for k=1:N-M+1
	RHO_full(:, k) = RHO_full(:, k) ./ norm(RHO_full(:, k));
end



for k = 1:N-M+1
	ek = reshape(RHO_full(:, k), M, D);
	A = Y1*RHO_full;
        for n = 1:N;
        	if (1 <= n) & (n <= M - 1);
        		M_n = n;
        		L_n = 1;
			U_n = n;
        	elseif (M <= n) & (n <= N - M + 1);
        		M_n = M;
        		L_n = 1;
        		U_n = M;
        	elseif (N - M + 2 <= n) & (n <= N);
        		M_n = N - n + 1;
        		L_n = n - N + M;
        		U_n = M;
        	end
            	for d=1:D;
			for m=L_n:U_n;
				RC_reduced(k, n, d) = RC_reduced(k, n, d) + 1/M_n*A(n - m + 1, k)*ek(m, d);
			end
            	end
	end
end

%the first index is the mode index, the second is time, and the third and fourth are the grid location.
RC_full = NaN*zeros(N-M+1,N,size(var_eofs, 1),size(var_eofs, 2));

for k=1:N-M+1
	for i=1:N
		tmp = reshape(RC_reduced(k, :, :), N, N)*filtered;
		tmp2 = NaN*zeros(size(var_eofs, 1),size(var_eofs, 2));
		tmp2(mask) = tmp(i, :);
		RC_full(k, i, :, :) = tmp2;
	end
end



RC_fullsst=RC_full;
Rcc_fullsst=RC_full(:,:,:);

%% 2nd variable reconstruction


var_eofs = eof_maps2;  % Change depending on variable
mask = ~isnan(var_eofs(:, :, 1)) & ~isinf(var_eofs(:, :, 1));

filtered = zeros(N, sum(sum(mask)));
for i=1:N
	slice = var_eofs(:, :, i);
	filtered(i, :) = slice(mask);
end

%var_eofs = permute(var_eofs, [3, 1, 2]);
%var_eofs = var_eofs(:, :);


RC_reduced=zeros(N-M+1,N,D);

% EVIATAR CHECK THIS 

Y2=Y(:,N*M+1:N*M*2);

RHO_full = Y2'*RHO_original;  % eigenvectors in D*M dimension

% Normalize eigenvectors
for k=1:N-M+1
	RHO_full(:, k) = RHO_full(:, k) ./ norm(RHO_full(:, k));
end

for k = 1:N-M+1
	ek = reshape(RHO_full(:, k), M, D);
	A = Y2*RHO_full;
        for n = 1:N
        	if (1 <= n) & (n <= M - 1);
        		M_n = n;
        		L_n = 1;
			U_n = n;
        	elseif (M <= n) & (n <= N - M + 1)
        		M_n = M;
        		L_n = 1;
        		U_n = M;
        	elseif (N - M + 2 <= n) & (n <= N);
        		M_n = N - n + 1;
        		L_n = n - N + M;
        		U_n = M;
        	end
            	for d=1:D
			for m=L_n:U_n
				RC_reduced(k, n, d) = RC_reduced(k, n, d) + 1/M_n*A(n - m + 1, k)*ek(m, d);
			end
            	end
	end
end

%the first index is the mode index, the second is time, and the third and fourth are the grid location.
RC_full = NaN*zeros(N-M+1,N,size(var_eofs, 1),size(var_eofs, 2));

for k=1:N-M+1
	for i=1:N
		tmp = reshape(RC_reduced(k, :, :), N, N)*filtered;
		tmp2 = NaN*zeros(size(var_eofs, 1),size(var_eofs, 2));
		tmp2(mask) = tmp(i, :);
		RC_full(k, i, :, :) = tmp2;
	end
end


RC_fullslp=RC_full;
Rcc_fullslp=RC_full(:,:,:);



%%%%%%%%%%%%%%%   some correlations      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 
% rho45=RHO_original(:,4)+RHO_original(:,5);
% 
% rho23=RHO_original(:,2)+RHO_original(:,3);
% 
% [rho,pval]=corr(enso_years(end-50:end,:),rho45)
% 
% 
% %enso_smooth=movmean(enso_years,
% 
% ensoy=squeeze(nmean(enso_years,2))
% 
% ee=enso_years';ee=ee(:)
% 
% 
% ensoyy=movmean(enso_years,6,2);

% 
% j=1;
% for i=7:12:length(ee)-6
%     eew(j)=nmean(ee(i:i+11));
%     j=j+1;
% end




% [rho,pval]=corr(ensoy(end-50:end),rho45)
% 
% 
% [rho,pval]=corr(enso_years(end-50:end,:),rho45)
% 
% [rho,pval]=corr(enso_years(:,end),rho45)

%%   THE PLOT BEGINS                    %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultfigurecolor','w')

set(0,'defaultAxesFontSize',11)

setappdata(gcf, 'DefaultAxesXLabelFontSize', 11)

%% example plot of time series of reconstructed components

figure('Renderer', 'painters', 'Position', [150 150 900 1100])
set(gcf,'name','Reconstructed components RCs for grid points 1 and 100')
clf;
for m=1:M
  subplot(M,2,2*m-1);
  plot(Rcc_fullsst(m,:,70),'r-');
  ylabel(sprintf('RC %d',m));
  axis tight  
  
  subplot(M,2,2*m);
  %plot(RC(100,:,m),'r-');
  plot(Rcc_fullsst(m,:,70),'r-');
  ylabel(sprintf('RC %d',m));
  axis tight
end;


yy=1959:2022;
random_point=40;

%% FIG 3 EXAMPLE TIMESERIES
figure('Renderer', 'painters', 'Position', [150 150 600 300])

subplot(2,1,1)
hold on
plot(yy,squeeze(nsum(Rcc_fullsst(:,:,random_point))),'b');axis tight;grid on; box on
plot(yy,squeeze(nsum(Rcc_fullsst(1:5,:,random_point))),'r');axis tight;grid on; box on

legend('All RCs','RCs 1 to 5','box','off','location','northwest')
text(0.01,0.1,'A)','Units','normalized','FontSize',12)

title('SST')
subplot(2,1,2)
hold on

plot(yy,squeeze(nsum(Rcc_fullslp(:,:,random_point))),'b');axis tight;grid on; box on
plot(yy,squeeze(nsum(Rcc_fullslp(1:5,:,random_point))),'r');axis tight;grid on; box on
%legend('All RCs','RCs 1 to 5','box','off','location','northwest')
title('SLP')
text(0.01,0.1,'B)','Units','normalized','FontSize',12)


set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal')

print(gcf,'-painters','-depsc2','-r600','fig_example_timeseries');%close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  PHASES PLOT  %%%%%%%%%%%%%%%%%%%%

%% Trend

pp=1:4;

RC1sst=squeeze(Rcc_fullsst(1,:,:));

RC1slp=squeeze(Rcc_fullslp(1,:,:));


%RC23=squeeze(sum(Rcc_fullsst(2:3,:,:)));

RC23=squeeze(Rcc_fullslp(1,:,:));

%RC23=squeeze(sum(Rcc_fullslp(2:6,:,:)));


%SST
%ff=RC23(:,:);   %pasa a tener dim  (lt,ly*lx)

ff=RC23(:,:);   %pasa a tener dim  (lt,ly*lx)

mapf=ff(1,:);    % Save geographical dimension to use in mask.
ff(:,any(isnan(ff)))=[];%Remove NaNs
[lt,lxyf_sinNANs]=size(ff);   %dim (lt,lx*ly sin NaNs)

RC5y=ff;RC5y=cat(2,RC5y,RC5y(:,end));

[U,S,V]=svd(RC5y);

%Instantaneaous phase
phin=angle(U(:,1)+1i*U(:,2));

%Follows Moron et al (1998) to diplay mode. Uses Instantaneous phase
%to find times when the phase is in a certain band. Divides [-pi,pi]
%into 8 bands


RC23=squeeze(RC_fullsst(1,:,:,:));
RC23slp=squeeze(RC_fullslp(1,:,:,:));

%RC23=squeeze(sum(RC_fullsst(2:6,:,:,:),1));

stdreg1=squeeze(nstd(regg1,1,3));


for i=1:N
RC23_with_values(i,:,:)=squeeze(RC23(i,:,:)).*stdreg1;
end



%RC23slp=squeeze(Rcc_fullslp(1,:,:,:));


%RC23slp=squeeze(sum(RC_fullslp(2:6,:,:,:),1));

stdreg2=squeeze(nstd(regg2,1,3));

for i=1:N
RC23slp_with_values(i,:,:)=squeeze(RC23slp(i,:,:)).*stdreg2;
end


aa=jet(4)-0.2;aa(aa<0)=0.15;

%aa=[0.5 0 0;.5 .5 0; 1 .83 0;.3 0 .5];
yy=1959:2022;

anticyclone=squeeze(nmean(nmean(RC23slp_with_values,3),2));
anticyclonesst=squeeze(nmean(nmean(RC23_with_values,3),2));
color_yyaxis=[0     0.8     0];%green

figure('Renderer', 'painters', 'Position', [50 50 1150 500]); hold on

hw =subplot(3,4,1:4); 

%just to acomodate the size of this subplot
pos1 = get(hw, 'Position'); % gives the position of current sub-plot
new_pos1 = pos1 +[+0.03 +0.05 -0.045 -0.045];
set(hw, 'Position',new_pos1 ); % set new position of current sub - plot

hold on
title(['Trend. ', num2str(explained_variance(1)),'% of explained variance.'])



plot(1959:2022,anticyclone,'k')
hold on

ylim([-0.5 .5])

for p=1:4
    ii=find(((p-1)*pi/2-pi)<phin & phin<(p*pi/2-pi));
pp=p;
plot(yy(ii),anticyclone(ii),'+','markersize',7,'linewidth',1,'color',aa(p,:))
end
%yyaxis right

grid on;box on;

ylabel(' Mean SLPA (hPa)') 


yyaxis right
plot(1959:2022,anticyclonesst,'color',color_yyaxis)
ylabel('Mean SSTA (ºC)')
axis tight
ylim([-0.5 .5])

ax = gca;
ax.YAxis(2).Color = color_yyaxis;

text(0.01,0.9,'a)','Units','normalized','FontSize',20,'Fontweight','bold')


for p=1:4
    ii=find(((p-1)*pi/2-pi)<phin & phin<(p*pi/2-pi));
pp=p;
%just to assign the correct temporal order;
if p==1;p=3;elseif p==3;p=1;end

    subplot(3,4,p+8)
m_proj('mercator','lon',[-63 19], 'lat',[-50 -10]);
%m_proj('mercator','lon',[-68 25], 'lat',[-52 -5]);
hold on
    %cont_netcdf(xax,yax,squeeze(mean(RC5y(ii,:,:)))',0,(-1:0.2:1))

%m_pcolor(lonf,latf,squeeze(mean(RC23(ii,:,:)))');shading interp;% axis ij

m_pcolor(lonf,latf,squeeze(mean(RC23_with_values(ii,:,:)))');shading interp;% axis ij

caxis([-0.7 .7]);

cmocean('balance',15);

%cmocean('balance','pivot',0);colorbar

%cmocean('balance');colorbar

title(['SSTA Phase ' num2str(p)],'color',aa(pp,:))
m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none')

end
colorbar
fg=colorbar;fg.Label.String ='(ºC)';
fg.Position = fg.Position +   [.03, 0, 0, 0] ;

for p=1:4
    ii=find(((p-1)*pi/2-pi)<phin & phin<(p*pi/2-pi));
pp=p;
%just to assign the correct temporal order;
if p==1;p=3;elseif p==3;p=1;end

    subplot(3,4,p+4)
m_proj('mercator','lon',[-63 19], 'lat',[-50 -10]);
%m_proj('mercator','lon',[-68 25], 'lat',[-52 -5]);
hold on




%[FX,FY] = gradient(aa);
%grads=sqrt(FX.^2+FY.^2);
%cont_netcdf(xax,yax,squeeze(mean(RC5y(ii,:,:)))',0,(-1:0.2:1))

%m_pcolor(lonf,latf,squeeze(mean(RC23(ii,:,:)))');shading interp;% axis ij

 m_pcolor(lonf,latf,squeeze(mean(RC23slp_with_values(ii,:,:)))');shading interp;% axis ij
% 
 caxis([-0.8 .8]);cmocean('balance',17);
%cmocean('balance','pivot',0);colorbar

%cmocean('balance');colorbar
%m_pcolor(lonf,latf,grads');shading interp;colorbar
% caxis([-0.8 .8]);cmocean('balance',17);


aaw=squeeze(mean(RC23slp_with_values(ii,:,:)));aaw=aaw*2;

bb=aaw+nmean(regg2,3);

M = max(bb,[],'all');[ix,ij]=find(bb==M);

m_plot(lonf(ix),latf(ij),'+k')

m_contour(lonf,latf,bb',[1012 1016 1020],'k','showtext','on','Labelspacing',1500)


title(['SLPA Phase ' num2str(p)],'color',aa(pp,:))

m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none')


end

colorbar('location','EastOutside')
fg=colorbar;fg.Label.String ='(hPa)';
%just to avoid colorbar making the figure smaller
fg.Position = fg.Position +   [.03, 0, 0, 0];
% 


subplot(3,4,5)
text(0.01,0.9,'b)','Units','normalized','FontSize',20,'Fontweight','bold')


subplot(3,4,9)
text(0.01,0.9,'c)','Units','normalized','FontSize',20,'Fontweight','bold')



  print(gcf,'-dpng','-r600','Trend'); close all

%% 12.8 yr period

%RC23=squeeze(sum(Rcc_fullsst(2:3,:,:)));

RC23=squeeze(sum(Rcc_fullslp(2:3,:,:)));

%RC23=squeeze(sum(Rcc_fullslp(2:6,:,:)));


%SST
ff=RC23(:,:);   %pasa a tener dim  (lt,ly*lx)

mapf=ff(1,:);    % Save geographical dimension to use in mask.
ff(:,any(isnan(ff)))=[];%Remove NaNs
[lt,lxyf_sinNANs]=size(ff);   %dim (lt,lx*ly sin NaNs)

RC5y=ff;RC5y=cat(2,RC5y,RC5y(:,end));

[U,S,V]=svd(RC5y);

%Instantaneaous phase
phin=angle(U(:,1)+1i*U(:,2));

%Follows Moron et al (1998) to diplay mode. Uses Instantaneous phase
%to find times when the phase is in a certain band. Divides [-pi,pi]
%into 8 bands



RC23=squeeze(sum(RC_fullsst(2:3,:,:,:),1));


ff=RC23(:,:);   %pasa a tener dim  (lt,ly*lx)
% mapf=ff(1,:);    % Save geographical dimension to use in mask.
% ff(:,any(isnan(ff)))=[];%Remove NaNs
% [lt,lxyf_sinNANs]=size(ff);   %dim (lt,lx*ly sin NaNs)


%V33 = clusterdata(ff,'linkage','ward','savememory','on','maxclust',3);

clusters= clusterdata(ff','linkage','ward','savememory','on','maxclust',2);

claster=reshape(clusters,length(lon),length(lat));

hold on

ind=find(claster==2);d2=squeeze(nmean(RC23(:,ind),2));
ind=find(claster==1);d1=squeeze(nmean(RC23(:,ind),2));


% clasterr=movmedian(claster,3,1);clasterr=movmedian(clasterr,3,2);
% clasterr=movmedian(clasterr,3,1);clasterr=movmedian(clasterr,3,2);

%figure;pcolor(lon,lat,claster')
%RC23=squeeze(sum(RC_fullsst(2:6,:,:,:),1));

stdreg1=squeeze(nstd(regg1,1,3));


for i=1:N
RC23_with_values(i,:,:)=squeeze(RC23(i,:,:)).*stdreg1;
end


RC23slp=squeeze(sum(RC_fullslp(2:3,:,:,:),1));


%RC23slp=squeeze(sum(RC_fullslp(2:6,:,:,:),1));

stdreg2=squeeze(nstd(regg2,1,3));

for i=1:N
RC23slp_with_values(i,:,:)=squeeze(RC23slp(i,:,:)).*stdreg2;
end

yy=1959:2022;

anticyclone=squeeze(nmean(nmean(RC23slp_with_values,3),2));
anticyclonesst=squeeze(nmean(nmean(RC23_with_values,3),2));




% figure;plot(d1); hold on; plot(d2)

dipole_index=d1-d2;

% plot(anticyclone,'k')
% 
% hold on; plot(dipole_index,'color',color_yyaxis); hold on
% 
% legend('pole1','pole2','dipole_index','mean_slp')
% 

figure('Renderer', 'painters', 'Position', [50 50 1150 500]); hold on

hw =subplot(3,4,1:4); 

%just to acomodate the size of this subplot
pos1 = get(hw, 'Position'); % gives the position of current sub-plot
new_pos1 = pos1 +[+0.03 +0.05 -0.045 -0.045];
set(hw, 'Position',new_pos1 ); % set new position of current sub - plot

hold on

title(['12.8-year mode. ', num2str(sum(explained_variance(2:3))),'% of explained variance.'])

plot(1959:2022,anticyclone,'k')
axis tight
ylim([-0.5 .5])

for p=1:4
    ii=find(((p-1)*pi/2-pi)<phin & phin<(p*pi/2-pi));

plot(yy(ii),anticyclone(ii),'+','markersize',7,'linewidth',1,'color',aa(p,:))
end
%yyaxis right

grid on;box on;

ylabel(' Mean SLPA (hPa)') 

yyaxis right
plot(1959:2022,dipole_index,'color',color_yyaxis)
ylabel('SSTA dipole index (ºC)')
axis tight
%ylim([-0.5 .5])

ax = gca;
ax.YAxis(2).Color = color_yyaxis;


text(0.01,0.9,'a)','Units','normalized','FontSize',20,'Fontweight','bold')


for p=1:4
    ii=find(((p-1)*pi/2-pi)<phin & phin<(p*pi/2-pi));
pp=p;
%just to assign the correct temporal order;
if p==1;p=3;elseif p==3;p=1;end
    subplot(3,4,p+8)
m_proj('mercator','lon',[-63 19], 'lat',[-50 -10]);
%m_proj('mercator','lon',[-68 25], 'lat',[-52 -5]);
hold on
    %cont_netcdf(xax,yax,squeeze(mean(RC5y(ii,:,:)))',0,(-1:0.2:1))

%m_pcolor(lonf,latf,squeeze(mean(RC23(ii,:,:)))');shading interp;% axis ij

m_pcolor(lonf,latf,squeeze(mean(RC23_with_values(ii,:,:)))');shading interp;% axis ij

caxis([-0.25 .25]);

cmocean('balance',15);
% 
% %plot cluster region as dots
% cc=flip(claster');cc=cc(:);ind=find(cc==1);
% [lonn,latt]=meshgrid(lon,lat);latt=flip(latt);
% lonn=lonn(:);latt=latt(:);
% m_plot(lonn(ind),latt(ind),'.g')
%clasterr=movmedian(claster,3,1);clasterr=movmedian(clasterr,3,2);
m_contour(lon,lat,claster',1,'color',color_yyaxis,'linewidth',1)

%contour(lon,lat,claster',1,'color',color_yyaxis,'linewidth',1)

%m_contour(lon,lat,clasterr',[2 2],color_yyaxis,'linewidth',1)
%cmocean('balance','pivot',0);colorbar

%cmocean('balance');colorbar

 title(['SSTA Phase ' num2str(p)],'color',aa(pp,:))
m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none')

end
colorbar
fg=colorbar;fg.Label.String ='(ºC)';
%just to avoid colorbar making the figure smaller
fg.Position = fg.Position +   [.03, 0, 0, 0];

for p=1:4
    ii=find(((p-1)*pi/2-pi)<phin & phin<(p*pi/2-pi));
pp=p;
%just to assign the correct temporal order;
if p==1;p=3;elseif p==3;p=1;end

    subplot(3,4,p+4)
m_proj('mercator','lon',[-63 19], 'lat',[-50 -10]);
%m_proj('mercator','lon',[-68 25], 'lat',[-52 -5]);
hold on




%[FX,FY] = gradient(aa);
%grads=sqrt(FX.^2+FY.^2);
%cont_netcdf(xax,yax,squeeze(mean(RC5y(ii,:,:)))',0,(-1:0.2:1))

%m_pcolor(lonf,latf,squeeze(mean(RC23(ii,:,:)))');shading interp;% axis ij

 m_pcolor(lonf,latf,squeeze(mean(RC23slp_with_values(ii,:,:)))');shading interp;% axis ij
% 
 caxis([-0.8 .8]);cmocean('balance',17);
%cmocean('balance','pivot',0);colorbar

%cmocean('balance');colorbar
%m_pcolor(lonf,latf,grads');shading interp;colorbar
% caxis([-0.8 .8]);cmocean('balance',17);


aaw=squeeze(mean(RC23slp_with_values(ii,:,:)));aaw=aaw*2;

bb=aaw+nmean(regg2,3);

M = max(bb,[],'all');[ix,ij]=find(bb==M);

m_plot(lonf(ix),latf(ij),'+k')

m_contour(lonf,latf,bb',[1012 1016 1020],'k','showtext','on','Labelspacing',1500)


title(['SLPA Phase ' num2str(p)],'color',aa(pp,:));

m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none')


end
colorbar('location','EastOutside')
fg=colorbar;fg.Label.String ='(hPa)';
%just to avoid colorbar making the figure smaller
fg.Position = fg.Position +   [.03, 0, 0, 0];

subplot(3,4,5)
text(0.01,0.9,'b)','Units','normalized','FontSize',20,'Fontweight','bold')


subplot(3,4,9)
text(0.01,0.9,'c)','Units','normalized','FontSize',20,'Fontweight','bold')


  print(gcf,'-dpng','-r600','12yr_coupled'); close all

%% 5.3 yr period

%RC23=squeeze(sum(Rcc_fullsst(2:3,:,:)));

RC45=squeeze(sum(Rcc_fullslp(4:5,:,:)));

%SST
ff=RC45(:,:);   %pasa a tener dim  (lt,ly*lx)

mapf=ff(1,:);    % Save geographical dimension to use in mask.
ff(:,any(isnan(ff)))=[];%Remove NaNs
[lt,lxyf_sinNANs]=size(ff);   %dim (lt,lx*ly sin NaNs)

RC5y=ff;RC5y=cat(2,RC5y,RC5y(:,end));

[U,S,V]=svd(RC5y);

%Instantaneaous phase
phin45=angle(U(:,1)+1i*U(:,2));

%Follows Moron et al (1998) to diplay mode. Uses Instantaneous phase
%to find times when the phase is in a certain band. Divides [-pi,pi]
%into 8 bands


RC45=squeeze(sum(RC_fullsst(4:5,:,:,:),1));

stdreg1=squeeze(nstd(regg1,1,3));


for i=1:N
RC45_with_values(i,:,:)=squeeze(RC45(i,:,:)).*stdreg1;
end


RC45slp=squeeze(sum(RC_fullslp(4:5,:,:,:),1));

stdreg2=squeeze(nstd(regg2,1,3));

for i=1:N
RC45slp_with_values(i,:,:)=squeeze(RC45slp(i,:,:)).*stdreg2;
end


yy=1959:2022;

anticyclone=squeeze(nmean(nmean(RC45slp_with_values,3),2));
anticyclonesst=squeeze(nmean(nmean(RC45_with_values,3),2));


RC45=squeeze(sum(RC_fullsst(4:5,:,:,:),1));


ff=RC45(:,:);   %pasa a tener dim  (lt,ly*lx)
% mapf=ff(1,:);    % Save geographical dimension to use in mask.
% ff(:,any(isnan(ff)))=[];%Remove NaNs
% [lt,lxyf_sinNANs]=size(ff);   %dim (lt,lx*ly sin NaNs)



clusters= clusterdata(ff','linkage','ward','savememory','on','maxclust',2);


claster2=reshape(clusters,length(lon),length(lat));

figure;pcolor(lon,lat,claster2')

figure;pcolor(lon,lat,(claster2-claster)');colorbar



ind=find(claster2==1);dd1=squeeze(nmean(RC45(:,ind),2));

ind=find(claster2==2);dd2=squeeze(nmean(RC45(:,ind),2));
dipole_index2=dd1-dd2;


%figure;plot(dd1); hold on; plot(dd2)
%hold on; plot(dipole_index2,'r'); hold on
%plot(anticyclone,'k')
%legend('pole1','pole2','dipole_index','mean_slp')



figure('Renderer', 'painters', 'Position', [50 50 1150 500]); hold on

hw =subplot(3,4,1:4); 

%just to acomodate the size of this subplot
pos1 = get(hw, 'Position'); % gives the position of current sub-plot
new_pos1 = pos1 +[+0.03 +0.05 -0.045 -0.045];
set(hw, 'Position',new_pos1 ); % set new position of current sub - plot

hold on
title(['5.3-year mode. ', num2str(sum(explained_variance(4:5))),'% of explained variance.'])

plot(1959:2022,anticyclone,'k')
axis tight
ylim([-0.25 .25])

for p=1:4
    ii=find(((p-1)*pi/2-pi)<phin45 & phin45<(p*pi/2-pi));

plot(yy(ii),anticyclone(ii),'+','markersize',7,'linewidth',1,'color',aa(p,:))
end
%yyaxis right

grid on;box on;

 ylabel(' Mean SLPA (hPa)') 

yyaxis right;
plot(1959:2022,anticyclonesst,'color',color_yyaxis)
ylabel('SSTA dipole index (ºC)');axis tight
ylim([-0.2 .2])

ax = gca;
ax.YAxis(2).Color = color_yyaxis;

text(0.01,0.9,'a)','Units','normalized','FontSize',20,'Fontweight','bold')

% subplot(3,4,1:4)
% 
% hold on
% title('5.3 yr mode')
% 
% plot(1959:2022,anticyclone,'k')
% 
% axis tight
% ylim([-0.25 .25])
% 
% for p=1:4
%     ii=find(((p-1)*pi/2-pi)<phin & phin<(p*pi/2-pi));
% 
% plot(yy(ii),anticyclone(ii),'+','markersize',7,'linewidth',1,'color',aa(p,:))
% end
% %yyaxis right
% 
% grid on;box on;
% 
%  ylabel(' Mean SLPA (hPa)') 

for p=1:4
    ii=find(((p-1)*pi/2-pi)<phin45 & phin45<(p*pi/2-pi));
%just to assign the correct temporal order;
if p==1;p=3;elseif p==3;p=1;end


pp=p;
    subplot(3,4,p+8)
m_proj('mercator','lon',[-63 19], 'lat',[-50 -10]);
%m_proj('mercator','lon',[-68 25], 'lat',[-52 -5]);
hold on
    %cont_netcdf(xax,yax,squeeze(mean(RC5y(ii,:,:)))',0,(-1:0.2:1))

%m_pcolor(lonf,latf,squeeze(mean(RC23(ii,:,:)))');shading interp;% axis ij

m_pcolor(lonf,latf,squeeze(mean(RC45_with_values(ii,:,:)))');shading interp;% axis ij

caxis([-0.25 .25]);

cmocean('balance',15);

m_contour(lon,lat,claster2',1,'color',color_yyaxis,'linewidth',1)

%cmocean('balance','pivot',0);colorbar

%cmocean('balance');colorbar

title(['Phase ' num2str(p)],'color',aa(p,:))
m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none')

end
colorbar
fg=colorbar;fg.Label.String ='(ºC)';
%just to avoid colorbar making the figure smaller
fg.Position = fg.Position +   [.03, 0, 0, 0];


for p=1:4
    ii=find(((p-1)*pi/2-pi)<phin45 & phin45<(p*pi/2-pi));
%just to assign the correct temporal order;
if p==1;p=3;elseif p==3;p=1;end

    subplot(3,4,p+4)
m_proj('mercator','lon',[-63 19], 'lat',[-50 -10]);
%m_proj('mercator','lon',[-68 25], 'lat',[-52 -5]);
hold on




%[FX,FY] = gradient(aa);
%grads=sqrt(FX.^2+FY.^2);
%cont_netcdf(xax,yax,squeeze(mean(RC5y(ii,:,:)))',0,(-1:0.2:1))

%m_pcolor(lonf,latf,squeeze(mean(RC23(ii,:,:)))');shading interp;% axis ij

 m_pcolor(lonf,latf,squeeze(mean(RC45slp_with_values(ii,:,:)))');shading interp;% axis ij
% 
 caxis([-0.8 .8]);cmocean('balance',17);
%cmocean('balance','pivot',0);colorbar

%cmocean('balance');colorbar
%m_pcolor(lonf,latf,grads');shading interp;colorbar
% caxis([-0.8 .8]);cmocean('balance',17);


aaw=squeeze(mean(RC45slp_with_values(ii,:,:)));aaw=aaw*2;

bb=aaw+nmean(regg2,3);

M = max(bb,[],'all');[ix,ij]=find(bb==M);

m_plot(lonf(ix),latf(ij),'+k')

m_contour(lonf,latf,bb',[1012 1016 1020],'k','showtext','on','Labelspacing',1500)


title(['Phase ' num2str(p)],'color',aa(p,:));

m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none')


end
colorbar('location','EastOutside');
fg=colorbar;fg.Label.String ='(hPa)';
%just to avoid colorbar making the figure smaller
fg.Position = fg.Position +   [.03, 0, 0, 0];

subplot(3,4,5)
text(0.01,0.9,'b)','Units','normalized','FontSize',20,'Fontweight','bold')


subplot(3,4,9)
text(0.01,0.9,'c)','Units','normalized','FontSize',20,'Fontweight','bold')


  print(gcf,'-dpng','-r600','5yr_coupled')

%% 12.8 yr SST WORLD COMPOSITE 
%load /Users/gaston/Desktop/mssa_tesis/regss_world
sstaa=sst_yearly-nmean(sst_yearly,3);
mslaa=msl_yearly-nmean(msl_yearly,3);

%%%%%%%%%% % 12yr mode      %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sst
figure('Renderer', 'painters', 'Position', [50 50 1350 500])
for p=1:4
    ii=find(((p-1)*pi/2-pi)<phin & phin<(p*pi/2-pi));
pp=p;
%just to assign the correct temporal order;
if p==1;p=3;elseif p==3;p=1;end

psia2no=sstaa(:,:,ii);psia2no=permute(psia2no,[3 1 2]);

C=setxor(ii,1:length(sstaa(1,1,:)));

psia2nn=sstaa(:,:,C);psia2nn=permute(psia2nn,[3 1 2]);

fase2=C;fase=ii;
nin_neup=squeeze(mean(psia2no))-squeeze(mean(psia2nn));

spp=sqrt( ( var(psia2no(:,:))*(length(fase)-1)+var(psia2nn(:,:))*(length(fase2)-1) )/ (length(fase)+length(fase2)-2) );

ttp=nin_neup(:)'./(spp*sqrt(1/length(fase)+1/length(fase2)));
ttp=reshape(ttp,length(lone),length(late))';
dof=length(fase)+length(fase2)-2;
jjp=find(abs(ttp)<=tinv(0.975,dof));%0.975 for 5%, 0.95 for 10%
ttp(jjp)=NaN.*ones(size(jjp));

[lonee,latee]=meshgrid(lone,late);%latee=latee';
gap=10;
ttpp=ttp./ttp;ttpp2=ttp./ttp;

lonee=lonee(1:gap:end,1:gap:end);
latee=latee(1:gap:end,1:gap:end);

ttpp2=ttpp2(1:gap:end,1:gap:end);

ind=isnan(ttpp2);ind=find(ind==1);
lonee(ind)=[];latee(ind)=[];ttpp(ind)=[];




subplot(2,2,p)
%m_proj('mercator'); m_proj('mercator','lon',[-68 25], 'lat',[-52 -5]);
hold on
%cont_netcdf(xax,yax,squeeze(mean(RC5y(ii,:,:)))',0,(-1:0.2:1))
%m_pcolor(lonf,latf,squeeze(mean(RC23(ii,:,:)))');shading interp;% axis ij
%pcolor(lone,late,squeeze(mean(sstaa(:,:,ii),3))');shading interp;% axis ij
m_proj('robinson','clongitude',-60);
hold on
m_pcolor(lone,late,squeeze(mean(sstaa(:,:,ii),3))');shading interp;% axis ij
hold on
m_pcolor(lone-358,late,squeeze(mean(sstaa(:,:,ii),3))');shading interp;% axis ij
caxis([-0.5 .5]);
cmocean('balance',11);colorbar;
m_plot(lonee-358,latee,'.k');m_plot(lonee,latee,'.k');

m_grid('xtick',-180:90:180,'ytick',-80:40:80,'tickdir','out');

m_coast('patch',[0.81818 0.77647 0.70909]);
caxis([-.6 .6]);
 
 title(['12.8-year mode SSTA composite Phase ' num2str(p)],'color',aa(pp,:))

%title(['SSTA composite Phase ' num2str(pp)]);
cmocean('balance',21);colorbar;
%cmocean('balance','pivot',0);colorbar
% hold on m_contour(lonf,latf,squeeze(mean(slp23(ii,:,:)))',[1010 1015
% 1020],'k');%shading interp;% axis ij
% m_contour(lonf,latf,squeeze(mean(regg2i))',[1010 1015
% 1020],'--k');%shading interp;% axis ij
% 
%cmocean('balance');colorbar
fg=colorbar;fg.Label.String ='(ºC)';
 %letters_to_subplots(pp);
 %title(['Phase ' num2str(p)])
% m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
%m_grid('linestyle','none')


end
% 
% 
set(gcf,'color','w');

 print(gcf,'-dpng','-r600','world_sst_composite_12yr_mode_p05'); close all






% %%% slp % figure('Renderer', 'painters', 'Position', [50 50 1350 500])
% for p=1:4
%     ii=find(((p-1)*pi/4-pi)<phin & phin<(p*pi/4-pi));
% 
% subplot(2,2,p) %m_proj('mercator'); %m_proj('mercator','lon',[-68 25],
% 'lat',[-52 -5]); hold on
% %cont_netcdf(xax,yax,squeeze(mean(RC5y(ii,:,:)))',0,(-1:0.2:1))
% %m_pcolor(lonf,latf,squeeze(mean(RC23(ii,:,:)))');shading interp;% axis
% ij %pcolor(lone,late,squeeze(mean(sstaa(:,:,ii),3))');shading interp;%
% axis ij m_proj('robinson','clongitude',-60); hold on
% m_pcolor(lone,late,squeeze(mean(mslaa(:,:,ii),3))');shading interp;% axis
% ij hold on
% m_pcolor(lone-358,late,squeeze(mean(mslaa(:,:,ii),3))');shading interp;%
% axis ij % caxis([-0.5 .5]); %cmocean('balance',11);colorbar
% 
% m_grid('xtick',-180:90:180,'ytick',-80:40:80,'tickdir','out');
% %caxis([-.6 .6]);
% 
% cmocean('balance',21);colorbar
% 
% caxis([-1.5 1.5]) %cmocean('balance','pivot',0);colorbar % hold on %
% m_contour(lonf,latf,squeeze(mean(slp23(ii,:,:)))',[1010 1015
% 1020],'k');%shading interp;% axis ij %
% m_contour(lonf,latf,squeeze(mean(regg2i))',[1010 1015
% 1020],'--k');%shading interp;% axis ij % %cmocean('balance');colorbar
% %fg=colorbar;fg.Label.String ='(ºC)'; m_coast('color','k')
%  title(['Phase ' num2str(p)])
% % m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
% %m_grid('linestyle','none')
% 
% 
% end % % set(gcf,'color','w'); % % %
% print(gcf,'-dpng','-r600','world_msl_composite_12yr_mode'); close all % %
% %
% 









%% 5.3 yr SST WORLD COMPOSITE   %%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sst

figure('Renderer', 'painters', 'Position', [50 50 1350 500])
for p=1:4
    ii=find(((p-1)*pi/2-pi)<phin45 & phin45<(p*pi/2-pi));
pp=p;
%just to assign the correct temporal order;
if p==1;p=3;elseif p==3;p=1;end


psia2no=sstaa(:,:,ii);psia2no=permute(psia2no,[3 1 2]);

C=setxor(ii,1:length(sstaa(1,1,:)));

psia2nn=sstaa(:,:,C);psia2nn=permute(psia2nn,[3 1 2]);

fase2=C;fase=ii;
nin_neup=squeeze(mean(psia2no))-squeeze(mean(psia2nn));

spp=sqrt( ( var(psia2no(:,:))*(length(fase)-1)+var(psia2nn(:,:))*(length(fase2)-1) )/ (length(fase)+length(fase2)-2) );

ttp=nin_neup(:)'./(spp*sqrt(1/length(fase)+1/length(fase2)));
ttp=reshape(ttp,length(lone),length(late))';
dof=length(fase)+length(fase2)-2;
jjp=find(abs(ttp)<=tinv(0.95,dof));%0.975 for 5%, 0.95 for 10%
ttp(jjp)=NaN.*ones(size(jjp));

[lonee,latee]=meshgrid(lone,late);%latee=latee';
gap=10;
ttpp=ttp./ttp;ttpp2=ttp./ttp;

lonee=lonee(1:gap:end,1:gap:end);
latee=latee(1:gap:end,1:gap:end);

ttpp2=ttpp2(1:gap:end,1:gap:end);

ind=isnan(ttpp2);ind=find(ind==1);
lonee(ind)=[];latee(ind)=[];ttpp(ind)=[];

subplot(2,2,p)
%m_proj('mercator'); m_proj('mercator','lon',[-68 25], 'lat',[-52 -5]);
hold on
%cont_netcdf(xax,yax,squeeze(mean(RC5y(ii,:,:)))',0,(-1:0.2:1))
%m_pcolor(lonf,latf,squeeze(mean(RC23(ii,:,:)))');shading interp;% axis ij
%pcolor(lone,late,squeeze(mean(sstaa(:,:,ii),3))');shading interp;% axis ij
m_proj('robinson','clongitude',-60);
hold on
m_pcolor(lone,late,squeeze(mean(sstaa(:,:,ii),3))');shading interp;% axis ij
hold on
m_pcolor(lone-358,late,squeeze(mean(sstaa(:,:,ii),3))');shading interp;% axis ij
caxis([-0.5 .5]);
cmocean('balance',11);colorbar

m_plot(lonee-358,latee,'.k');m_plot(lonee,latee,'.k');

m_grid('xtick',-180:90:180,'ytick',-80:40:80,'tickdir','out');

m_coast('patch',[0.81818 0.77647 0.70909]);
caxis([-.6 .6]);

cmocean('balance',21);colorbar
%cmocean('balance','pivot',0);colorbar
% hold on m_contour(lonf,latf,squeeze(mean(slp23(ii,:,:)))',[1010 1015
% 1020],'k');%shading interp;% axis ij
% m_contour(lonf,latf,squeeze(mean(regg2i))',[1010 1015
% 1020],'--k');%shading interp;% axis ij
% 
%cmocean('balance');colorbar
fg=colorbar;fg.Label.String ='(ºC)';

 title(['5.3-year mode SSTA composite Phase ' num2str(p)],'color',aa(pp,:))
% m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
%m_grid('linestyle','none') letters_to_subplots(pp);

end
% 
% 
set(gcf,'color','w');



 print(gcf,'-dpng','-r600','world_sst_composite_5yr_mode_p010'); close all







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(lone)

for j=1:length(late)


variable=squeeze(sstaa(i,j,:));

[rhol(i,j),pvall(i,j)]=corr(variable,dipole_index);


end
end

ind=find(pvall>0.05);rhol(ind)=NaN;


figure;

hold on
m_proj('robinson','clongitude',-60);
hold on
m_pcolor(lone,late,rhol');shading interp;% axis ij
hold on
m_pcolor(lone-358,late,rhol');shading interp;% axis ij
caxis([-0.5 .5]);
cmocean('balance',11);colorbar


m_grid('xtick',-180:90:180,'ytick',-80:40:80,'tickdir','out');

m_coast('patch',[0.81818 0.77647 0.70909]);




























psia2no=sstaa(:,:,ii);psia2no=permute(psia2no,[3 1 2]);

C=setxor(ii,1:length(sstaa(1,1,:)));

psia2nn=sstaa(:,:,C);psia2nn=permute(psia2nn,[3 1 2]);

fase2=C;fase=ii;
nin_neup=squeeze(mean(psia2no))-squeeze(mean(psia2nn));

spp=sqrt( ( var(psia2no(:,:))*(length(fase)-1)+var(psia2nn(:,:))*(length(fase2)-1) )/ (length(fase)+length(fase2)-2) );

ttp=nin_neup(:)'./(spp*sqrt(1/length(fase)+1/length(fase2)));
ttp=reshape(ttp,length(lone),length(late))';
dof=length(fase)+length(fase2)-2;
jjp=find(abs(ttp)<=tinv(0.975,dof));
ttp(jjp)=NaN.*ones(size(jjp));

[lonee,latee]=meshgrid(lone,late);%latee=latee';
gap=8;
ttpp=ttp./ttp;ttpp2=ttp./ttp;

lonee=lonee(1:gap:end,1:gap:end);
latee=latee(1:gap:end,1:gap:end);

ttpp2=ttpp2(1:gap:end,1:gap:end);

ind=isnan(ttpp2);ind=find(ind==1);
lonee(ind)=[];latee(ind)=[];ttpp(ind)=[];






figure;

hold on
m_proj('robinson','clongitude',-60);
hold on
m_pcolor(lone,late,ttp);shading interp;% axis ij
hold on
m_pcolor(lone-358,late,ttp);shading interp;% axis ij

m_contour(lone,late,ttpp2,[1 1],'k');

m_plot(lonee,latee,'.k');
m_plot(lonee-358,latee,'.k');




% m_scatter(lonee(:),latee(:),ttpp(:),'.k');%shading interp
% m_scatter(lonee(:)-358,latee(:),ttpp(:),'.k');%shading interp








%caxis([-0.5 .5]);
cmocean('balance',11);colorbar

m_grid('xtick',-180:90:180,'ytick',-80:40:80,'tickdir','out');

m_coast('patch',[0.81818 0.77647 0.70909]);




disp('Check the figures in the folder')











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [x]=generate_surrogate_ts(X);
% 
% % Generate surrogate ts
% % The ts are AR(1), with temporal lag1 corr rho
% 
% %d = 5;          %number of ts (= number of spatial grid points)
% %n = 100;        %longitude of ts (=number of years of the observations)
% %rho is a vector with longitude d, tells the AR(1) coeff at each grid point
% 
% %x has n rows (#observations) & d columns (#spatial grid points) 
% 
% n = size(X,1);
% d = size(X,2);
% for i=1:d
%    aux=corrcoef(X(1:end-1,i),X(2:end,i));
%    rho(i)=aux(1,2);
% end
% 
% 
% %Generate the AR(1)
% for i=1:d;
% w=randn(n,1);
% x(1,i)=w(1);
%     for t=2:n
%         x(t,i) = rho(i).*x(t-1,i)+(1-rho(i)^2)^(1/2).*w(t);
%     end
% end
% 
% %make anomalies
% for i=1:d;
% x(:,i) = (x(:,i)-mean(x(:,i)));
% end;
% 
% %adjust std
% for i=1:d;
% x(:,i) = x(:,i).*std(X(:,i));
% end;
% end








function [x]=generate_surrogate_ts(X);

% Generate surrogate ts
% The ts are AR(1), with temporal lag1 corr rho

%d = 5;          %number of ts (= number of spatial grid points)
%n = 100;        %longitude of ts (=number of years of the observations)
%rho is a vector with longitude d, tells the AR(1) coeff at each grid point

%x has n rows (#observations) & d columns (#spatial grid points) 

n = size(X,1);
d = size(X,2);
for i=1:d
   aux=corrcoef(X(1:end-1,i),X(2:end,i));
   rho(i)=aux(1,2);
end


%Generate the AR(1)
for i=1:d;
w=randn(n,1);
x(1,i)=w(1);
    for t=2:n
        x(t,i) = rho(i).*x(t-1,i)+(1-rho(i)^2)^(1/2).*w(t);
    end
end

%make anomalies
for i=1:d;
x(:,i) = (x(:,i)-mean(x(:,i)));
end;

%adjust std
for i=1:d;
x(:,i) = x(:,i).*std(X(:,i));
end;
