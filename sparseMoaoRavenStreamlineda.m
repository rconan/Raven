%% ADAPTIVE OPTICS TOMOGRAPHY HOWTO
% Demonstrate how to build a tomographic adaptive optics system

%% Atmosphere 
Cn2 = [6.39 3.94 1.46 1.73 3.11 2.69 2.81];
fr0 = Cn2/sum(Cn2);
% atm = atmosphere(photometry.V,0.15,30,...
%     'altitude',[0,0.5,1,2,4,8,16]*1e3,...
%     'fractionnalR0',fr0,...
%     'windSpeed',[5.6,5.8,6.3,7.6,13.3,19.1,12.1],...
%     'windDirection',[0,5,15,30,60,90,180]*pi/180);
atm = atmosphere(photometry.V,0.15,30,...
    'altitude',[0]*1e3,...
    'fractionnalR0',[1],...
    'windSpeed',[5],...
    'windDirection',[0]);
atm.wavelength = photometry.R;

%% Telescope
tel = telescope(8,...%7.98
    'fieldOfViewInArcMin',2,...
    'samplingTime',1/500);

%% Calibration source
cSrc = source;

%% Wavefront sensor
nLenslet = 16;
wfs = shackHartmann(nLenslet,nLenslet*15);
wfs.tag = 'OPEN-LOOP WFS';
cSrc = cSrc.*tel*wfs;
wfs.lenslets.throughput = 0.75;
wfs.camera.exposureTime = 1/500;
wfs.camera.frameRate    = 500;
wfs.camera.quantumEfficiency = 0.8;
wfs.lenslets.minLightRatio = 0.8;
setValidLenslet(wfs)
+wfs;
wfs.referenceSlopes = wfs.slopes;
+wfs;
figure
imagesc(wfs.camera)
nPx = wfs.camera.resolution(1);

%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',25/100);
dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',wfs.validActuator);

% %% Zernike definition 
% maxRadialDegree = 12;
% zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
% zern = zernike(2:zernModeMax,tel.D,'resolution',nPx);

%% Sources
 % Calibration source
ngs = source('wavelength',photometry.R);
gs = source('asterism',{[3,60*constants.arcsec2radian,0]},...
     'wavelength',photometry.R,'magnitude',0,'tag','GUIDE STARS');
% gs = source('zenith',[60,60,60]*constants.arcsec2radian,'azimuth',[0,2*pi/3,0],...
%     'wavelength',photometry.R,'magnitude',0,'tag','GUIDE STARS');
scs = source('wavelength',photometry.R);
% scs = source('asterism',{[0,0],[30*cougarConstants.arcsec2radian,pi/4]},'wavelength',photometry.H);
% scs = source('asterism',{[0,0],...
%     [12, 30*cougarConstants.arcsec2radian, 0],...
%     [12, 60*cougarConstants.arcsec2radian, 0]},'wavelength',photometry.H,...
%     'tag','SCIENCE STARS');
% scs = source('asterism',{[0,0],...
%     [8, 45*cougarConstants.arcsec2radian, 0]},'wavelength',photometry.H,...
%     'tag','SCIENCE STARS');
nScs = length(scs);
nGs = length(gs);

%% Building the system
ngs=ngs.*tel*dm*wfs;
slopesAndFrameDisplay(wfs)
% % pause
% +ngs;
% slopesAndFrameDisplay(wfs)

%   %% DM/WFS calibration
% dm.coefsDefault = 0;
% stroke = 3;
% dm.coefs = eye(dm.nValidActuator)*stroke;
% +ngs;
% dm.coefs = 0;
% calibrationMatrix = wfs.slopes./stroke;
% figure(10)
% subplot(1,2,1)
% imagesc(calibrationMatrix)
% xlabel('DM actuators')
% ylabel('WFS slopes [px]')
% ylabel(colorbar,'slopes/actuator stroke')
% 
% %% Command matrix derivation
% [nS,nC] = size(calibrationMatrix);
% [U,S,V] = svd(calibrationMatrix);
% nThresholded = 4;
% eigenValues = diag(S);
% subplot(1,2,2)
% semilogy(eigenValues,'.')
% xlabel('Eigen modes')
% ylabel('Eigen values')
% iS = diag(1./eigenValues(1:end-nThresholded));
% iS(nC,nS) = 0;
% commandMatrix = V*iS*U';


% %% Zernike measurement
% zern.lex = false;
% % figure(10)
% % imagesc(zern.phase)
% zern.c = eye(zern.nMode);
% ngs=ngs.*zern*wfs;
% z = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree),tel.D)\wfs;
% Dz = z.c;
% 
%% With noise
withNoise = true;
ngs.magnitude = gs(1).magnitude;
wfs.camera.readOutNoise = 3;
wfs.camera.photonNoise = true;
wfs.framePixelThreshold = 3;
ngs=ngs.*tel*wfs;
slopesAndFrameDisplay(wfs)

%% noise convariance matrix
nMeas = 250;
slopes = zeros(wfs.nSlope,nMeas);
for kMeas=1:nMeas
    +wfs;
    slopes(:,kMeas) = wfs.slopes;
end
Cn = slopes*slopes'/nMeas;
iCn = sparse(1:wfs.nSlope,1:wfs.nSlope, 1./diag(Cn) );
figure(5)
% subplot(1,2,1)
imagesc(Cn)
axis equal tight
colorbar
% wfs.slopes = slopes;
% z = z\wfs;
% Czn = z.c*z.c'/nMeas;
% subplot(1,2,2)
% imagesc(Czn)
% axis equal tight
colorbar

%% Phase reconstruction
tel = tel+atm;
% %% wavefront reconstruction least square fit
% ngs=ngs.*tel;
% ps = ngs.meanRmPhase;
% ngs=ngs*wfs;
% z = z\wfs;
% zern.c = Dz\z.c(2:end);
% phaseLS = zern.phase;
% 
% %% wavefront reconstruction minimum variance
% Cz = phaseStats.zernikeCovariance(zern,atm,tel);
% M = Cz*Dz'/(Dz*Cz*Dz'+Czn);
% zern.c = M*z.c(2:end);
% phaseMV = zern.phase;
% figure(11)
% subplot(2,1,1)
% imagesc([ps,phaseLS,phaseMV])
% axis equal tight
% colorbar
% subplot(2,1,2)
% imagesc([ps-phaseLS,ps-phaseMV])
% axis equal tight
% colorbar

%% TOMOGRAPHY
% % Zernike section expansion
% modes = 2;%:zernike.nModeFromRadialOrder(30);
% zern = zernike(modes,tel.D,'resolution',atm.layer.nPixel);
% Gz = zernikeStats.covariance(zern,atm);
% L = chol(Gz,'lower');
% ps = zern.p*(L*randn(length(modes),1));
% ps = reshape(ps,atm.layer.nPixel,atm.layer.nPixel);
% atm.layer.phase = ps;
% ps = tel.pupil.*( ps - mean(ps(tel.pupilLogical)) );
%%
figure(3)
subplot(1,atm.nLayer+1,[1,atm.nLayer])
imagesc(tel,[gs,scs])
subplot(1,atm.nLayer+1,1+atm.nLayer)
polar(scs)
hold on
polar(gs,'ro')
hold off
% figure(4)
% imagesc(ps)
% axis square
% colorbar

%% The matrices
[Gamma,gridMask] = sparseGradientMatrix(wfs);
nGeom     = length(gridMask);
layerGrid = cell(1,atm.nLayer);
for kLayer=1:atm.nLayer
    if atm.layer(kLayer).altitude==0
        layerGrid{kLayer} = gridMask;
    else
        layerGrid{kLayer} = utilities.piston( round(nGeom.*atm.layer(kLayer).D./tel.D), 'type', 'logical');
    end
end
[H,inG, outG , interpPhase] = sparseInterpMatrix(layerGrid,gridMask,tel,gs,atm);
Hscs = sparseInterpMatrix(layerGrid,gridMask,tel,scs,atm);
L2 = phaseStats.sparseInverseCovarianceMatrix(layerGrid,atm);

G   = repmat( {Gamma}, 1, length(gs));
G   = blkdiag( G{:} ) * cell2mat(H);
L2  = blkdiag( L2{:} );
iGn = repmat( {iCn}, 1, length(gs));
iGn = blkdiag( iGn{:} );
regul = 1e5; % Tikonov regularisation

gs = gs.*tel*wfs;

d = tel.D/nLenslet;
q = 2*pi/d/(wfs.lenslets.nyquistSampling*2)/2;

A = G'*iGn*G+L2+speye(size(L2)).*regul;
b = G'*iGn*(wfs.slopes(:)*q);
% M   = (G'*iGn*G+L2+speye(size(L2)).*regul)\(G'*iGn);
buffer = A\b;
% buffer = (G'*iGn*G+L2)\(G'*iGn*wfs.slopes(:));

nPxLayer = cellfun( @(x) sum(x(:)), layerGrid);
buffer   = mat2cell( buffer , nPxLayer , 1 );
psEst    = cell(1,atm.nLayer);
for kLayer=1:atm.nLayer
    psEst{kLayer} = zeros(size(layerGrid{kLayer}));
    psEst{kLayer}(layerGrid{kLayer}) = buffer{kLayer};
end

for kLayer = 1:atm.nLayer
    n = length(atm.layer(kLayer).phase);
    [x,y] = meshgrid( (1-n:2:n-1)/(n-1)*length(inG{1,kLayer}{1})/2 );
    F = TriScatteredInterp(x(:),y(:),atm.layer(kLayer).phase(:));
    fullScreen = F(inG{1,kLayer}{1},inG{1,kLayer}{2});
    cutScreen = F(outG{1,kLayer}{1},outG{1,kLayer}{2});
    figure
    subplot(1,3,1)
    p = utilities.piston(length(fullScreen),'type','logical');
    imagesc(p.*(fullScreen - mean(fullScreen(p))))
    axis xy square
    clim = get(gca,'clim');
    subplot(1,3,2)
    p = utilities.piston(length(cutScreen),'type','logical');
    imagesc(p.*(cutScreen - mean(cutScreen(p))))
    axis xy square
    set(gca,'clim',clim)
    subplot(1,3,3)
    imagesc(layerGrid{kLayer}.*(psEst{kLayer} - mean(psEst{kLayer}(layerGrid{kLayer}))))
    axis xy square
    set(gca,'clim',clim)
end

%%
% %% Turbulence covariance matrix
% % tic
% % S = phaseStats.zernikeAngularCovariance(zern,atm,gs);
% % toc
% % S = cell2mat(S);
% load('S12')
% %% Data/Target covariance
% % C = phaseStats.zernikeAngularCovariance(zern,atm,gs,scs);
% load('C12_45')
% %% tomographic matrices
% CznAst = blkdiag( Czn , Czn , Czn );
% DzAst = blkdiag( Dz , Dz , Dz );
% 
% % %% wavefront reconstruction least square fit
% % gs=gs.*tel;
% % ps = [gs.meanRmPhase];
% % gs=gs*wfs
% % z = z\wfs;
% % zern.c = Dz\z.c(2:end,:);
% % phaseLS = reshape(zern.phase,nPx,nPx*length(gs));
% % 
% % %% wavefront reconstruction minimum variance
% % M = S/(S+CznAst);
% % M = S*DzAst'/(DzAst*S*DzAst'+CznAst);
% % z = z - 1; % piston removed
% % zern.c = reshape(M*z.c(:),z.nMode,[]);
% % phaseMV = reshape(zern.phase,nPx,nPx*length(gs));
% % figure(12)
% % imagesc([ps;phaseLS;phaseMV])
% % axis equal tight xy
% 
% %% Target command matrix
% M = cell(nScs,1);
% denom = DzAst*S*DzAst'+CznAst;
% for kScs = 1:nScs
%     M{kScs,1} = cell2mat(C(:,kScs))'*DzAst'/denom;
% end
% M = cell2mat(M);
% 
% %% Command matrix
% lambdaRatio = 1;%gs(1).wavelength/scs(1).wavelength;
% gs=gs.*tel*wfs;
% z = zernike(1:zernModeMax)\wfs;
% zern.c = reshape(M*(lambdaRatio*z.c(:)),z.nMode,nScs);
% ngsWavelength  = ngs.wavelength;
% % ngs.wavelength = scs(1).wavelength;
% ngs = ngs.*zern;
% % ngs.wavelength = ngsWavelength;
scs = scs.*tel;
% turbPhase = [scs.meanRmPhase];
nIt =50;
turbPhaseStd = zeros(nIt,nScs);
turbPhaseStd(1,:) = scs.var;
% figure
% imagesc([scs.meanRmPhase;reshape(ngs.phase,nPx,[])*ngs.wavelength/scs(1).wavelength])
% axis equal tight xy
% colorbar
% 
% %% DM
% zern.lex = true;
% zernC = zern.c;
% dm.coefs = eye(dm.nValidActuator);
% tel = tel - atm;
% ngs = ngs.*tel*dm; % influence functions
% % buf = utilities.toggleFrame(ngs.phase);
% % buf = buf(ngs.mask,:);
% % % zern2dm = dm.modes.modes(tel.pupilLogical,:)\zern.p(tel.pupilLogical,:)/2;%\zPoly/2;
% % zern2dm = -buf\zern.p(tel.pupilLogical,:);
% zern2dm = -(ngs\zern);
% %zern2dm = -zern2dm*ngs.wavelength/scs(1).wavelength;
% dm.coefs = zern2dm*zernC;
% scs = scs.*tel*dm;
% correction = [scs.meanRmPhase];
% tel = tel + atm;
% scs = scs.*tel*dm;
% residue = [scs.meanRmPhase];
% % turbRes = zeros(nPx,nPx*nScs,nIt);
% turbRes = turbPhase;
turbResStd = turbPhaseStd;
% figure
% % plot([turbPhaseStd(1:k,:),turbResStd(1:k,:)],'.');
% % set(h(1),'YDataSource',turbPhaseStd(:,1))
% h = imagesc([turbPhase;residue;-correction]);
% axis equal tight xy
% colorbar
% 
%% Open-Loop
% DM with down--sampled influence functions
bifGeom = influenceFunction('monotonic',25/100);
dmGeom = deformableMirror(nActuator,...
    'modes',bifGeom,...
    'resolution',nGeom,...
    'validActuator',wfs.validActuator);

k = 1;
% scs(1).saveImage = true;
% scs(2).saveImage = true;
log = logBook.checkIn;
log.verbose=false;
normTel = sum(tel.pupil(:));
nOtf = 2*nPx;
otfTel = fftshift(ifft2(abs(fft2(tel.pupil,nOtf,nOtf)).^2))/normTel;
meanOtfPd = 0;
warning off MATLAB:rankDeficientMatrix
fprintf(' --> Open-loop started at %s, %d iterations\n',datestr(now),nIt);
tic
while k<nIt
    
    % propagation of the guide stars to the WFSs
    gs=gs.*tel*wfs;
    % Atmospheric phasescreens estimation and science phase derivation
    b = G'*iGn*(wfs.slopes(:)*q);
    dmPhase = blkdiag(Hscs{:})*(A\b);
    
        dm.coefs = (dmGeom.modes.modes(gridMask(:),:)\reshape(dmPhase,sum(gridMask(:)),nScs))/ngs.waveNumber/2;


    % slopes projection onto the Zernikes, piston removed
%     z = z\wfs;
%     z.c = wfs.zernCoefs;
    % DMs command vectors
%     dm.coefs = lambdaRatio*zern2dm*reshape(M*wfs.zernCoefs(:),z.nMode,nScs);
    % Atmosphere update
    +tel;
    % propagation of science star to the telescope through the atmosphere
    scs = scs.*tel;
    k = k + 1 ;
    turbPhaseStd(k,:) = scs.var;
    % propagation of science star resumes to the DMs
    scs = scs*dm;
    turbRes = cat(3,scs.meanRmPhase);
    turbResStd(k,:) = scs.var;
    if k>1
        residualWave = bsxfun(@times,tel.pupil,exp(1i.*turbRes));
        meanOtfPd = meanOtfPd + abs(fft( fft( residualWave, nOtf, 1) , nOtf, 2)).^2;
    end

%     set(h,'Ydata',[turbPhaseStd,turbResStd])
%     set(h,'Cdata',[turbPhase;turbRes(:,:,k);reshape(-dm.phase,nPx,[])])
%     drawnow
    
end
meanOtfPd = fftshift( fftshift( ifft (ifft( meanOtfPd, nOtf, 1), nOtf, 2)/normTel, 1), 2)/(nIt-1);
meanOtfPd = reshape(meanOtfPd,nOtf,nOtf*nScs);
meanOtfPd = mat2cell(meanOtfPd,nOtf,nOtf*ones(1,nScs));
toc

%%
u = [1 nIt];
atmWavelength = atm.wavelength;
atm.wavelength = scs(1).wavelength;
rmsMicron = @(x) 1e6*sqrt(x).*atm.wavelength/2/pi;
figure
plot(rmsMicron(turbPhaseStd),'.')
hold all
plot(rmsMicron(turbResStd),'.')
hold off
line(u,rmsMicron(phaseStats.zernikeResidualVariance(1,atm,tel))*ones(2,1),'color','k','LineWidth',2)
line(u,ones(2,1)*mean(rmsMicron(turbPhaseStd)),'lineStyle','--')
% line(u,phaseStats.zernikeResidualVariance(zern.nMode,atm,tel)*ones(2,1),'color','r')
grid
xlabel('Iterations')
ylabel('Rms [\mum]')
atm.wavelength = atmWavelength;
% 
% %%
% % % Optical transfer function
% % normTel = sum(tel.pupil(:));
% % nOtf = 2*nPx;
% % otfTel = fftshift(ifft2(abs(fft2(tel.pupil,nOtf,nOtf)).^2))/normTel;
% % bigRam = false;
% % if bigRam
% %     turbRes = reshape(turbRes,nPx,nPx*nScs*nIt);
% %     turbRes = mat2cell(turbRes,nPx,nPx*ones(1,nScs*nIt));
% %     residualWave = cellfun(@(x)tel.pupil.*exp(1i*x),turbRes,'UniformOutput',false);
% %     tic
% %     otfPd = cellfun(@(x)fftshift(ifft2(abs(fft2(x,nOtf,nOtf)).^2))/normgTel,residualWave,'UniformOutput',false);
% %     toc
% %     clear residualWave
% %     for kScs = 1:nScs
% %         otfPd{kScs} = mean(reshape(cell2mat(otfPd((nScs+kScs):nScs:end)),nOtf,nOtf,nIt-1),3);
% %     end
% %     meanOtfPd = otfPd(1:nScs);
% %     clear otfPd
% % else
% %     meanOtfPd = cell(1,nScs);
% %     h = waitbar(0,'Computing OTFs ...');
% %     tic
% %     for kScs = 1:nScs
% %         otfPd = bsxfun(@times,tel.pupil,exp(1i.*turbRes(:,:,(nScs+kScs):nScs:end)));
% %         otfPd = mean( abs(fft( fft( otfPd, nOtf, 1) , nOtf, 2)).^2 , 3);
% %         meanOtfPd{kScs} = fftshift(ifft2( otfPd, nOtf, nOtf)/normTel);
% %         waitbar(kScs/nScs)
% %     end
% %     toc
% %     close(h)
% %     clear residualWave otfPd
% % end
% % Strehl ratio
% u = linspace(-tel.D,tel.D,nOtf);
% strehlRatioFun = @(x)real(trapz(u,trapz(u,x)))/tel.area;
% strehlRatio = cellfun(strehlRatioFun,meanOtfPd);
% % entraped energy
% [x,y] = meshgrid(u);
% ircsSlitWidthInArcsec = 0.14';
% a = (ircsSlitWidthInArcsec/(photometry.H/tel.D*constants.radian2arcsec))/tel.D; % diameter
% eeFilter = a^2*(sin(pi.*x.*a)./(pi.*x.*a)).*...
%     (sin(pi.*y.*a)./(pi.*y.*a));
% eNrgFun = @(x) real(trapz(u,trapz(u,x.*eeFilter)));
% fprintf(' > Diffraction limited nrg (fine/coarse) : %4.1f%%/%4.1f%%\n',...
%     entrappedEnergy(tel-atm,a/2,'square','otf')*100,eNrgFun(otfTel)*100)
% tel = tel + atm;
% eNrg = cellfun(eNrgFun,meanOtfPd);
% 
% %%
% [x,y] = pol2cart([scs.azimuth],[scs.zenith]*cougarConstants.radian2arcsec);
% z = zeros(size(x));
% tri = delaunay(x,y);
% figure
% subplot(1,2,1)
% trisurf(tri,x,y,z,strehlRatio)
% view(2)
% shading interp
% axis square
% colorbar
% hold on
% polar(scs,'k*')
% polar(gs,'wo')
% hold off
% title('Strehl ratio')
% set(gca,'View',[0 90],'Box','on','xlim',[-1,1]*60,'ylim',[-1,1]*60)
% xlabel('arcsec')
% ylabel('arcsec')
% subplot(1,2,2)
% trisurf(tri,x,y,z,eNrg)
% view(2)
% shading interp
% axis square
% colorbar
% hold on
% polar(scs,'k*')
% polar(gs,'wo')
% hold off
% title('entrapped energy')
% set(gca,'View',[0 90],'Box','on','xlim',[-1,1]*60,'ylim',[-1,1]*60)
% xlabel('arcsec')
% ylabel('arcsec')
% % [o,r] = meshgrid([0:12]*pi/6,[0 30 60]);
% % [x,y] = pol2cart(o,r);
% % u = reshape(strehlRatio(2:end),[],2)';
% % u = [repmat(strehlRatio(1),1,13); u, u(:,1)];
% % v = reshape(eNrg(2:end),[],2)';
% % v = [repmat(eNrg(1),1,13); v, v(:,1)];
% % figure
% % subplot(1,2,1)
% % pcolor(x,y,u)
% % shading interp
% % axis square
% % colorbar
% % subplot(1,2,2)
% % pcolor(x,y,v)
% % shading interp
% % axis square
% % colorbar
% 
% %%
% filename = sprintf('raven-%s-%dscs-%dgs%dMag-#it%d',datestr(now,30),nScs,nGs,gs(1).magnitude,nIt);
% save(filename)
% fprintf(' >> Run saved in %s\n',filename)