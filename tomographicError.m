%% MODAL MOAO
% tomographic error estimation
fromSavedRun = false;
figLabel = 1;
if ~fromSavedRun && matlabpool('size')==0
    matlabpool open
else
    latexInclude = [];
end

%% Atmosphere
Cn2 = [6.39 3.94 1.46 1.73 3.11 2.69 2.81];
fr0 = Cn2/sum(Cn2);
atm = atmosphere(photometry.V,0.15,30,...
    'altitude',[0,0.5,1,2,4,8,16]*1e3,...
    'fractionnalR0',fr0,...
    'windSpeed',[5.6,5.8,6.3,7.6,13.3,19.1,12.1],...
    'windDirection',[0,5,15,30,60,90,180]*pi/180);
atm.wavelength = photometry.H;

%% Zernike for sensing and correction
for radialDegree = 4
    
    fovInArcmIn = 3.5;
    zern = zernike(2:zernike.nModeFromRadialOrder(radialDegree),8,'fieldOfViewInArcmin',fovInArcmIn);
    pistonFreeVariance = phaseStats.zernikeResidualVariance(1,atm,zern);
    
    %% Science direction
    logs = logBook.checkIn;
    for nGs = 4

        scDirZenith = (0:7.5:fovInArcmIn*30);
        if fromSavedRun
            
            nGsRing =  6;
            nDir    = 15;
            files = dir('*.mat');
            fileEnd  = sprintf('-%dGS-%drd-%drings-%ddirs',nGs,radialDegree,nGsRing,length(scDirZenith));
            files    = files( cellfun( @(x)~isempty(x) , regexp({files.name},fileEnd) ) );
            filename = files(end).name;
            fprintf(' ==> Loading %s ...!\n',filename)
            logs.verbose = false;
            load(filename)
            logs.verbose = true;
            fileEnd = sprintf('-1+3GS-%drd-%drings-%ddirs',radialDegree,nGsRing,length(scDirZenith));
            fromSavedRun = true;
        else
            
            nDir = length(scDirZenith);
            scDir1 = source('zenith',scDirZenith*constants.arcsec2radian,...
                'azimuth',zeros(size(scDirZenith)),...
                'wavelength',photometry.H);
            scDir2 = source('zenith',scDirZenith*constants.arcsec2radian,...
                'azimuth',ones(size(scDirZenith))*pi/(nGs-1),...
                'wavelength',photometry.H);
            
            %% Guide Stars
            gsZenith = 0.25:0.25:(fovInArcmIn/2-0.25);
            nGsRing = length(gsZenith);
            count = 0;
            residualVariance = zeros(nDir,2*nGsRing);
            residualVarianceSCAO = zeros(nDir,2*nGsRing);
            logs.verbose = false;
            
            fprintf('_____________________________________________________\n')
            fprintf('Tomographic error for:\n . %d GS,\n . %d Zernike radial degrees,\n . %d GS rings,\n . %d science dir.\n',...
                nGs,radialDegree,nGsRing,length(scDirZenith))
            fileEnd = sprintf('-%dGS-%drd-%drings-%ddirs',nGs,radialDegree,nGsRing,length(scDirZenith));            
            
            for kGs = 1:nGsRing
                tic
                
                fprintf(' ==> GS Ring Radius: %4.2f''',gsZenith(kGs))
                gs = source('asterism',{[0,0],[nGs-1,gsZenith(kGs)*constants.arcmin2radian,0]},'wavelength',photometry.H);
                %% Sensing covariance
                S = phaseStats.zernikeAngularCovariance(zern,atm,gs);
                S1 = S{1};
                S = cell2mat(S);
                R = chol(S);
                
                %% Sensing/Correction Covariance
                C1 = phaseStats.zernikeAngularCovariance(zern,atm,gs,scDir1);
                C2 = phaseStats.zernikeAngularCovariance(zern,atm,gs,scDir2);
                
                %% Optimal command matrix
                count = count + 1;
                for kDir = 1:nDir
                    CDir = cell2mat(C1(:,kDir));
                    residualVariance(kDir,count)     = ...
                        pistonFreeVariance - trace(CDir'/S*CDir);
                    CDir = cell2mat(C1(1,kDir));
                    residualVarianceSCAO(kDir,count) = ...
                        pistonFreeVariance - trace(CDir'/S1*CDir);
                    CDir = cell2mat(C2(:,kDir));
                    residualVariance(kDir,count+1)   = ...
                        pistonFreeVariance - trace(CDir'/S*CDir);
                    CDir = cell2mat(C2(1,kDir));
                    residualVarianceSCAO(kDir,count+1) = ...
                        pistonFreeVariance - trace(CDir'/S1*CDir);
                end
                count = count + 1;
                
                elapsedTime = toc;
                fprintf(' - Elapsed time: %6.2f\n',elapsedTime)
            end
            fprintf('_____________________________________________________\n')
            
            logs.verbose = true;
            
        end
        %%
        var2StrokeInMicron = @(x) 1e6*sqrt(x)*atm.wavelength/2/pi;
        rmsMicron = var2StrokeInMicron(residualVariance);
        rmsMicronSCAOmin = var2StrokeInMicron(phaseStats.zernikeResidualVariance(zern.j(end),atm,zern));
        rmsMicronSCAO = var2StrokeInMicron(residualVarianceSCAO);

        figure(1)
        plot(scDirZenith,rmsMicron(:,1:2:end))
        hold on
        plot(scDirZenith,rmsMicronSCAO(:,1:2:end),':')
        hold off
        line(scDirZenith([1,end]),ones(1,2)*rmsMicronSCAOmin,'color','k')
        set(gca,'xtick',gsZenith*60,'xlim',[0,fovInArcmIn/2]*60)
        grid
        xlabel('Zenith [arcsec]')
        ylabel('Residual rms [micron]')
        legend(num2str(gsZenith'*60),'location','EastOutside')
        graphFile1 = ['inGs-',fileEnd];
        
        figure(2)
        plot(scDirZenith,rmsMicron(:,2:2:end))
        hold on
        plot(scDirZenith,rmsMicronSCAO(:,2:2:end),':')
        hold off
        set(gca,'xtick',gsZenith*60,'xlim',[0,fovInArcmIn/2]*60)
        grid
        xlabel('Zenith [arcsec]')
        ylabel('Residual rms [micron]')
        legend(num2str(gsZenith'*60),'location','EastOutside')
        graphFile2 = ['midGs-',fileEnd];
        
        figure(3)
        plot(gsZenith*60,rmsMicron(1,1:2:end),':o',...
            gsZenith*60,rmsMicron(1,2:2:end),':x',...
            gsZenith*60,rmsMicron(9,1:2:end),':p',...
            gsZenith*60,rmsMicron(9,2:2:end),':h',...
            gsZenith*60,rmsMicron(end,1:2:end),':d',...
            gsZenith*60,rmsMicron(end,2:2:end),':s')
        grid
        set(gca,'xtick',gsZenith*60)
        xlabel('GS ring radius [arcmin]')
        ylabel('Residual rms [micron]')
        legend('0''/In-GS','0''/Mid-GS',...
            '1''/In-GS','1''/Mid-GS',...
            '1.75''/In-GS','1.75''/Mid-GS',...
            'location','EastOutside')
        graphFile3 = ['ringGs-',fileEnd];
        figLabel = figLabel + 1;

        if fromSavedRun
            saveas(1,graphFile1,'png')
            latexInclude = sprintf(['%s\n%%%s\n\\begin{figure}\n\\centering\n',...
                '\\includegraphics[width=0.45\\linewidth]{%s}\n'],...
                latexInclude,graphFile1,graphFile1);
            saveas(2,graphFile2,'png')
            latexInclude = sprintf('%s%%%s\n\\includegraphics[width=0.45\\linewidth]{%s}\n',...
                latexInclude,graphFile2,graphFile2);
            saveas(3,graphFile3,'png')
            latexInclude = sprintf(['%s%%%s\n\\includegraphics[width=0.45\\linewidth]{%s}\n',...
                '\\caption{%d GSs and %d Zernike radial order}\n\\label{fig:%d}\n\\end{figure}\n'],...
                latexInclude,graphFile3,graphFile3,nGs,radialDegree,figLabel);
            fid = fopen('graphs.tex','w');
            fprintf(fid,'%s',latexInclude);
            fclose(fid);
        end
        
        %%
%         filename = sprintf('raven-%s-%dGS-%drd-%drings-%ddirs',datestr(now,30),nGs,radialDegree,nGsRing,length(scDirZenith));
%         save(filename)
%         fprintf(' >> Run saved in %s\n',filename)
        
    end % NGS loop
    
end % Zernike radial degree loop
