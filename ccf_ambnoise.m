% inpu sac files should removed instrument response and down sample

% Calculate ambient noise cross correlation record from multiple stationpairs
%

clear

% parameter_AmbNoise;
setup_parameters;
% comp = parameters.comp;
isfigure2 = 0;
IsFigure = 1;
IsOutput = 1;
% input path
datadir = parameters.datapath;
dist_min = 20;
% output path
ccf_path = './ccf/'

if ~exist(ccf_path)
    mkdir(ccf_path);
end
%% ------------------- loop through center station station-------------------

% PZpath = ['/Volumes/Fraser/SEGMeNT_ambnoise/complete_array/',comp,'_comp/polezero/'];
PZpath = './polezero/';

stalist = parameters.stalist;
nsta=parameters.nsta; % number of target stations to calculate for

for ista1=1:nsta
    
    sta1=char(stalist(ista1,:));
    %ccf_cstadir = [ccf_path,sta1]
    if ~exist([ccf_path,sta1])
        mkdir([ccf_path,sta1]);
    end
    
    %cd([ccf_path,sta1]);
    fpair=fopen([ccf_path,sta1,'/stationpair.txt'],'w');
    fprintf(fpair,'station      pair    lat1     lon1        dep1        lat2      lon2        dep2       distance      azimuth     back azimuth \n');
    fclose(fpair);
    
    % Find the number of days with data
    daylist1 = dir([datadir,sta1,'/*HZ.sac']);
    
    %     temp = strsplit(daylist1(1).name,'.');
    %
    %     netid_temp = temp(2);
    %
    for ista2=1:nsta
        clear lat1 lat2 lon1 lon2 dist az baz
        % -------------------read station coord -------------------
        % sta2
        sta2=char(stalist(ista2,:));
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
        
        % check to see if we've already done this ccf
        if(exist([ccf_path,sta1,'/',sta1,'_',sta2,'_f.mat']))
            display('CCF already exist, skip this pair');
            continue
        elseif exist([ccf_path,sta1,'/',sta2,'_',sta1,'_f.mat'])
            display('CCF already exist, skip this pair');
            continue
        end
        % ---------- get poles zeros values from SAC_pole_Zero file for Z comp ----------
        %         SACPoleZero = dir([INSTRUMENTdir,'SAC_PZs_ZA*',sta1,'*',chz,'*']);
        %         SACPoleZero_Z = [INSTRUMENTdir,SACPoleZero.name];
        %         [zzeros1,ppoles1,ggain1] =  read_sac_pole_zero(SACPoleZero_Z);
        %
        %         SACPoleZero = dir([INSTRUMENTdir,'SAC_PZs_ZA*',sta2,'*',chz,'*']);
        %         SACPoleZero_Z = [INSTRUMENTdir,SACPoleZero.name];
        %         [zzeros2,ppoles2,ggain2] =  read_sac_pole_zero(SACPoleZero_Z);
        %
        
        display(['performing cross-correlation for staion pair : ',sta1,'  ', sta2]);
        % -------------loop through each day--------------------
        nday_stack=0;
        coh_sum = 0;
        coh_num = 0;
        
        % Get a list of all available data
        ihday = 0;
        while ihday <= length(daylist1)-1
            
            ihday = ihday +1;
            clear temp
            temp = strsplit(daylist1(ihday).name,'.');
            
            hdayid = char(temp(1));
            disp(['Looking at ',hdayid,' ',sta2]);
            %------------------- read DATA------------------------
            data1=dir([datadir,sta1,'/',hdayid,'.*.',sta1,'.*HZ.sac']);
            data2=dir([datadir,sta2,'/',hdayid,'.*.',sta2,'.*HZ.sac']);
            
            data1 = [datadir,sta1,'/',data1.name];
            data2 = [datadir,sta2,'/',data2.name];
            
            %------------------- TEST IF DATA EXIST------------------------
            if length(data1) < length(datadir)+20
                disp(['no data sac files for ',hdayid,' ',sta1]);
                continue
            elseif length(data2) < length(datadir)+20
                disp(['no data sac files for ',hdayid,' ',sta2]);
                continue
            end
            
            % Read the sac file
            [vec_tz,Z1raw]=readsac(data1);
            [vec_tz2,Z2raw]=readsac(data2);
            
            % Remove instrument resopnse
            pzfile1 = dir([PZpath,'*',sta1,'*sacpz']);
            pzfile2 = dir([PZpath,'*',sta2,'*sacpz']);
            if length(pzfile1) ~= 1
                pzfile = pzfile1;
%                                 disp('More than one reseponse found!')
                
                % Figure out which response to read
                for ii = 1:length(pzfile)
                    pzdate(ii) = doy2date(str2num(pzfile(ii).name(19:21)),str2num(pzfile(ii).name(14:17)));
                end
                otime = datenum(hdayid,'yyyymmddHHMMSS');
                ind = find(abs(otime-pzdate) == min(abs(otime-pzdate)));
                if ind > length(pzfile) & ind > 1
                    ind = ind(end)-1;
                end
                pzfile = pzfile(ind);
%                 disp(['Using ',pzfile.name])
                pzfile1 = pzfile;
            elseif length(pzfile2) ~= 1
                pzfile = pzfile2;
            
%                 disp('More than one reseponse found!')
                
                % Figure out which response to read
                for ii = 1:length(pzfile)
                    pzdate(ii) = doy2date(str2num(pzfile(ii).name(19:21)),str2num(pzfile(ii).name(14:17)));
                end
                otime = datenum(hdayid,'yyyymmddHHMMSS');
                ind = find(abs(otime-pzdate) == min(abs(otime-pzdate)));
                if ind > length(pzfile) & ind > 1
                    ind = ind(end)-1;
                end
                pzfile = pzfile(ind);
%                 disp(['Using 's,pzfile.name])
                pzfile2 = pzfile;
            end
            dt_new = parameters.dt;
        % Read sacpz file
        [p,z,c] = read_SACPZ([PZpath,pzfile1.name]);  
        
        dt1 = abs(vec_tz(1)-vec_tz(2));
        dt2 = abs(vec_tz2(1)-vec_tz2(2));
        
        % Remove instrument response
        Z1raw = rm_SACPZ(Z1raw,z,p,c,dt1);
        
        [p,z,c] = read_SACPZ([PZpath,pzfile2.name]);
        Z2raw = rm_SACPZ(Z2raw,z,p,c,dt2);
%         [Z2raw,vec_tz2] = resample(odata,dt_new,dt2);
        
            dt = dt_new;
            % Determine the time span to cut to ... this will change with
            % different segments
            clear tcut
            minT1 = min(vec_tz);
            minT2 = min(vec_tz2);
            
            Nstart = 50;
            npts= 86000;
            if (minT1 > minT2) || (minT1 == minT2)
                tcut = [minT1+Nstart*dt:dt:minT1+Nstart*dt+npts];
            elseif minT2 > minT1
                tcut = [minT2+Nstart*dt:dt:minT2+Nstart*dt+npts];
            else
                error('Cant tell which data segment starts first!')
            end
            

            if length(Z1raw) < 20000
                disp(['Sta1 ',sta1,' : ',num2str(length(Z1raw)),' is too short!'])
                continue
            elseif length(Z2raw) < 20000
                disp(['Sta2 ',sta2,' : ',num2str(length(Z2raw)),' is too short!'])
                continue
            end
            
            if(~exist('lat2','var'));
                
                S1 = readsac(data1);
                S2 = readsac(data2);
                
                lat1=S1.STLA;
                lon1=S1.STLO;
                dep1=S1.STEL; % depth is negative for OBS and positive for land stations
                
                
                lat2=S2.STLA;
                lon2=S2.STLO;
                dep2=S2.STEL; % depth is negative for OBS and positive for land stations
                
                
                
                [delta,az]=distance(lat1,lon1,lat2,lon2);
                [delta,baz]=distance(lat2,lon2,lat1,lon1);
                
                dist=deg2km(delta);
                
                Delta=S1.DELTA;
%                 if(abs(Delta-dt) >= 0.01*dt )
%                     error('sampling interval does not match data! check dt');
%                 end
                
                if(dist < dist_min)
                    display('distance shorter than 80 km, skip');
                    break
                end
            end % if lat variables
            
            stapairsinfo.stanames = {sta1,sta2};
            stapairsinfo.lats = [lat1,lat2];
            stapairsinfo.lons = [lon1,lon2];
            
            % cut in time
            Z1=interp1(vec_tz,Z1raw,tcut);
            Z1(isnan(Z1))=0;
            
            Z2=interp1(vec_tz2,Z2raw,tcut);
            Z2(isnan(Z2))=0;
            
            %detrend
            Z1=detrend(Z1);
            Z2=detrend(Z2);
            
            if isfigure2
                figure(1)
                clf
                subplot(2,1,1)
                plot(tcut,Z1,'-r')
                ylim([-0.15e-5 0.15e-5])
                xlim([0 86400])
                hold on
                subplot(2,1,2)
                plot(tcut,Z2,'-b')
                ylim([-0.15e-5 0.15e-5])
                xlim([0 86400])
                hold on
                pause
            end
            
            %             %remove instrumet respsonse
            %             DISP1= rm_SACPZ(Z1,zzeros1,ppoles1,ggain1,dt);
            %             DISP2= rm_SACPZ(Z2,zzeros2,ppoles2,ggain2,dt);
            %             Z1 = DISP1;
            %             Z2 = DISP2;
            
            %despike
            Z1=runwin_norm(Z1);
            Z2=runwin_norm(Z2);
            
            if isfigure2
                figure(2)
                clf
                subplot(2,1,1)
                plot(tcut,Z1,'-r')
                title(['Station 1 : ',sta1])
                hold on
                subplot(2,1,2)
                plot(tcut,Z2,'-b')
                hold on
                title(['Station 2 : ',sta2])
                disp(['Distance : ',num2str(dist)])
                pause
            end
            
            %fft
            fftZ1 = fft(Z1);
            fftZ2 = fft(Z2);
            
            %Whiten
            fftZ1 = spectrumwhiten(fftZ1);
            fftZ2 = spectrumwhiten(fftZ2);
            
            % calculate the cross-correlation
            coh_trace = fftZ1 .* conj(fftZ2);
            coh_trace = coh_trace ./ abs(fftZ1) ./ abs(fftZ2);
            
            % add the cross-correlation to the stack
            coh_sum = coh_sum + coh_trace;
            
            % keep track of the number of cross-correlations in the stack
            coh_num = coh_num + 1;
        end % hday
        
        if coh_num > 1
            if IsFigure
                figure(101);clf;
                set(gcf,'position',[400 400 600 300]);
                dt = 1;
                T = length(coh_sum);
                faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
                ind = find(faxis>0);
                plot(faxis(ind),smooth(real(coh_sum(ind)/coh_num),100));
                title(sprintf('%s %s coherency,station distance: %f km',sta1,sta2,dist));
                xlim([0.04 0.16])
                drawnow
                print(figure(101),'-dpsc',[sta1,'_',sta2,'spec.ps'])
%                 pause
            end
            if IsOutput
                save(sprintf('%s%s/%s_%s_f.mat',ccf_path,sta1,sta1,sta2),'coh_sum','coh_num','stapairsinfo');
                
                fpair=fopen([ccf_path,sta1,'/stationpair.txt'],'a');
                fprintf(fpair,'%s  %5f   %5f  %5f  %5f  %5f  %5f   %5f   %5f   %5f  \n',[sta1,'_',sta2],lat1,lon1,dep1,lat2,lon2,dep2,dist,az,baz);
                fclose(fpair);
            end
        end
    end % ista2
    
    
    
    
    
end % ista1
