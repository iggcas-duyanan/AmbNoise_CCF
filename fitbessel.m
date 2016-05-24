% Program to fit the xcorf and get the travel time for each frequency for two station pairs

clear

global tN
global waxis
global twloc
global weight


setup_parameters;
% npts = parameters.npts;
npts = parameters.npts; % number of hours in the windows
IsFigure = 1;
isplotinit = 0;

isfigure2 = 0;


frange = [0.04 0.125];
% frange = [0.033 0.2];

tN = 9;

isoutput = 1;
nearstadist = 20;

%refc1 = 4.2;
%refc2 = 1.5;
%refc = refc1 + (refc2-refc1)/(tN)*[1:tN];

wholesec = npts;

twloc = frange(1):(frange(2)-frange(1))/(tN-1):frange(2);
t_vec=1./twloc;


%% Make an initial model
% From GOC_CC card
vec_h = [5 10 20 20 20];
vec_vs = [2.5 3.3 3.6 4.2 4.4];
% vec_vs = [2.5 3.3 3.6 4.1 4.336];
vec_vp = vec_vs.*1.8;
vec_rho = [2.7 3.027 3.027 3.352 3.359];

refmod(:,1) = vec_h(:);
refmod(:,2) = vec_vp(:);
refmod(:,3) = vec_vs(:);
refmod(:,4) = vec_rho(:);
refmod

% Calculate the predicted dispersion curve
[c, vg, grad_c, grad_vg] = Calc_Ray_dispersion(t_vec,refmod,1,0);
wvec1 = (2*pi)./t_vec;
wvec1 = wvec1';
refc = c';

if isplotinit
t_vec1 = interp(t_vec,4);
[c1,vg1,grad_c1,grad_vg1] = Calc_Ray_dispersion(t_vec1,refmod,1,0);
wvec1 = (2*pi)./t_vec1;
wvec1 = wvec1';
end

% input path
ccf_path = './ccf/';

% output path
XSP_path =  './Xsp/';

if ~exist(XSP_path)
    mkdir(XSP_path)
end


warning off; %#ok<WNOFF>




% Get your axis correct
twloc = twloc*2*pi;
waxis = (frange(1):1/wholesec:frange(2))*2*pi;

faxis = [0:wholesec/2 -1]*1/wholesec;

stalist = parameters.stalist;
nsta=parameters.nsta; % number of target stations to calculate for

%%% --- Loop through station 1 --- %%%
for ista1=1:nsta
    
    sta1=char(stalist(ista1,:));
    sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
    
    %%% --- Loop through station 2 --- %%%
    for ista2 = 1: nsta % length(v_sta)
        sta2 = char(stalist(ista2,:));
        
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
        
        % Check to see if we have already done this
        if exist([XSP_path,sta1,'_',sta2,'_xsp.mat'])
            disp('Already fit this one!')
            continue
        end
        clear data1 xcorf1 xsp1 filename
        
        %%% --- Load in the ccf --- %%%
        filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);
        
        if ~exist(filename,'file')
            disp(['not exist ',filename])
            continue;
        end
        
        
        data1 = load(filename);
        delta = distance(data1.stapairsinfo.lats(1),data1.stapairsinfo.lons(1),data1.stapairsinfo.lats(2),data1.stapairsinfo.lons(2));
        r1    = deg2km(delta); % distance
        
        if r1 < nearstadist
            continue;
        end
        
        
        %%% - Get the normalized ccf - %%%
        xcorf1 = data1.coh_sum./data1.coh_num;
        dumnan = find(isnan(xcorf1)==1);
        
        if length(dumnan) > 10
            disp([sta1,' and ',sta2,'is NaN! Moving on']);
            continue
        end
        
        N = 10000;
        if length(xcorf1) < N
            disp('Dataset is too short! Moving on')
            continue
        end
        
        xcorf1 = real(xcorf1(1:N));
        xcorf1(1) = 0;
        
        if isfigure2 
            figure(1)
            T = length(xcorf1);
            dt = 1;
            temp_faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
            ind = find(temp_faxis>0);
            subplot(2,1,1)
            plot(temp_faxis(ind),smooth(real(xcorf1(ind)),100));
            xlim([frange(1) frange(2)])
            hold on
            subplot(2,1,2)
            plot(temp_faxis(ind),smooth(real(data1.coh_sum(ind)/data1.coh_num),100),'-r')
            xlim([frange(1) frange(2)])
            
        end

        %%% - Convert xcorf into spherical frequency - %%%
        faxis = [0:N-1]*1/wholesec;
        xsp1 = interp1(faxis*2*pi,xcorf1,waxis);

        xsp1 = smooth(xsp1,100);

        tw1 = ones(1,tN)*r1./refc;
        
        %%% - Invert for the bessel function 2x - %%%
        weight  = 1./waxis;
        tw2 = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw1],[tw1]*0.8,[tw1]*1.2);
        
        weight(:) = 1;
        tw = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw2],[tw2]*0.8,[tw2]*1.2);
        
        
        %%% - Set up the variable structure - %%%
        xspinfo.sta1 = sta1;
        xspinfo.sta2 = sta2;
        xspinfo.lat1 = data1.stapairsinfo.lats(1);
        xspinfo.lon1 = data1.stapairsinfo.lons(1);
        xspinfo.lat2 = data1.stapairsinfo.lats(2);
        xspinfo.lon2 = data1.stapairsinfo.lons(2);
        
        xspinfo.r = r1;
        xspinfo.tw = tw;
        xspinfo.xsp = xsp1;
        xspinfo.coherenum = data1.coh_num;
        err = besselerr(tw,xsp1);
        err = err(1:length(waxis));
        xspinfo.sumerr = sum(err.^2)./sum((xsp1./weight(:)).^2);
        xspinfo.err = err./weight(:);
        xspinfo.tw1 = tw1;
        xspinfo.twloc = twloc;
        
        data = r1./tw;
        
        xcorf1 = data1.coh_sum./data1.coh_num;

        %%% Calculate SNR %%%
        snrdata = real(ifft((xcorf1)));
        snrdata = fftshift(snrdata);
        groupv_max = 4.5;
        groupv_min = 2.0;
        NN= length(snrdata);
        lag = [-NN/2:NN/2];
        win_min = r1./groupv_max;
        win_max = r1./groupv_min;
        if win_min < 15; win_min = 0; end
        if win_max < 50; win_max = 50; end
        signal_ind = find((lag>-win_max & lag<-win_min) | (lag>win_min & lag<win_max));
        signal_amp = sum(snrdata(signal_ind).^2)/length(signal_ind);
        noise_amp = sum(snrdata.^2)/length(snrdata);
        snr = signal_amp/noise_amp;

        xspinfo.filename = filename;
        xspinfo.snr = snr;
        
        % Calculate the predicted bessel function from the initial model
        A = 1;
        binit = besselj(0,(wvec1.*r1)./c)*A;
        binit = binit./mean(abs(binit)).*mean([abs(xsp1)]);
        disp([filename,' fitted'])
        if IsFigure
            besselerr(tw,xsp1,3);
            if isplotinit
                plot(wvec1/2/pi,binit,'-k')
            end
            hold on
            subplot(2,1,2);
            plot(1./(twloc/2/pi),r1./tw1,'ko-');hold on;
            plot(1./(twloc/2/pi),r1./tw,'ro-');
            title([sta1,'-',sta2])
            
            subplot(2,1,1)
            xlim([frange(1) frange(2)])
            
            if isfigure2
            figure(12)
            clf
            T = length(data1.coh_sum);
            dt = 1;
            faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
            ind = find(faxis>0);
            plot(faxis(ind),smooth(real(data1.coh_sum(ind)/data1.coh_num),100));
            xlim([frange(1) frange(2)])
            end
            psfile = [XSP_path,'Xsp_',sta1,'_',sta2,'.ps'];
            print('-dpsc2',psfile);
            
            
%             pause;
        end
        if isoutput
            save(sprintf('%s/%s_%s_xsp.mat',XSP_path,sta1,sta2),'xspinfo','twloc','waxis');
        end
        
        
    end %end of station j
end  %end of station i
%stapairn
%%soundsc(rand(2000,1),1000,8)
