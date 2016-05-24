% Plot the cross-spectra in the time domain for the individual station pairs
%
% NJA, 04/08/2016

clear all;

setup_parameters;
stalist = parameters.stalist;
nsta = parameters.nsta;
% input path
% ccf_path = './ccf_6hour_stack/';
LBLFNT = 17;


coperiod = [ 8 30 ];

% Plot the velocities of the fundamental mode and overtones
% vel_1st = 3.67;
% per_1st = 9;
% vel_0st = 2.79;
% per_0st = 10.3;

% 6 seconds
% vel_1st = 2.96;
% per_1st = 6.66;
% vel_0st = 2.81;
% per_1st = 6.66;

% 19 s
% vel_1st = 4.0;
% per_1st = 19;
% vel_0st = 3.14;
% per_1st = 19;

f1 = 1/coperiod(2);
f2 = 1/coperiod(1);
[b a] = butter(2,[f1 f2]);

% range of group velocity to use for windoing CCF in km/s
cmax=7.5;
cmin=2.5;

ccf_path = './ccf/';
% Get the color information
nsta=length(stalist); % number of target stations to calculate for
npairall = 0;
% nsta1 =
for ista1=1:6
    
    sta1=char(stalist(ista1,:));
    
    sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
    
    nstapair = 0;
    for ista2 = 1: nsta
        sta2 = char(stalist(ista2,:));
        
        % Get color of station pair depending on the station pair
        clr = color_station(sta2);
        
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
        %clear data1 xcorf1 xsp1 filename
        
        filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);
        
        if ~exist(filename,'file')
            disp(['not exist ',filename])
            continue;
        end
        nstapair = nstapair + 1;
        data1 = load(filename);
        xcorf1 = data1.coh_sum./data1.coh_num;
        
        N = length(xcorf1);
        
        
        
        xcorifft1 = real(ifft(2*xcorf1([1:N/2+1]),N));
        %rearrange and keep values corresponding to lags: -(len-1):+(len-1)
        xcorifft2 = [xcorifft1(end-N+2:end) ; xcorifft1(1:N)];
        
        
        xcor(nstapair,:) =  filtfilt(b,a,xcorifft2);
        
        % keep track of colors
        clrall(nstapair,:) = clr;
        xcor(nstapair,:) = xcor(nstapair,:)/max(abs(xcor(nstapair,:)));
        stapairdist(nstapair) = deg2km(distance(data1.stapairsinfo.lats(1),data1.stapairsinfo.lons(1),data1.stapairsinfo.lats(2),data1.stapairsinfo.lons(2)));
        
        stapairinv = [sta2,'_',sta1];
        
        if exist('existpair','var')
            if find(strncmp(stapairinv,existpair,length(stapairinv)))
                continue
            end
        end
        dumsta2{nstapair} = sta2;
        npairall = npairall + 1;
        xcor_all(npairall,:) = xcor(nstapair,:) ;
        %xcorf_all(npairall,:) = real(xcorf1);
        pairdistace_all(npairall) = stapairdist(nstapair);
        existpair(npairall) = {[sta1,'_',sta2]};
        
        
    end % ista2
    
    
    
    
    IsFigure = 1;
    if IsFigure
        f101 = figure(101);
        clf
        hold on;
        N= length(xcorifft2);
        time = [-N/2:N/2];
        amp = 1e1;
        indtime = find(abs(time)<=500);
        set(gca,'YDir','reverse');
        for istapair = 1: nstapair
            plot(time(indtime(1):indtime(end)),xcor(istapair,indtime(1):indtime(end))*amp+stapairdist(istapair),'color',clrall(istapair,:)); hold on;
            %             text(0,stapairdist(istapair),dumsta2{istapair})
            %             return
            %             pause
        end
        xlim([-300 300])
        title(['reference station:',sta1,'  filtered in ',num2str(coperiod(1)), ' - ',num2str(coperiod(2)),'(s)'],'fontname','Times New Roman', 'fontsize',LBLFNT);
        %         return
        print(f101,'-dpsc',['ccf_',sta1,'.ps']);
        %         pause
    end
    % Plot lines of group velocity for hte first overtone and fundamental mode.
    % x1st = stapairdist/vel_1st;
    % x0st = stapairdist/vel_0st;
    % ydum = stapairdist;
%     figure(101)
%     hold on
    % plot(x1st,ydum,'-r','linewidth',1)
    % plot(x0st,ydum,'-b','linewidth',1)
    % plot(x1st*-1,ydum,'-r','linewidth',1)
    % plot(x0st*-1,ydum,'-b','linewidth',1)
    %     pause
end % ista1




N= length(xcorifft2);
time = [-N/2:N/2];
amp = 1e1;
indtime = find(abs(time)<=500);

f102 = figure(102);
clf
hold on;
set(gca,'YDir','reverse');
for istapair = 1: npairall
    plot(time(indtime(1):indtime(end)),xcor_all(istapair,indtime(1):indtime(end))*amp+pairdistace_all(istapair),'k');
end
xlim([-200 200])
title(['All non-repeated pairs, filtered in ',num2str(coperiod(1)), ' -',num2str(coperiod(2)),'(s)'],'fontname','Times New Roman', 'fontsize',LBLFNT);

print(f102,'-dpsc','all_ccf.ps');
