% Script to do the ray theory tomography based on the ambient noise measurement
% Written by Ge Jin, jinwar@gmail.com
% Nov 2012


clear

% input files
% load stainfo_BHZ.mat
% load xspinfo.mat
load seiscmap.mat

% Set up geometry parameters

setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
r = parameters.r;


xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx = length(xnode);
Ny = length(ynode);
isoutput = 1;

savefile = parameters.savefile;
% Xsp_path = './XspTT/';

temp = load(savefile);
raytomo = temp.raytomo;

% Get periods
for ip = 1:length(raytomo)
    Tperiods(ip) = raytomo(ip).period;
end

lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);


for ip=1:length(Tperiods)
    
    figure(17)
    clf
    subplot(2,2,1)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,raytomo(ip).GV);
    drawlocal
    title(['PhaseV: ',num2str(Tperiods(ip))],'fontsize',15)
    avgv = nanmean(raytomo(ip).GV(:));
    caxis([avgv*(1-r) avgv*(1+r)])
    colorbar
    colormap(seiscmap)
    
    subplot(2,2,2)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,raytomo(ip).raydense);
    drawlocal
    title(['Rays: ',num2str(Tperiods(ip))],'fontsize',15)
    colorbar
    caxis([0 3000])
    
    subplot(2,2,3)
        ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    clear rays
    rays = raytomo(ip).rays;

    scatterm((rays(:,1)+rays(:,3))./2,(rays(:,2)+rays(:,4))./2,30,raytomo(ip).fiterr,'filled')
    drawlocal    
    title(['RayErr: ',num2str(Tperiods(ip))],'fontsize',15)
    colorbar
    caxis([ 0 5])
    
    if Tperiods(ip) < 10
        print(figure(17),'-dpsc',['Raytomo_0',num2str(Tperiods(ip)),'.ps']);
    else
       
    print(figure(17),'-dpsc',['Raytomo_',num2str(Tperiods(ip)),'.ps']);
    end
end



