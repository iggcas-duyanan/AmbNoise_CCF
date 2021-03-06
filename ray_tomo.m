% Script to do the ray theory tomography based on the ambient noise measurement
% Written by Ge Jin, jinwar@gmail.com
% Nov 2012
%
% Modified by NJA, April 2016
clear

% Load color scale
load seiscmap.mat

% Set up geometry parameters
setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx = length(xnode);
Ny = length(ynode);

% Save results?
isoutput = 1;

% Set up error parameters
errlevel = parameters.errlevel;
snrtol = parameters.snrtol;
mincoherenum = parameters.mincoherenum;
fiterrtol = parameters.fiterrtol;

refv = parameters.refv;
distrange = parameters.distrange;

maxerrweight = parameters.maxerrweight;
polyfit_dt_err = parameters.polyfit_dt_err;
smweight0 = parameters.smweight0;


dterrtol = parameters.dterrtol;

raydensetol = parameters.raydensetol;
r = parameters.r;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
Nx = length(xnode);
Ny = length(ynode);

savefile = parameters.savefile;

% read in bad station list, if existed
if exist('badsta.lst')
    badstnms = textread('badsta.lst','%s');
    badstaids = find(ismember({stainfo.staname},badstnms));
    disp('Found Bad stations:')
    disp(badstnms)
end

% Set up initial smoothing kernel
[i,j] = ndgrid(1:Nx,2:(Ny-1));
ind = j(:) + Ny*(i(:)-1);
dy = diff(ynode);
dy1 = dy(j(:)-1);
dy2 = dy(j(:));
Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
    [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), -2./(dy2.*(dy1+dy2))],Nx*Ny,Nx*Ny);
[i,j] = ndgrid(2:(Nx-1),1:Ny);
ind = j(:) + Ny*(i(:)-1);
dx = diff(xnode);
dx1 = dx(i(:)-1);
dx2 = dx(i(:));
Areg = [Areg;sparse(repmat(ind,1,3),[ind-Ny,ind,ind+Ny], ...
    [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), -2./(dx2.*(dx1+dx2))],Nx*Ny,Nx*Ny)];
F=Areg;

% Initialize the xsp structure
Xsp_path = './Xsp/';
xspfiles = dir([Xsp_path,'*_*xsp.mat']);

disp('Looking at Xsp Files')
for ixsp = 1:length(xspfiles)
    
    temp = load([Xsp_path,xspfiles(ixsp).name]);
    xspinfo = temp.xspinfo;
    
    if ixsp ==1
        Tperiods = (2*pi)./temp.twloc;
        waxis = temp.waxis;
        twloc = temp.twloc;
        xspinfo.isgood = 0;
        xspsum = xspinfo;
    else
        xspinfo.isgood = 0;
        xspsum = [xspsum;xspinfo];
    end
    clear temp

    
    % 	xspinfo(ixsp).isgood = 0;
    if xspsum(ixsp).sumerr < errlevel ...
            && xspsum(ixsp).snr > snrtol && xspsum(ixsp).coherenum > mincoherenum
        xspsum(ixsp).isgood = 1;
    end
    
    if rem(ixsp,500)==0
        disp(['Looking at #',num2str(ixsp),' of ',num2str(length(xspfiles))])
    end
end % end of loop ixsp'


% Loop through periods
for ip=1:length(Tperiods)
    disp(' ');
    disp(['Inversing Period: ',num2str(Tperiods(ip))]);
    clear rays dt fiterr mat phaseg err raydense dist
    raynum = 0;

    for ixsp = 1:length(xspsum)
        if xspsum(ixsp).isgood ==0;
            continue;
        end
        if xspsum(ixsp).r > refv*Tperiods(ip)*distrange(2)...
                || xspsum(ixsp).r < refv*Tperiods(ip)*distrange(1)
            continue;
        end
        
        raynum = raynum+1;
        rays(raynum,1) = xspsum(ixsp).lat1;
        rays(raynum,2) = xspsum(ixsp).lon1;
        rays(raynum,3) = xspsum(ixsp).lat2;
        rays(raynum,4) = xspsum(ixsp).lon2;
        
        dist(raynum) = deg2km(distance(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4)));
        dt(raynum) = xspsum(ixsp).tw(ip);
        err = smooth((abs(xspsum(ixsp).err)./mean(abs(xspsum(ixsp).xsp))).^2,round(length(waxis)/length(twloc)));
        fiterr(raynum) = interp1(waxis(:),err(:),twloc(ip)); 
        csnum(raynum) = xspsum(ixsp).coherenum;
        snr(raynum) = xspsum(ixsp).snr;
        errays(raynum,1) = xspsum(ixsp).lat1;
        errays(raynum,2) = xspsum(ixsp).lon1;
        errays(raynum,3) = xspsum(ixsp).lat2;
        errays(raynum,4) = xspsum(ixsp).lon2; 
        errays(raynum,5) = fiterr(raynum);
    end
    if size(dt,1) ~=raynum
        dt = dt';
    end
    
    % Building the data kernel
    disp('Start building the kernel');
    tic
    mat=ray_kernel_build(rays,xnode,ynode);  
    toc
    % Calculate the weighting matrix
    W = sparse(length(dt),length(dt));
    for i=1:length(dt)
        W(i,i)=1./fiterr(i);
    end
    ind = find(W > maxerrweight);
    W(ind) = maxerrweight;
    ind = find(W < 1/fiterrtol);
    W(ind) = 0;
    for i=1:length(dt)
        W(i,i)=W(i,i).*(csnum(i).^0.5);
    end
    para = polyfit(dist(:),dt,1);
    polyerr = polyval(para,dist(:)) - dt;
    errind = find(abs(polyerr) > polyfit_dt_err);
    for i = errind
        W(i,i) = 0;
    end
    
    % calculate the smoothing weight
    smweight = smweight0;
    NR=norm(F,1);
    NA=norm(W*mat,1);
    smweight = smweight0*NA/NR;
    
    disp('start inverse');
    A=[W*mat;smweight*F];
    rhs=[W*dt;zeros(size(F,1),1)];

    phaseg=(A'*A)\(A'*rhs);
    %        toc
    %        disp('Done');
    
    
    % Iteratively down weight the measurement with high error
    niter=1;
    
    while niter < 2
        niter=niter+1;
        err = mat*phaseg - dt;

        stderr=std(err);
        if stderr > dterrtol
            stderr = dterrtol;
        end
        ind = find(diag(W)==0);
        disp('Before iter:');
        disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
        disp(['Bad Measurement Number: ', num2str(length(ind))]);
        for i=1:length(err)
            if abs(err(i)) > 2*stderr
                W(i,i)=0;
            end
        end
        ind = find(diag(W)==0);
        disp('After iter:');
        disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
        disp(['Bad Measurement Number: ', num2str(length(ind))]);
        
        % Rescale the smooth kernel
        NR=norm(F,1);
        NA=norm(W*mat,1);
        smweight = smweight0*NA/NR;
        
        % Invert
        A=[W*mat;smweight*F];
        rhs=[W*dt;zeros(size(F,1),1)];

        phaseg=(A'*A)\(A'*rhs);

        
    end
    
    %        disp(' Get rid of uncertainty area');
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            raydense(i,j) = sum(mat(:,n));
            if raydense(i,j) < raydensetol
                phaseg(n)=NaN;
            end
        end
    end
    
    % Convert into phase velocity
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            GV(i,j)= 1./phaseg(n);
        end
    end

    raytomo(ip).GV = GV;
    raytomo(ip).mat = mat;
    raytomo(ip).raydense = raydense;
    raytomo(ip).period = Tperiods(ip);
    raytomo(ip).w = diag(W);
    raytomo(ip).err = err;
    raytomo(ip).rays = rays;
    raytomo(ip).fiterr = fiterr;
    raytomo(ip).dt = dt;
    raytomo(ip).smweight0 = smweight0;
    
    
    if isFigure
        figure(1)
        clf
        ax = worldmap(lalim, lolim);
        set(ax, 'Visible', 'off')
        surfacem(xi,yi,raytomo(ip).GV);
        drawlocal
        title([num2str(Tperiods(ip))],'fontsize',15)
        avgv = nanmean(raytomo(ip).GV(:));
        caxis([avgv*(1-r) avgv*(1+r)])
        colorbar
        colormap(seiscmap)
    end
    
end % end of period loop

lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
isoutput = 1;
if isoutput
    save(savefile,'raytomo','xnode','ynode');
    save('coor.mat','xi','yi','xnode','ynode','gridsize','lalim','lolim');
end

Mp = 4; Np = 3;
figure(17)
clf
for ip=1:length(Tperiods)
    subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,raytomo(ip).GV);
    drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
    avgv = nanmean(raytomo(ip).GV(:));
    caxis([avgv*(1-r) avgv*(1+r)])
    colorbar
    colormap(seiscmap)

end

figure(18)
clf

for ip=1:length(Tperiods)
subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,raytomo(ip).raydense);
    drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
    colorbar
    caxis([0 3000])

end

figure(19)
clf


for ip=1:length(Tperiods)
    subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    clear rays
    rays = raytomo(ip).rays;
%     surfacem(xi,yi,raytomo(ip).err);
%     drawpng
scatterm((rays(:,1)+rays(:,3))./2,(rays(:,2)+rays(:,4))./2,30,raytomo(ip).fiterr,'filled')
drawlocal
title([num2str(Tperiods(ip))],'fontsize',15)
colorbar
caxis([ 0 5])


end

