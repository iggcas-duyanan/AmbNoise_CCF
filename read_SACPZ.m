% Read SACPZ files
%
% Expects a SAC polezero formatted file
% NJA, 3/30/2016
%
% p = poles
% z = zeros
% c = constant
function [p,z,c] = read_SACPZ(filename)

% filename = './polezero/YQ.103B..HH1.2015_064_14_55_00.sacpz';
fid = fopen(filename,'r');

D = textscan(fid,'%s%f\n',1);

z = zeros(D{2},1);
for iz = 1:D{2}
    temp = textscan(fid,'%f%f\n',1);
    z(iz) = double(temp{1})+i*double(temp{2});
end
D = textscan(fid,'%s%f\n',1);

p = zeros(D{2},1);

for ip = 1:D{2}
    temp = textscan(fid,'%f%f\n',1);
    
    p(ip) = double(temp{1})+i*double(temp{2});
end

D = textscan(fid,'%s%f',1);

c = double(D{2});

fclose(fid);