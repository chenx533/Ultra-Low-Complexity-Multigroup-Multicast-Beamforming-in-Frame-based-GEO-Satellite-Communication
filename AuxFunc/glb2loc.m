function [laz,lel] = glb2loc(gaz,gel)
% 全局坐标系下的俯仰角转换为局部坐标系下的俯仰角
gx=cosd(gel)*cosd(gaz);
gy=cosd(gel)*sind(gaz);
gz=sind(gel);
lx=-gz;
ly=gy;
lz=gx;
laz=mod(rad2deg(atan(ly/lx))+180*(lx<0),360);
lel=rad2deg(atan(lz/sqrt(lx^2+ly^2)));
end