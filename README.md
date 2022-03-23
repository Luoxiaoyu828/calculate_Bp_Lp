# calculate_Bp_Lp
计算手绘黑子图的经纬度
"""% 计算Bp，Lp
% 输入手绘黑子图太阳圆半径radii，以圆心为原点的坐标系的坐标(Xa,Ya),
% 及三个参数P,B0,L0(单位：度)，太阳角半径appdiam（单位：度）
% 输出坐标(Xa,Ya)对应的经度Lp,纬度Bp;
% 算法来源：http://www.petermeadows.com/html/location.html
% 算例验证通过
% 算例：for example, if on the 1st January 1999 at 11h 10m, when Bo = -3.0°, Lo = 139.5°, P = +2.1°
% and the apparent solar diameter = 32' 35", a sunspot was measured to be at x = -27mm and y = -22mm on a 150mm diameter disk drawing,
% then rou1 = 0.13°, rou = 27.5° and = -129.2° giving B = -20.6° and L = 161.3°.
%
% [B0,L0,P]=GetB0L0PbyDate(1999,1,1,11,10);
% appdiam= 32/60+35/3600;
% Xa=-27;
% Ya=-22
% radii=150/2;
%  [Bp,Lp]=GetBLbyXYPBLR(Xa,Ya,P,B0,L0,radii,appdiam)
