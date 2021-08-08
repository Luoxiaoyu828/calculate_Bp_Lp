from math import *
import time

def get_bl_by_xyblr(xa, ya, P, B0, L0, radii, appdiam):
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
"""
    x = xa
    y = ya
    sita = degrees(atan(xa/ya)) - 180
    diskdiam = 2 * radii
    # 单位均为：度
    rou1 = (appdiam/diskdiam)*sqrt(xa * xa + ya * ya)
    ra = 2 * rou1 / appdiam
    if ra > 1:
        ra = 1
    elif ra < -1:
        ra = -1

    rou = degrees(asin(ra) - rou1)
    # if x < 0 and y < 0:
    #     sita = math.degrees(math.atan(xa/ya)) - 180
    # if x < 0 and y > 0:
    #     sita = math.degrees(math.atan(xa / ya))
    # if x > 0 and y > 0:
    #     sita = math.degrees(math.atan(xa / ya))
    # if x > 0 and y < 0:
    #     sita = math.degrees(math.atan(xa / ya)) - 180

    if y < 0:
        sita = degrees(atan(xa/ya)) - 180
    else:
        sita = degrees(atan(xa / ya))
    B0_ = radians(B0)
    rou_ = radians(rou)
    P_ = radians(P)
    site_ = radians(sita)
    Bp = degrees(asin(sin(B0_) * cos(rou_) + cos(B0_) * sin(rou_) * cos(P_ - site_)))
    Bp_ = radians(Bp)
    # Lp = (degrees(asin(sin(rou_)*sin(P_-site_)/cos(Bp_))) + L0) % 360
    Lp = degrees(asin(sin(rou_) * sin(P_ - site_) / cos(Bp_))) % 360

    return [Bp, Lp]


def get_b0l0p_by_data(year, month, day, hour, minute):
    """% 根据日期/时间计算日照坐标B0,L0和方位角P ；
    % 输出:  日照坐标B0,L0和方位角P(单位：度)
    %   P:从太阳圆盘北点向东测量的旋转轴北端的位置角;
    %     B0:日照纬度
    %   L0:日照经度
    % 输入：年，月，日，时，分 (世界标准时间，UT)

    % 参考文献：
    % Astronomical Algorithms 第22、25和29章
    """
    # 将时和分转换为天
    day = day + hour/24 + minute/24/60
    # 计算儒略日数
    JD = GetTbyDate(year, month, day)
    T = (JD - 2451545.0) / 36525
    LL = (280.46646 + 36000.76983 * T + 0.0003032 * T **2) % 360
    M = (357.52911 + 35999.05029 * T - 0.000153 * T ** 2) % 360

    C = (1.914602 - 0.004817 * T - 0.000014 * T ** 2) * sin(radians(M)) + (0.019993 - 0.000101 * T) * sin(2 * radians(M)) + 0.000289 * sin(3 * radians(M))

    # % 好像是地球日心黄经
    O_ = LL + C
    Womiuga = 125.04 - 1934.136 * T
    # % 需要的两个参数
    Lamda = O_ - 0.00569 - 0.00478 * sin(radians(Womiuga))
    Lamda_ = O_ - 1.397 * T - 0.00031 * T ** 2
    # % Astronomical Algorithms中公式22.2
    epxlong0 = 23 + 26 / 60 + 21.448 / 3600 - 46.8150 / 3600 * T - 0.00059 / 3600 * T ** 2 + 0.001813 / 3600 * T ** 3

    L_ = 280.4665 + 36000.27698 * T
    L__ = 218.3165 + 481267.28813 * T
    Q = (125.04452 - 1934.136261 * T + 0.0020708 * T ** 2 + T ** 3 / 450000) % 360
    # % Astronomical Algorithms P144
    depxlong = 9.20 / 3600 * cos(Q * pi / 180) + 0.57 / 3600 * cos(2 * radians(L_)) + 0.1 / 3600 * cos( 2 * radians(L__)) - 0.09 / 3600 * cos(2 * radians(Q))
    # % 需要的参数2，黄道的倾斜
    epxlong = epxlong0 + depxlong

    sita = ((JD - 2398220) * 360 / 25.38) % 360     #% 单位：度

    I = 7.25       #% 单位：度
    K = 73.6667 + 1.3958333 * (JD - 2396758) / 36525   #% 单位：度

    x = degrees(atan(-cos(radians(Lamda_)) * tan(radians(epxlong))))     # % 单位：度
    y = degrees(atan(-cos(radians(Lamda - K)) * tan(radians(I))))     # % 单位：度

    P = x + y      #% 单位：度
    # % P: 从太阳圆盘北点向东测量的旋转轴北端的位置角;
    # % 当太阳旋转轴的北极向东倾斜时P为正，如果向西倾斜则为负。
    # % 天体和太阳北极可相差26度.
    # % P在4月7日左右达到最小值 - 26.3度，在10月11日左右达到最大值为 + 26.3度，
    # % ，在1月5日和7月7日附近为零。

    B0 = degrees(asin(sin(radians(Lamda - K)) * sin(radians(I))))        # % 单位：度
    # % % B0: the heliographic latitude of the center of the solar disk
    # % B0代表太阳的北极向地球的（+）或远离（ - ）倾斜。
    # % 它在6月6日和12月7日为零，
    # % 并在3月6日（-7.25度）和9月8日（+7.25度）。

    eina = degrees(atan(tan(radians(Lamda - K)) * cos(radians(I))))       # % 单位：度
    # % eina应与Lamda - K - 180
    #     # 处于同一象限 ？？？？？？？？？
    #     # % LK(i, 1) = Lamda - K - 180;
    LK = Lamda - K - 180
    eina_1 = getQuadrantby_eina(eina)[0]
    LK_1 = getQuadrantby_eina(LK)[0]
    if eina_1 == LK_1:
        eina = eina - 180
    # % 日面经度L0每天减少约13.2度。
    # % 平均朔望期为27.2752天。
    # % 每个“同步旋转”的开始是L0经过0°的瞬间。 1号旋转于1853年11月9日开始。
    L0 = (eina - sita - 180) % 360      #% 单位：度
    return [B0, L0, P]


def getQuadrantby_eina(eina):
    # % 根据度数判断所在象限
    # % 1, 2, 3, 4分别代表一，二，三，四象限
    eina_ = eina % 360
    eina_1 = []
    if eina_ > 0 and eina_ <= 90:
        eina_1.append(1)
    if eina_ > 90 and eina_ <= 180:
        eina_1.append(2)
    if eina_ > 180 and eina_ <= 270:
        eina_1.append(3)
    if eina_ > 270 and eina_ <= 360:
        eina_1.append(4)
    return eina_1


def GetTbyDate(year, month, day):
    """% 计算儒略日
% 输入参数：公历年，月，日(日期可以带小数点)
% 输出参数：以儒略日表示的观测时间
% 公式参考：http://www.360doc.com/content/17/0814/11/460866_679072040.shtml

% 测试代码
% 1858年11月17日世界时0时的儒略日数为：2400000.5
% 2000年1月1日世界时12时的儒略日数为：2451545.0
% 1992年10月13日零时的儒略日数为：2448908.5

"""
    if month <= 2:
        year = year - 1
        month = month + 12
    A = floor(year / 100)
    B = 2 - A + floor(A / 4)
    T = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + B - 1524.5
    return T

def XY_to_NE(imageName, Xa, Ya):
    # 提取时间
    Utime = imageName[11:26]
    # 转为时间数组
    timeArray = time.strptime(Utime, "%Y%m%d_%H%M%S")
    [B0, L0, P] = get_b0l0p_by_data(timeArray.tm_year,
                                    timeArray.tm_mon,
                                    timeArray.tm_mday,
                                    timeArray.tm_hour,
                                    timeArray.tm_min);
    appdiam = 32 / 60 + 35 / 3600
    #radii = 150/2  # 4096*4096  图像半径
    radii = 1900 * 0.5
    P = 0
    # Xa = Xa * 150 / 4096
    # Ya = Ya * 150 / 4096

    [Bp, Lp] = get_bl_by_xyblr(Xa, Ya, P, B0, L0, radii, appdiam)
    return [Bp, Lp, L0]


if __name__ == '__main__':
    #[B0, L0, P] = get_b0l0p_by_data(1999, 1, 1, 11, 10)
    [B0, L0, P] = get_b0l0p_by_data(2012, 11, 18, 4, 48)
    appdiam = 32 / 60 + 35 / 3600
    Xa = (312.5 - 2039) // 2
    Ya = (2051 - 2843.5) // 2
    radii = 1900 * 0.5
    P = 0
    [Bp, Lp] = get_bl_by_xyblr(Xa, Ya, P, B0, L0, radii, appdiam)
    
    # import  math
    # 通过投影校正后计算面积不准确的问题未解决
    # R = 6960
    # adjustArea = abs(36266 * 1e5 / (8 * math.pi * R * R * math.cos(math.radians(Bp)) * math.cos(math.radians(Lp-L0))))