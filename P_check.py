from math import *


def get_coordinate(xa, ya, P, B0, L0, radius, appdiam):
    # ------------------------------------------------------------------
    # --功能：从日面图坐标换算到太阳经纬度
    # radius：太阳圆半径；(Xa,Ya)：以图片中心为原点的坐标；
    # B0（L0）：圆盘中心的日面纬度（经度）；P:从太阳圆盘北点向东测量的旋转轴北端的位置角;
    # appdiam：太阳视直径，约在31.50'–32.60'之间变化，取平均值32.05'
    # 输出坐标(Xa,Ya)对应的经度lon,纬度lat;
    # lon：经度；lat：纬度。
    # ------------------------------------------------------------------
    # sita = degrees(atan(xa/ya)) - 180
    diskdiam = 2 * radius
    # 单位均为：度
    rou1 = (appdiam / diskdiam) * sqrt(xa * xa + ya * ya)
    ra = 2 * rou1 / appdiam
    if ra > 1:
        ra = 1
    elif ra < -1:
        ra = -1

    rou = degrees(asin(ra)) - rou1
    if ya < 0:
        sita = degrees(atan(xa / ya)) - 180
    else:
        sita = degrees(atan(xa / ya))
    B0_ = radians(B0)
    rou_ = radians(rou)
    P_ = radians(P)
    site_ = radians(sita)
    Bp = degrees(asin(sin(B0_) * cos(rou_) + cos(B0_) * sin(rou_) * cos(P_ - site_)))
    Bp_ = radians(Bp)
    Lp = (degrees(asin(sin(rou_) * sin(P_ - site_) / cos(Bp_))) + L0) % 360
    Lp = L0-Lp
    return Bp, Lp


def get_blp(year, month, day, hour, minute):
    # -------------------------------------------------------------------------------
    # 输入：年，月，日，时，分 (世界标准时间，UT)
    # 参考文献：Astronomical Algorithms 第22、25和29章
    # -------------------------------------------------------------------------------
    day = day + hour / 24 + minute / 24 / 60
    JD = get_JulianDay(year, month, day)  # 计算儒略日数
    T = (JD - 2451545.0) / 36525  # P151
    LL = (280.46646 + 36000.76983 * T + 0.0003032 * T ** 2) % 360  # 角度，P171
    M = (357.52911 + 35999.05029 * T - 0.000153 * T ** 2) % 360  # 角度；the mean anomaly(近点角) of the sun

    C = (1.914602 - 0.004817 * T - 0.000014 * T ** 2) * sin(radians(M)) \
        + (0.019993 - 0.000101 * T) * sin(2 * radians(M)) + 0.000289 * sin(
        3 * radians(M))  # sun's equation of the center

    O_ = LL + C  # sun's true longitude, P172
    # omega = 125.04 - 1934.136 * T
    omega = 125.04452 - 1934.136261 * T + 0.0020708 * T ** 2 + T ** 3 / 450000  # P152, chap22

    Lamda = O_ - 0.00569 - 0.00478 * sin(radians(omega))
    Lamda_ = O_ - 1.397 * T - 0.00031 * T ** 2  # 单位：度，P174
    epsilon0 = (84381.448 - 46.815 * T - 0.00059 * T ** 2 + 0.001813 * T ** 3) / 3600  # 单位：度

    L = 280.4665 + 36000.27698 * T
    L_ = 218.3165 + 481267.28813 * T
    # = (125.04452 - 1934.136261 * T + 0.0020708 * T ** 2 + T ** 3 / 450000) % 360
    # % Astronomical Algorithms P144
    delta_epsilon = (9.20 * cos(omega * pi / 180) +
                     0.57 * cos(2 * radians(L)) +
                     0.1 * cos(2 * radians(L_)) -
                     0.09 * cos(2 * radians(omega))) / 3600  # P152
    epsilon = epsilon0 + delta_epsilon

    sita = ((JD - 2398220) * 360 / 25.38) % 360  # 单位：度
    I = 7.25  # 单位：度
    K = 73.6667 + 1.3958333 * (JD - 2396758) / 36525  # 单位：度

    x = degrees(atan(-cos(radians(Lamda_)) * tan(radians(epsilon))))  # 单位：度
    y = degrees(atan(-cos(radians(Lamda - K)) * tan(radians(I))))  # 单位：度

    P = x + y  # 单位：度
    # 当太阳旋转轴的北极向东倾斜时P为正，如果向西倾斜则为负。
    # 天体和太阳北极可相差26度，P在4月7日左右达到最小值 - 26.3度，在10月11日左右达到最大值为 + 26.3度，在1月5日和7月7日附近为零。

    B0 = degrees(asin(sin(radians(Lamda - K)) * sin(radians(I))))  # 单位：度
    # B0代表太阳的北极向地球的（+）或远离（ - ）倾斜。
    # 它在6月6日和12月7日为零，并在3月6日（-7.25度）和9月8日（+7.25度）。

    eta = degrees(atan(tan(radians(Lamda - K)) * cos(radians(I))))  # 单位：度

    LK = Lamda - K - 180
    eta_1 = get_quadrant(eta)[0]
    LK_1 = get_quadrant(LK)[0]
    if eta_1 == LK_1:
        eta = eta - 180
    # 日面经度L0每天减少约13.2度，平均朔望期为27.2752天。
    # 每个“同步旋转”的开始是L0经过0°的瞬间，1号旋转于1853年11月9日开始。
    L0 = (eta - sita - 180) % 360  # 单位：度
    return B0, L0, P


def get_quadrant(eta):
    # % 根据度数判断所在象限
    # % 1, 2, 3, 4分别代表一，二，三，四象限
    eta_ = eta % 360
    eta_1 = []
    if 0 < eta_ <= 90:
        eta_1.append(1)
    if 90 < eta_ <= 180:
        eta_1.append(2)
    if 180 < eta_ <= 270:
        eta_1.append(3)
    if 270 < eta_ <= 360:
        eta_1.append(4)
    return eta_1


def get_JulianDay(year, month, day):
    # --------------------------------------------
    # --计算儒略日
    # 输入参数：公历年，月，日(日期可以带小数点)；输出参数：以儒略日表示的观测时间
    # 公式参考：http://www.360doc.com/content/17/0814/11/460866_679072040.shtml
    # ---------------------------------------------
    """
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
    JD = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + B - 1524.5
    return JD


def check_xy(x, y, P):

    if P >= 0:
        beta = radians(P)
        X1 = x*cos(beta) + y*sin(beta)
        Y1 = y*cos(beta) - x*sin(beta)
    else:
        Z = abs(P)
        beta = radians(Z)
        X1 = (x / sin(beta) - y / cos(beta)) / (tan(beta) + 1 / tan(beta))
        Y1 = (y / sin(beta) + x / cos(beta)) / (tan(beta) + 1 / tan(beta))
    return X1, Y1

if __name__ == '__main__':
    print('算例结果：')
    # B0, L0, P = get_blp(1999, 1, 1, 11, 10)
    # B0, L0, P = get_blp(2004, 6, 1, 12, 0)
    B0, L0, P = get_blp(2022, 3, 1, 19, 0)
    appdiam = 32 / 60 + 20 / 3600
    Xa, Ya = check_xy(52.3, 28.4, P)
    radii = 150 / 2
    # radii = 950
    Bp, Lp = get_coordinate(Xa, Ya, P, B0, L0, radii, appdiam)
    print("纬度：%f" % Bp)
    print("经度：%f" % Lp)


