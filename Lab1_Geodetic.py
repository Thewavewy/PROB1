import pandas as pd
import pyproj
import pygeodesy
from geopy.distance import geodesic
import pygeodesy.dms as dms
from pygeodesy.ellipsoidalKarney import LatLon
import math
from math import pi, sqrt, pow, atan, sin, cos, radians


# แปลงองศา เป็น องศา ลิบดา ฟิลิบดา
def dd2dms(dd, PREC=5):
    # arcsecond is 5-digit decimal fraction
    return dms.toDMS(dd, prec=5)


def dms2dd(dms_str, SEP=' '):
    return dms.parseDMS(dms_str, sep=SEP)


# คำนวณรัศมีโลก ณ Latitude ที่กำหนด
def radius(B):
    # The geocentric radius is the distance from the Earth's center to a point on the spheroid surface at the geodetic latitude.
    # Ref: https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius
    # input B is Latitude in Degree unit
    # output radius in Km unit
    B = radians(B)  # converting into radians
    a = 6378.137  # Radius at sea level at equator
    b = 6356.752  # Radius at poles
    c = (a ** 2 * cos(B)) ** 2
    d = (b ** 2 * sin(B)) ** 2
    e = (a * cos(B)) ** 2
    f = (b * sin(B)) ** 2
    R = sqrt((c + d) / (e + f))

    return R


# แปลงค่าพิกัด xyz เป็น ละติจูด ลองจิจูด และ ความสูงเหนือทรงรี  หน่วยดีกรี
def xyz2lla(x, y, z):
    ecef = '+proj=geocent +ellps=WGS84'  # Cartisian
    lla = "epsg:4326"  # WGS84 Geodetic
    transproj = pyproj.Transformer.from_crs(ecef, lla)
    lat, lon, ell_h = transproj.transform(x, y, z, radians=False)
    return lat, lon, ell_h


def lla2utm(lat, lon, ell_h, zone):
    '''
        I = lat ,long unit Degree
        I = ell_h unit metre
        O = E,N,h unit metre
    '''
    lla = "epsg:4326"  # WGS84 Geodetic
    utm = "epsg:326" + str(zone)  # WGS84 UTM Zone 47N
    transproj = pyproj.Transformer.from_crs(lla, utm, always_xy=True)  # always_xy หมายถึงให้เรียง lon,lat
    E, N, h = transproj.transform(lon, lat, ell_h, radians=False)
    return E, N, h


# คำนวนสเกลแฟกเตอร์ในระบบพิกัด UTM โดยระบบพิกัดละติจูด ลองจิจูด และ โซน
def utm_scale_factor(lat, lon, zone):
    '''
        I = lat ,long unit Degree
        I = zone (UTM)
        O = scaleX, scaleY, scaleArea, mer_conv
    '''
    prj = pyproj.Proj('epsg:326' + zone)
    fc = prj.get_factors(lon, lat, radians=False)
    scaleX = fc.parallel_scale
    scaleY = fc.meridional_scale
    scaleArea = fc.areal_scale
    mer_conv = fc.meridian_convergence
    return scaleX, scaleY, scaleArea, mer_conv


# คำนวนความสูงเหนือระดับน้ำทะเลปานกลาง อ้างอิงแบบจำลอง TGM2017
def MSL(lat, lon, h):
    #  Download: tgm2017-1_2.pgm size 469 MB
    #  https://www.priabroy.name/archives/sdm_downloads/tgm2017-1
    geoid_model = r'C:\Users\USER\Desktop\StudySection\Semester 4.1\2108415 Geodetic\tgm2017-1_2.pgm'
    geoid_interpolator = pygeodesy.GeoidKarney(geoid_model)
    # Get the geoid height
    pos = LatLon(lat, lon)
    N = geoid_interpolator(pos)
    print(N)
    return h - N


def xyz2tm(x, y, z):
    tm = '+proj=tmerc +lat_0=0.0 +lon_0=0.0 +k=1.00000000 +x_0=500000 +y_0=0 +ellps=WGS84 +units=m +no_defs'
    ecef = {"proj": 'geocent', "ellps": 'WGS84'}
    transproj = pyproj.Transformer.from_crs(ecef, tm)
    E, N, h = transproj.transform(x, y, z)
    return (E, N, h)


# แปลงค่าพิกัด xyz เป็นพิกัด Lambert Conformal Conic (LCC)
def xyz2LCC(x, y, z):
    # Ref: https://proj.org/operations/projections/index.html
    # Ref: https://proj.org/operations/projections/lcc.html

    # Lambert Conformal Conic 1 standard parallel Conversion
    llc_1sp = '+proj=lcc +lon_0=99 +lat_0=0 +x_0=500000 +y_0=0 +ellps=WGS84 +lat_1=33'
    # Lambert Conformal Conic conversion with scale factor
    llc_2sp_1 = '+proj=lcc +lon_0=99 +lat_0=0 +x_0=500000 +y_0=0 +ellps=WGS84 +k=0.9996'
    # Lambert Conformal Conic 2 standard parallel Conversion
    llc_2sp_2 = '+proj=lcc +lon_0=99 +lat_0=0 +x_0=500000 +y_0=0 +ellps=WGS84 +lat_1=9 +lat_2=12'

    ecef = {"proj": 'geocent', "ellps": 'WGS84'}
    transproj = pyproj.Transformer.from_crs(ecef, llc_1sp)
    E, N, h = transproj.transform(x, y, z)
    return (E, N, h)




# --------------------------------------------------------------------------------
# 1) แปลงค่าพิกัด ECEF เหล่านี้เป็นค่าพิกัด Geodetic และ UTM Zone 47, 48 และความสูงเหนือ MSL ด้วย TGM2017
# สร้าง dataframe เก็บข้อมูลพิกัดในข้อที่ 1
'''
data1 = {
    'Name': ['AKSN', 'AMKO', 'APKN'],
    'x': [-1482251.762, -883143.728, -1037944.776],
    'y': [5925272.273, 6010905.970, 6165673.259],
    'z': [1831475.943, 1937627.074, 1255828.562]
}
df1 = pd.DataFrame(data1)

# สร้าง loop ที่คำนวณทุกอย่างแล้วแทรกค่าใน dataframe เดิม

Lat = list()
Lon = list()
Hi = list()
e47 = list()
n47 = list()
hi47 =list()
e48 = list()
n48 = list()
hi48 =list()
mesl = list()
for i in range(df1.shape[0]): #rows = 3
    (x,y,z) = df1['x'][i], df1['y'][i], df1['z'][i]
    La,Lo,H = xyz2lla(x, y, z)
    E47, N47, h47 = lla2utm(La, Lo, H, 47)
    E48, N48, h48 = lla2utm(La, Lo, H, 48)
    msl = MSL(La, Lo, H)
    # แทรกค่าเข้าไปใน list + formatting
    Lat.append(dd2dms(La, PREC=5))
    Lon.append(dd2dms(Lo, PREC=5))
    Hi.append(round(H,3))
    e47.append(round(E47,3))
    n47.append(round(N47,3))
    hi47.append(round(h47,3))
    e48.append(round(E48,3))
    n48.append(round(N48,3))
    hi48.append(round(h48,3))
    mesl.append(round(msl,3))


print(mesl)
# สร้าง column ใหม่ แสดงผลลัพธ์ใน dataframe

df1['Latitude'] = Lat
df1['Longitude'] = Lon
df1['H-Geodetic'] = Hi
df1['E-47'] = e47
df1['N-47'] = n47
df1['h-47'] = hi47
df1['E-48'] = e48
df1['N-48'] = n48
df1['h-48'] = hi48
df1['MSL'] = mesl

# ส่งออกผลลัพธ์เป็น excel
# Export the DataFrame to an Excel file
excel_filename = 'result1.xlsx'
df1.to_excel(excel_filename, index=False)

print(f"DataFrame exported to {excel_filename}")
'''
# 2) พัฒนาฟังก์ชันสำหรับแปลงค่าพิกัด ปรับปรุงจาก code เดิม

def utm2lla(E,N,h,zone):
    P = pyproj.Proj(proj = 'utm',zone = zone, ellps = 'WGS84')
    La, Lo, hi = P(E,N,inverse = True)
    return E, N, hi

def lla2xyz(La,Lo,Hi):
    # ความสัมพันธ์ของ Geodetic และ ECEF
    # กำหนดค่าคงที่และค่าที่ใช้ในสมการ
    a = 6317137 ;f = 1/298.257223563 ;e = sqrt((2*f)-(f**2))
    Phi = radians(La) ;Lam = radians(Lo) ;H = Hi

    # คำนวณ X,Y,Z
    N = (a)/sqrt(1-((e**2)*(sin(Phi)**2)))
    X = (N+H)*(cos(Phi))*(cos(Lam))
    Y = (N+H)*(cos(Phi))*(sin(Lam))
    Z = ((1-e**2)*(N+H))*(sin(Phi))
    return X,Y,Z

def utm2xyz(E,N,H,zone):
    La,Lo,Hi = utm2lla(E,N,H,zone)
    X,Y,Z = lla2xyz(La,Lo,Hi)
    return X,Y,Z

def xyz2utm(X,Y,Z,zone):
    La,Lo,H = xyz2lla(X,Y,Z)
    E,N,Hi = lla2utm(La,Lo,H,zone)
    return E,N,Hi


# 3) แปลงค่าพิกัด ECEF เหล่านี้เป็นค่าพิกัด Geodetic และ UTM Zone 47
# แต่นำข้อมูลพิกัดจาก CSV File มาแปลง

path = r"C:\Users\USER\Desktop\StudySection\Semester 4.1\2108415 Geodetic\Lab01_MapProjection-513652-16927710625159\Lab01_MapProjection\coordinate.csv"

df3 = pd.read_csv(path)

# เปลี่ยน z+space ให้เป็น z
df3.rename(columns = {'z ':'z'}, inplace = True)
print(df3)
# สร้าง loop ที่คำนวณทุกอย่างแล้วแทรกค่าใน dataframe เดิม

Lat = list()
Lon = list()
Hi = list()
e47 = list()
n47 = list()
hi47 =list()

for i in range(df3.shape[0]): #rows = 134
    (x,y,z) = df3['x'][i], df3['y'][i], df3['z'][i]
    La,Lo,H = xyz2lla(x, y, z)
    E47, N47, h47 = xyz2utm(x,y,z,47) #ใช้ฟังก์ชันที่พัฒนาแล้ว
    # แทรกค่าเข้าไปใน list + formatting
    Lat.append(dd2dms(La, PREC=5))
    Lon.append(dd2dms(Lo, PREC=5))
    Hi.append(round(H,3))
    e47.append(round(E47,3))
    n47.append(round(N47,3))
    hi47.append(round(h47,3))
# สร้าง column ใหม่ แสดงผลลัพธ์ใน dataframe

df3['Latitude'] = Lat
df3['Longitude'] = Lon
df3['H-Geodetic'] = Hi
df3['E-47'] = e47
df3['N-47'] = n47
df3['h-47'] = hi47

# ส่งออกผลลัพธ์เป็น excel
# Export the DataFrame to an Excel file
excel_filename = 'result3.xlsx'
df3.to_excel(excel_filename, index=False)
print(f"DataFrame exported to {excel_filename}")

'''
# 4)คำนวณระยะทางราบบนจุด A-B โดยใช้ Scale Factor
# zone 47N
# ระยะทางกริด UTM คำนวณจาก xyz2utm
A = [-1482251.762, 5926272.273, 1831475.943]
B = [-1037944.776, 6165673.259, 1255828.562]
(xA, yA, zA) = A[0], A[1], A[2]
(xB, yB, zB) = B[0], B[1], B[2]
eA,nA,hA = xyz2utm(xA,yA,zA,47)
eB,nB,hB = xyz2utm(xB,yB,zB,47)
Grid_47 = sqrt(((eA-eB)**2)+((nA-nB)**2))
# ระยะบนทรงรีคำนวณจาก Geodesic Distance
LA,LOA,HA = xyz2lla(xA,yA,zA)
LB,LOB,HB = xyz2lla(xB,yB,zB)
Alist = [LA,LOA]
Blist = [LB,LOB]
geo_dis_47 = geodesic(Alist, Blist, ellipsoid='WGS-84').m
# ระยะทางจริงบนพื้นโลก
h = (HA+HB)/2 ; R = 6378137
ESF = (R)/(R+h)
real_dis_47 = geo_dis_47/ESF
print(real_dis_47)

#zone 48N
eA,nA,hA = xyz2utm(xA,yA,zA,48)
eB,nB,hB = xyz2utm(xB,yB,zB,48)
Grid_48 = sqrt(((eA-eB)**2)+((nA-nB)**2))
# ระยะบนทรงรีคำนวณจาก Geodesic Distance
LA,LOA,HA = xyz2lla(xA,yA,zA)
LB,LOB,HB = xyz2lla(xB,yB,zB)
Alist = [LA,LOA]
Blist = [LB,LOB]
geo_dis_48 = geodesic(Alist, Blist, ellipsoid='WGS-84').m
# ระยะทางจริงบนพื้นโลก
h = (HA+HB)/2 ; R = 6378137
ESF = (R)/(R+h)
real_dis_48 = geo_dis_48/ESF
print(real_dis_48)

# 5)ออกแบบการฉายแผนที่ TM
# เฉลี่ย scale factor ใน zone ให้เหมาะสม
lat_min = float(input('ใส่ lat_min:'))
lat_max = float(input('ใส่ lat_max:'))
lon_min = float(input('ใส่ lon_min:'))
lon_max = float(input('ใส่ lon_max:'))
zone = int(input('ใส่เลขโซน:'))
n = int(input('จำนวนสุ่ม:'))

import math
import random

def utm_sf(lat,lon,k0,zone):
    a = 6378137 #m in WGS84
    b = 6356752.3142 #m in WGS84
    cen = ((zone-30)*6)-3
    lon = cen - lon
    e2 = (a**2 - b**2) / b**2
    n = e2 * (math.cos(math.radians(lat))**2)
    F2 = (1 + (n**2)) / 2
    L2 = (math.radians(lon) * math.cos(math.radians(lat))) ** 2
    t = math.tan(math.radians(lat))
    F4 = (5 - (4 * (t**2)) + (n**2) * (9 - 24 * (t**2))) / 12
    k = k0 * (1 + (F2 * L2 * (1 + F4 * L2)))
    return k

# random n จากขอบเขต
sum_sf = []
for j in range(n):
    lat_v =list()
    lon_v =list()
    sum_r = 0
    for i in range(1,n+1):
        lat_r = float(round(lat_min + (lat_max - lat_min)*random.uniform(0,1),3))
        lon_r = float(round(lon_min + (lon_max - lon_min)*random.uniform(0,1),3))
        lat_v.append(lat_r)
        lon_v.append(lon_r)
        i += 1
    for i in range(len(lat_v)):
        sf_r = utm_sf(lat_v[i],lon_v[i],0.9996,zone)
        se_r = 1 - sf_r
        sum_r = sum_r+se_r
    mean_sf_r = round(abs(1-(sum_r/n)),8)
    sum_sf.append(mean_sf_r)
result = round(sum(sum_sf)/n,8)
print('ค่า k ที่ต้องใส่ไปใน QGIS ='+' '+str(result)+'\n')
print('pyproj code')
print('+proj=tmerc +lat_0=0.0 +lon_0='+str(((zone-30)*6)-3)+' k0='+str(result)+' x_0=500000 y_0=0 +a=6378137.0 +b=6356752.3142 +units=m +no_defs')
'''