from math import *
from numpy import array
from shapely.geometry import Point
import csv
import geopandas as gpd
from matplotlib.colors import ListedColormap

# Define the input points
phi1_WGS = 50.0780391
la1_WGS = 14.4488108

phi2_WGS = 50.0785864
la2_WGS = 14.4509351

def uv_to_sd(u, v, uk, vk):
    """Function that converts latitude and longitude"""
    dv = vk - v
    s = asin(sin(u)*sin(uk) + cos(u)*cos(uk)*cos(dv))
    d = atan2(cos(u)*sin(dv), cos(u)*sin(uk)*cos(dv) - sin(u)*cos(uk))
    return array([s, d])

def WGSToJTSK(phi_WGS, la_WGS):
    """Function to convert WGS84 coordinates to JTSK coordinates"""
    
    #Radian conversion
    phir_WGS = phi_WGS * pi / 180
    lar_WGS = la_WGS * pi / 180


    #WGS84, parameters
    a_WGS = 6378137.0000
    b_WGS = 6356752.3142
    
    e2_WGS = (a_WGS * a_WGS - b_WGS * b_WGS) / (a_WGS * a_WGS)
    W_WGS = sqrt(1 - e2_WGS * (sin(phir_WGS))**2)
    N_WGS = a_WGS / W_WGS
    
    #Geocentric coordinates, WGS84
    X_WGS = N_WGS * cos(phir_WGS) * cos(lar_WGS)
    Y_WGS = N_WGS * cos(phir_WGS) * sin(lar_WGS)
    Z_WGS = N_WGS * (1 - e2_WGS) * sin(phir_WGS)
    
    #7 parameters of Helmert transformation
    om_x = 4.9984 / 3600 * pi / 180
    om_y = 1.5867 / 3600 * pi / 180
    om_z = 5.2611 / 3600 * pi / 180
    
    m = 1 - 3.5623e-6
    
    dx = -570.8285
    dy = -85.6769
    dz = -462.8420
    
    #Performing Helmert transformation
    R = array([[1, om_z, -om_y], [-om_z, 1, om_x], [om_y, -om_x, 1]])
    
    delta = array([dx, dy, dz])
    xyz_wgs = array([X_WGS, Y_WGS, Z_WGS])
    xyz_bess = m * R.dot(xyz_wgs) + delta
    
    #Geocentrinc coordinates, Bessel
    X_B, Y_B, Z_B = xyz_bess

    #Bessel, parameters
    a_B = 6377397.1550
    b_B = 6356078.9633
    e2_B = (a_B * a_B - b_B * b_B) / (a_B * a_B)
    
    #Converting geocentric coordinates to geodetic coordinates
    lar_B = atan2(Y_B, X_B)
    phir_B = atan2(Z_B, (1 - e2_B) * sqrt(X_B**2 + Y_B**2))
    
    #Reduction to Ferro meridian
    larf_B = lar_B + (17 + 2 / 3) * pi / 180

    #Gaussian conformal projection, parameters
    phir_0 = 49.5 * pi / 180
    alpha = sqrt(1 + (e2_B * (cos(phir_0))**4) / (1 - e2_B))
    u0_r = asin(sin(phir_0) / alpha)
    kn = (tan(phir_0 / 2 + pi / 4))**alpha * ((1 - sqrt(e2_B) * sin(phir_0)) / (1 + sqrt(e2_B) * sin(phir_0)))**(alpha * sqrt(e2_B) / 2)
    kd = tan(u0_r / 2 + pi / 4)
    k = kn / kd
    R = (a_B * sqrt(1 - e2_B)) / (1 - e2_B * (sin(phir_0))**2)
    
    #Gaussian conformal projection
    arg = 1 / k * (tan(phir_B / 2 + pi / 4) * ((1 - sqrt(e2_B) * sin(phir_B)) / (1 + sqrt(e2_B) * sin(phir_B)))**(sqrt(e2_B) / 2))**alpha
    u = 2 * (atan(arg) - pi / 4)
    v = alpha * larf_B
    
    #Transformation to oblíque aspect
    uk = (59 + 42 / 60 + 42.6969 / 3600) * pi / 180
    vk = (42 + 31 / 60 + 31.41725 / 3600) * pi / 180
    s, d = uv_to_sd(u, v, uk, vk)
    
    #LCC
    s0r = 78.5 * pi / 180
    c = sin(s0r)
    rho0 = 0.9999 * R * (1/tan((s0r)))
    
    rho = rho0 * ((tan(s0r / 2 + pi / 4)) / (tan(s / 2 + pi / 4)))**c
    eps = c * d
    
    #Conversion to Cartesian coordinates
    x = rho * cos(eps)
    y = rho * sin(eps)
    
    #Local linear scale
    mr1 = (c * rho) / (R * cos(s))
    delta_rho = (rho - rho0) / 100000
    mr2 = 0.9999 + 0.00012282 * delta_rho**2 - 0.00000315 * delta_rho**3 + 0.00000018 * delta_rho**4
    
    #Convergence
    ksi = asin(cos(uk) * sin(pi - d) / cos(u))
    conv1 = (eps - ksi) * 180 / pi
    conv2 = 0.008257 * y / 1000 + 2.373 * y / x
    return round(x, 2), round(y, 2), round(mr1, 6), round(conv1, 6)

# Convert input points from WGS-84 to JTSK
x1, y1, mr1, conv1 = WGSToJTSK(phi1_WGS, la1_WGS)
x2, y2, mr2, conv2 = WGSToJTSK(phi2_WGS, la2_WGS)

# print the results
print("WGS-84 to JTSK")
print("Point 1: ", x1, y1, mr1, conv1)
print("Point 2: ", x2, y2, mr2, conv2)

d1 = sqrt((x1 - x2)**2 + (y1 - y2)**2)
print("Distance between points: ", round((d1),2))
#round to 5 dihgits


def BessToJTSK(phi_B, la_B):
    """Function to convert Bessel coordinates to JTSK coordinates"""
    
    # Radian conversion
    phir_B = phi_B * pi / 180
    lar_B = la_B * pi / 180
    
    # Bessel, parameters
    a_B = 6377397.1550
    b_B = 6356078.9633
    e2_B = (a_B * a_B - b_B * b_B) / (a_B * a_B)

    #Reduction to Ferro meridian
    larf_B = lar_B + (17 + 2 / 3) * pi / 180

    #Gaussian conformal projection, parameters
    phir_0 = 49.5 * pi / 180
    alpha = sqrt(1 + (e2_B * (cos(phir_0))**4) / (1 - e2_B))
    u0_r = asin(sin(phir_0) / alpha)
    kn = (tan(phir_0 / 2 + pi / 4))**alpha * ((1 - sqrt(e2_B) * sin(phir_0)) / (1 + sqrt(e2_B) * sin(phir_0)))**(alpha * sqrt(e2_B) / 2)
    kd = tan(u0_r / 2 + pi / 4)
    k = kn / kd
    R = (a_B * sqrt(1 - e2_B)) / (1 - e2_B * (sin(phir_0))**2)
    
    #Gaussian conformal projection
    arg = 1 / k * (tan(phir_B / 2 + pi / 4) * ((1 - sqrt(e2_B) * sin(phir_B)) / (1 + sqrt(e2_B) * sin(phir_B)))**(sqrt(e2_B) / 2))**alpha
    u = 2 * (atan(arg) - pi / 4)
    v = alpha * larf_B
    
    #Transformation to oblíque aspect
    uk = (59 + 42 / 60 + 42.6969 / 3600) * pi / 180
    vk = (42 + 31 / 60 + 31.41725 / 3600) * pi / 180
    s, d = uv_to_sd(u, v, uk, vk)
    
    #LCC
    s0r = 78.5 * pi / 180
    c = sin(s0r)
    rho0 = 0.9999 * R * (1/tan((s0r)))
    rho = rho0 * ((tan(s0r / 2 + pi / 4)) / (tan(s / 2 + pi / 4)))**c
    eps = c * d
    
    #Conversion to Cartesian coordinates
    x = rho * cos(eps)
    y = rho * sin(eps)
    
    #Local linear scale
    mr1 = (c * rho) / (R * cos(s))
    delta_rho = (rho - rho0) / 100000
    mr2 = 0.9999 + 0.00012282 * delta_rho**2 - 0.00000315 * delta_rho**3 + 0.00000018 * delta_rho**4
    
    #Convergence
    ksi = asin(cos(uk) * sin(pi - d) / cos(u))
    conv1 = (eps - ksi) * 180 / pi
    conv2 = 0.008257 * y / 1000 + 2.373 * y / x
    return round(x, 2), round(y, 2), round(mr1, 6), round(conv1, 6)

# Convert input points from Bessel to JTSK
x3, y3, mr3, conv3 = BessToJTSK(phi1_WGS, la1_WGS)
x4, y4, mr4, conv4 = BessToJTSK(phi2_WGS, la2_WGS)

# print the results
print("Bessel to JTSK")
print("Point 1: ", x3, y3, mr3, conv3)
print("Point 2: ", x4, y4, mr4, conv4)

# Distance between points, Bessel to JTSK
d2 = sqrt((x4 - x3)**2 + (y4 - y3)**2)
print("Distance between points: ", round((d2),2))

def SphereToJTSK(u, v):
    """Fucntion to convert spherical coordinates to JTSK coordinates"""
    
    #Radian conversion
    u = u * pi / 180
    v = v * pi / 180
    
    #Reduction to Ferro meridian
    vf = v + (17 + 2 / 3) * pi / 180
    
    #Transformation to oblíque aspect
    uk = (59 + 42 / 60 + 42.6969 / 3600) * pi / 180
    vk = (42 + 31 / 60 + 31.41725 / 3600) * pi / 180
    s, d = uv_to_sd(u, vf, uk, vk)
    
    #LCC
    R = 6380703.6104635
    s0r = 78.5 * pi / 180
    c = sin(s0r)
    rho0 = 0.9999 * R * (1/tan((s0r)))
    
    rho = rho0 * ((tan(s0r / 2 + pi / 4)) / (tan(s / 2 + pi / 4)))**c
    eps = c * d
    
    #Conversion to Cartesian coordinates
    x = rho * cos(eps)
    y = rho * sin(eps)
    
    #Local linear scale
    mr1 = (c * rho) / (R * cos(s))
    delta_rho = (rho - rho0) / 100000
    mr2 = 0.9999 + 0.00012282 * delta_rho**2 - 0.00000315 * delta_rho**3 + 0.00000018 * delta_rho**4
    
    #Convergence
    ksi = asin(cos(uk) * sin(pi - d) / cos(u))
    conv1 = (eps - ksi) * 180 / pi
    conv2 = 0.008257 * y / 1000 + 2.373 * y / x
    return round(x, 2), round(y, 2), round(mr1, 6), round(conv1, 6)

# Convert input points from Sphere to JTSK
x5, y5, mr5, conv5 = SphereToJTSK(phi1_WGS, la1_WGS)
x6, y6, mr6, conv6 = SphereToJTSK(phi2_WGS, la2_WGS)

print("Sphere to JTSK")
print("Point 1: ", round(x5, 2), round(y5, 2), round(mr5, 6), round(conv5, 6))
print("Point 2: ", round(x6, 2), round(y6, 2), round(mr6, 6), round(conv6, 6))

# Distance between points, Sphere to JTSK   
d3 = sqrt((x6 - x5)**2 + (y6 - y5)**2)
print("Distance between points: ", round((d3),2))


# Points coordinates by ČÚZK
x1cz = 1044339.83
y1cz = 741019.02
x2cz = 1044275.97
y2cz = 740793.83

# Distance between points, sphere to ČÚZK
d4 = sqrt((x2cz - x1cz)**2 + (y2cz - y1cz)**2)

print("Sphere to ČÚZK")
print('Distance between points: ', round(d4, 2))

# Differences in coordinates

#difference between distances in wgs and bessel point 1
dx1wb = x3 - x1 
dy1wb = y3 - y1

#difference between distances in bessel and wgs point 2
dx2wb = x4 - x2 
dy2wb = y4 - y2 # Bessel

#difference between distances in sphere and wgs
dx1ws = x5 - x1
dy1ws = y5 - y1

dx2ws = x6 - x2  
dy2ws = y6 - y2

#difference between distances in cuzk and wgs
dx1cz = x1cz - x1
dy1cz = y1cz - y1

dx2cz = x2cz - x2
dy2cz = y2cz - y2


print('Difference in coordinates between wgs and bessel')
print('p1',round(dx1wb,2), round(dy1wb,2))  
print('p2',round(dx2wb,2), round(dy2wb,2))
print('Difference in coordinates between wgs and sphere')
print('p1',round(dx1ws,2), round(dy1ws,2))
print('p2',round(dx2ws,2), round(dy2ws,2))
print('Difference in coordinates between wgs and cuzk')
print('p1',round(dx1cz,2), round(dy1cz,2))
print('p2',round(dx2cz,2), round(dy2cz,2))

#difference in distances between wgs and bessel
dd_wb = d2 - d1
dd_ws = d3 - d1
dd_wcz = d4 - d1
print('Difference in distances between wgs and bessel', round(dd_wb,2))
print('Difference in distances between wgs and sphere', round(dd_ws,2))
print('Difference in distances between wgs and cuzk', round(dd_wcz,2))

#Create a csv file to store the results
with open('results.csv', mode='w') as file:
    
    writer = csv.writer(file)
    writer.writerow(['Method', 'Point 1 x', 'Point 1 y', 'M', 'Conv', 'Distance'])
    writer.writerow(['WGS-84 to JTSK ', x1, y1, mr1, conv1])
    writer.writerow(['Bessel to JTSK ', x3, y3, mr3, conv3])
    writer.writerow(['Sphere to JTSK ', x5, y5, mr5, conv5])
    writer.writerow(['ČÚZK ', x1cz, y1cz])
    writer.writerow(['Method', 'Point 2 x', 'Point 2 y', 'M', 'Conv', 'Distance'])
    writer.writerow(['WGS-84 to JTSK ', x2, y2, mr2, conv2])
    writer.writerow(['Bessel to JTSK ', x4, y4, mr4, conv4])
    writer.writerow(['Sphere to JTSK ', x6, y6, mr6, conv6])
    writer.writerow(['ČÚZK ', x2cz, y2cz])

#Create a csv file to store the differences
with open('differences_distance.csv', mode='w') as file:
    
    writer = csv.writer(file)
    writer.writerow(['Method', 'Distances','Difference in distance'])
    writer.writerow(['WGS-84 and JTSK', round(d1, 2), '-'])
    writer.writerow(['WGS-84 and Bessel', round(d2, 2), round(dd_wb, 2)])
    writer.writerow(['WGS-84 and Sphere', round(d3, 2), round(dd_ws, 2)])
    writer.writerow(['WGS-84 and ČÚZK', round(d4, 2), round(dd_wcz, 2)]) 
    
with open('differences_coordinates.csv', mode='w') as file:
    
    writer = csv.writer(file)
    writer.writerow(['point', 'x coordinate', 'y coordinate'])
    writer.writerow(['point 1 WGS-84 - Bessel',round(dx1wb,2), round(dy1wb,2)])
    writer.writerow(['point 2 WGS-84 - Bessel',round(dx2wb,2), round(dy2wb,2)])
    writer.writerow(['point 1 WGS-84 - Sphere',round(dx1ws, 2), round(dy1ws, 2)])
    writer.writerow(['point 2 WGS-84 - Sphere',round(dx2ws, 2), round(dy2ws, 2)])
    writer.writerow(['point 1 WGS-84 - ČÚZK',round(dx1cz, 2), round(dy1cz, 2)])
    writer.writerow(['point 2  WGS-84 - ČÚZK',round(dx2cz, 2), round(dy2cz, 2)])
        
# Create points
point1_W = Point(-y1, -x1)
point2_W = Point(-y2, -x2)
point1_B = Point(-y3, -x3)
point2_B = Point(-y4, -x4)
point1_S = Point(-y5, -x5)
point2_S = Point(-y6, -x6)
point1_C = Point(-y1cz, -x1cz)
point2_C = Point(-y2cz, -x2cz)


# Create a GeoDataFrame
Points = gpd.GeoDataFrame(geometry=[point1_W, point2_W, point1_B, point2_B, point1_S, point2_S],
                       crs="EPSG:5514")
Points['point'] = ['P1', 'P2', 'P1', 'P2', 'P1', 'P2']
Points['type'] = ['WGS-84', 'WGS-84', 'Bessel', 'Bessel', 'Sphere', 'Sphere']

# Plot the points and save the map
cmap = ListedColormap(['red', 'blue', 'green', 'gold'])
Points.explore('type', cmap=cmap, tooltip='point',tiles='cartodbpositron').save('Points.html')
