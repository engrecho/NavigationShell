#! /usr/bin/env python
# -*- coding: utf-8 -*-
import gzip,os,sys
import glob,datetime


Jpath="16o" + os.sep + "*.16o.gz"


# /********1*********2*********3*********4*********5*********6*********7**
#  * Name:        xyz2llp
#  * Version:     1.0
#  * Author:      Jaylon.Wang
#  * Purpose:     Converts XYZ geocentric coordinates to Phi (latitude),
#  *              Lambda (longitude), H (height) referred to an
#  *              ellipsoid of semi-major axis A and flattening FL.
#  * Input:
#  * xyz           geocentric Cartesian coordinates [units are of distance]
import math
def xyz2llh(x, y, z) :

    A   = 6378137.0
    FL  = 0.00335281068118 
    B = A * (1.0 - FL) 
    if( z < 0 ):
        B= -1.0 * B
    r = math.sqrt( x*x + y*y )
    e = ( B*z - (A*A - B*B) ) / ( A*r )
    f = ( B*z + (A*A - B*B) ) / ( A*r )
    p= (4.0 / 3.0) * (e*f + 1.0) 
    q= 2.0 * (e*e - f*f) 
    d= p*p*p + q*q 

    if ( d >= 0 ) :
        v= pow( (math.sqrt( d ) - q), (1.0 / 3.0) ) - pow( (math.sqrt( d ) + q), (1.0 / 3.0) ) 
    else: 
        v= 2.0 * math.sqrt( -p ) * math.cos( math.acos( q/(p * math.sqrt( -p )) ) / 3.0 ) 
    
#  *   4.0 improve v
#  *       NOTE: not really necessary unless point is near pole
    if( v*v < math.fabs(p) ) :
        v= -(v*v*v + 2.0*q) / (3.0*p) 
        
    g   = (math.sqrt( e*e + v ) + e) / 2.0 
    t   = math.sqrt( g*g  + (f - v*g)/(2.0*g - e) ) - g 
    lat = math.atan( (A*(1.0 - t*t)) / (2.0*B*t) ) 

#  *   5.0 compute height above ellipsoid
    height= (r - A*t)*math.cos( lat ) + (z - B)*math.sin( lat ) 

#  *   6.0 compute longitude east of Greenwich
    zlong = math.atan2( y, x )
    if( zlong < 0 ):
        zlong= zlong + 2.0 * math.pi 

    lo = zlong 

#  *   7.0 convert latitude and longitude to degrees
    lat = lat * 360.0 / (2 * math.pi) 
    lo  = lo * 360.0 / (2 * math.pi) 
    if lo > 180:
        lo = lo - 360
    return (lat,lo,height) 
#End DEF

# print xyz2llh(-775858.6194, -4903039.9874, 3991748.6190)

import csv
print "This program gets Posion of the NOAA NGS's Osbservertion Station all over USA"
print "All the infomation is stored into OsbverStationPostion.csv"
print "The Format is ",'StationID', 'Latitude','Logitude','Height'
dir16o=glob.glob(Jpath)
if dir16o == [] :
    print "Please maka sure you have downloaded .16o.gz files into directory named 16o" 
    raw_input("Press Any Key to exit...")
    sys.exit()
# print dir16o
csvfile = open('OsbverStationPostion.csv', 'wb') 
txtfile = open('HTML.txt','wb')
spamwriter = csv.writer(csvfile,dialect='excel')
spamwriter.writerow(['StationID', 'Latitude','Logitude','Height'])
for f16o in dir16o:
    StaID=f16o[4:8]
    fid=gzip.GzipFile(f16o,'r')
    content = fid.readline()
    while content.find('END OF HEADER') < 0:
        if content.find('APPROX POSITION XYZ') > 0:
            x,y,z = content[0:60].split()
            lat,lon,h = xyz2llh(float(x), float(y), float(z))
            break
        content = fid.readline()
    #end while
    spamwriter.writerow([StaID, lat, lon , h])
    txtfile.write('        new BMap.Point(%.6f,%.6f),\n' %(lon,lat))
    print 'Dealing',StaID , 'lat,lon,h is',lat,lon,h 
#end for
spamwriter.writerow(['Created by Jaylon.Wang'])
csvfile.close()
txtfile.close()

