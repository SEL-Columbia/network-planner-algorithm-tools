#!/usr/bin/env python


import csv
import os
 

outputDir=r"C:\Users\Selin\Desktop\stat" 
	    

statsFile= outputDir + os.sep + "statTemp.csv"
statsFile2= outputDir + os.sep + "statTemp2.csv"
csvReader= csv.reader(open(statsFile, 'rb'))

for row in csvReader:
    for a in row:
        print a


#csvReader.readrow([x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12])

#csvWriter = csv.writer(open(statsFile2, 'wb'))

#csvWriter.writerow([x, y, size, MID])


           
        
	        
	
        
     
     
        
       
