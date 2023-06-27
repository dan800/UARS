# UARS
An archive of programs to open PROD files. These are available on the earthdata.nasa.gov website:

https://search.earthdata.nasa.gov/search?fi=HRDI

Some software is available here:

https://disc.gsfc.nasa.gov/information/documents?title=UARS%20Mission%20Preservation%20Documents

Including the UARS C Language Read Program: readuars.c
(http://docserver.gesdisc.eosdis.nasa.gov/repository/Mission/UARS/3.8_ScienceDataSoftwareTools/ReadSW/readuars.c)

To dump the contents of a binary PROD file to an ascii file, first compile the readuars.c code:
`
gcc -o readuars readuars.c -lm
`

Then run the executable:
`
./readuars -b HRDI_L3AT_SVOLER_A_D0057.V0011_C01_PROD > prod.txt
`

This produces the following:

```
-------------------------------
    SFDU Descriptor
-------------------------------
T[z] field : CCSD1Z000001
L[z] field : 00175132
T[i] field : NURS1I00HR04
L[i] field : 00175112

-------------------------------
    Header Record
-------------------------------
Satellite ID  : UARS
Record type   :  1
Instrument ID : HRDI        
Data subtype or species : VOLER_A     
Format version number   :    1
Physical record count   :        1
Number of continuation records for file label :    0
Number of physical records in file :      533
File Creation Time : 22-APR-1996 21:25:17.68
Year of 1st data record  :  91
Day of 1st data record   : 311
Milliseconds of 1st data record  :  1015808
Year of last data record :  91
Day of last data record  : 311
Milliseconds of last data record : 84115456
Data level : 3AT
UARS day number :   57
Number of data points per record :   33
Base index of data point value :    1
Record length in bytes :   328
CCB version number :        11
File cycle number      :     0
Virtual file flag :  
Total number of time/version entries in file :    0
Number of time/version entries in record :    0

---------------------------------------------
    Data Record Number        2
---------------------------------------------
Satellite ID  : UARS
Record type   :  3
Instrument ID : HRDI        
Physical record count :        2
Total number points in record : 33
Number of Actual points : 8
Starting index of first actual point : 20
Time in UDTF format :  91311  1015808 (yyddd, millisec)
Latitude  : 45.156013 (degrees)
Longitude : 137.727814 (degrees)
Local solar time   : 9.736342 (hours)
Solar zenith angle : 68.383324 (degrees)

Index  Altitude (km)  Data (1/cm**3/s)    Quality (1/cm**3/s)**2
-----  -------------  ------------------  ---------------------
   20             84        93232.5             2250.87
   21             87        113152.             2009.63
   22             90        120382.             2115.14
   23             93        83493.2             1644.25
   24             96        59886.0             1511.72
   25             99        43193.4             1201.75
   26            102        34115.7             975.170
   27            105        33333.5             1106.12
 ...
```

# useful links

https://catalog.data.gov/dataset/uars-high-resolution-doppler-imager-hrdi-level-3at-v011-uarhr3at-at-ges-disc
https://disc.gsfc.nasa.gov/datasets/UARHR3AT_011/summary
