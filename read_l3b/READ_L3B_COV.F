CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c
c   Program name: READ_L3B_COV.F
c
c   Purpose:      read uars level 3b covariance prod files.
c
c   Author:       Ke - Jun Sun   Hughes STX
c
c   E-mail:       sun@daac.gsfc.nasa.gov
c
c   Usage:        READ_L3B_COV
c
c   Description:  This program allows a user to read the Fourier covariances
c                 stored in the UARS level 3B data files, and transform these
c                 to geophysically meaningful values defined on 4 degree
c                 latitude intervals. The program allows a user to input a
c                 UARS level 3B covariance data filename in the current
c                 directory, after which he has the option of dumping the
c                 results to the screen or saving the results in an ascii disk
c                 file. Options are provided for allowing selection of the
c                 latitude range of interest, and minimum and maximum pressure
c                 levels.
c
c   To compile :  f77 -o READ_L3B_COV READ_L3B_COV.F
c
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        CHARACTER FLNAME*70
        CHARACTER*3   OLDNEW
        CHARACTER*6   INSTR,PARAM
        CHARACTER*12  UNITS
        BYTE  SFDU_BYTE(86)
        INTEGER*4 LUN,INTIME(2),IORB_MODE,INSTR_MODE,IPR_MIN,IPR_MAX
        INTEGER*4 ITWO, IWR_TIME(2),IOS,IPR,NCOEF
        INTEGER*4 INTIME1(2),INTIME2(2),REC
        REAL*4 RLAT_MIN,RLAT_MAX,DY_ALT,RNT_ALT,STD_FAC,VER_NUM
        REAL*4 RLAT,SLAT,NLAT
        REAL*4 CE(169)
        REAL*4 THETA(360),Z(360)
        INTEGER*4 RESOL,NTHETA
        REAL*4 MINLAT,MAXLAT,MINPRS,MAXPRS,REQ_MIN_PRS,REQ_MAX_PRS,REC_PRS
        REAL*4 ZM
        INTEGER*4 MIN_PLV,MAX_PLV,REQ_MIN_PLV,REQ_MAX_PLV
        CHARACTER*70 OUTFILE
        CHARACTER*1 SVFILE,HITMISS
        INTEGER*4 COUNT/1/,I,SELECT,NMLEN,IRETURN,STATS(13)
        INTEGER*4 MULTI
        INTEGER*4 AS_DS,YR,M,D,JDAY
        CHARACTER*3 MONTHS(12)/
     &  'JAN','FEB','MAR','APR','MAY',
     &  'JUN','JUL','AUG','SEP','OCT','NOV','DEC'/
        CHARACTER*8 ALT_PRS
        CHARACTER*2 MB_KM


        PRINT *,' Enter UARS L3B Corvariance file name :'
        READ(5,'(A)') FLNAME

        REC=1
        OLDNEW='OLD'
        LUN=10
        CALL OPN3B_TWO(LUN,FLNAME,OLDNEW,IOS)
        IF(IOS .NE. 0) THEN
          PRINT *,' Can not open ',FLNAME
          GOTO 999
        ENDIF

c 	PRINT *, '********************************************'
c 	PRINT *, 'READING   ', FLNAME
c 	PRINT *, '********************************************'
c        PRINT *, FLNAME
        CALL RD3BCV_TWO_HDR(LUN,INSTR,PARAM,  
     &  IORB_MODE,INSTR_MODE,
     &  INTIME,SFDU_BYTE,ITWO,IWR_TIME,VER_NUM,IOS)

        IF((INTIME(1)/1000).EQ.0)THEN
         WRITE(6,*)'This is not a UARS L3B COVARIANCE file.'
         WRITE(6,*)'Please check file name and try again.'
         GOTO 999
        ENDIF
 
        PRINT *,'Do you wish to save result to an ASCII file (y/n) ?'
        READ(5,15) SVFILE
 15     FORMAT(A)
        IF(SVFILE.EQ.'Y' .OR. SVFILE.EQ.'y') THEN
          PRINT *,'Enter output file name, or <CR> for default name'
          READ(5,'(A)')OUTFILE
          IF(OUTFILE .EQ. '')THEN
            NMLEN=INDEX(FLNAME,' ')
            NMLEN=NMLEN-1
            OUTFILE=FLNAME(1:NMLEN)//'.DUMP'
            NMLEN=INDEX(OUTFILE,' ')
          ENDIF
          WRITE(6,*)'The result will be stored in ',OUTFILE
          WRITE(6,*)
        ENDIF

       MIN_PLV=100
       MAX_PLV=0
       MINLAT=90.0
       MAXLAT=-90.0

c       WRITE(6,*) 'IORB_MODE  : ',IORB_MODE
c       WRITE(6,*) 'INSTR_MODE : ',INSTR_MODE
c       WRITE(6,*) 'INTIME     : ',INTIME(1),INTIME(2)
C       WRITE(6,*) 'SFDU_BYTE  : ',SFDU_BYTE
c       WRITE(6,*) 'ITWO       : ',ITWO
c       WRITE(6,*) 'INSTR      : ',INSTR
c       WRITE(6,*) 'PARAM      : ',PARAM
c       WRITE(6,*) 'IWR_TIME   : ',IWR_TIME(1),IWR_TIME(2)
c       WRITE(6,*) 'VER_NUM    : ',VER_NUM

       MULTI = 1
       IF(INSTR.EQ.'HRDI'.OR.INSTR.EQ.'WINDII')THEN
         UNITS='(m/s)'
       ELSEIF(PARAM.EQ.'TEMP') THEN
         UNITS='(degrees K)'
       ELSEIF(PARAM.EQ.'CF2CL2'.OR.PARAM.EQ.'CFCL3'.OR.
     %        PARAM.EQ.'CH4'.OR.PARAM.EQ.'CO'.OR.
     %        PARAM.EQ.'CLONO2'.OR.PARAM.EQ.'H2O'.OR.
     %        PARAM.EQ.'HNO3'.OR.PARAM.EQ.'N2O5'.OR.
     %        PARAM.EQ.'N2O'.OR.PARAM.EQ.'NO2'.OR.
     %        PARAM.EQ.'NO'.OR.PARAM.EQ.'HCL'.OR.
     %        PARAM.EQ.'HF'.OR.PARAM.EQ.'CLO'.OR.
     %        PARAM.EQ.'O3_183'.OR.PARAM.EQ.'O3_205'.OR.
     %        PARAM.EQ.'O3B8'.OR.PARAM.EQ.'O3B9'.OR.
     %        PARAM.EQ.'O3') THEN
         UNITS = '(ppmv)'
         MULTI = 1.0E6 
       ELSE
         UNITS = '(1/km)'
       ENDIF

       JDAY = MOD(INTIME(1),1000)
       YR = INTIME(1)/1000 + 1900
       CALL JUL(YR,JDAY,M,D)

       IF(((INSTR.EQ.'HRDI').OR.(INSTR.EQ.'WINDII')).AND.
     %     ((PARAM.EQ.'MWIN_A').OR.(PARAM.EQ.'ZWIN_A')).OR.
     %     ((INSTR.EQ.'WINDII').AND.(PARAM.EQ.'TEMP'))) THEN
          ALT_PRS = 'altitude'
          MB_KM = 'km'
       ELSE
          ALT_PRS = 'pressure'
          MB_KM = 'mb'
       ENDIF


C    SCAN THRU THE WHOLE FILE TO GET MIN/MAX LAT AND PRESSURE

       DO 123 WHILE (1 .EQ. 1)
         CALL RD3BCV_TWO(LUN,RLAT,IPR,NCOEF,CE,INTIME1,INTIME2,IOS)
 
        IF (IOS .EQ. 0) THEN

         IF (IPR .GT. MAX_PLV) THEN
           MAX_PLV=IPR
          ELSEIF (IPR .LT. MIN_PLV) THEN
           MIN_PLV=IPR
         ENDIF
         IF (RLAT .GT. MAXLAT) THEN
          MAXLAT=RLAT
         ELSEIF (RLAT .LT. MINLAT) THEN
          MINLAT=RLAT
         ENDIF

         REC=REC + 1
        ELSE
         GOTO 888 
        ENDIF

  123  CONTINUE

 888    MINPRS=1000.0*10.0**(MAX_PLV/(-6.0))
        MAXPRS=1000.0*10.0**(MIN_PLV/(-6.0))
        WRITE(6,*)'****************************************************'
        WRITE(6,*)' Instrument    : ',INSTR
        WRITE(6,*)' Parameter     : ',PARAM
        WRITE(6,875)D,MONTHS(M),YR
 875    FORMAT('  Date          : ',I2.2,'-',A3,'-',I4)
        WRITE(6,*)' Time          : 12:00:00.00 GMT'
        IF(IORB_MODE .EQ. 0) THEN
         WRITE(6,*)' Orbit Node    : Ascending'
        ELSE
         WRITE(6,*)' Orbit Node    : Descending'
        END IF
        WRITE(6,501)MINLAT
 501    FORMAT("  Min. Latitude : ",F5.1," deg")
        WRITE(6,502)MAXLAT
 502    FORMAT("  Max. Latitude : ",F5.1," deg")
        WRITE(6,503)ALT_PRS,MINPRS,MB_KM
 503    FORMAT("  Min. ",A8," : ",G9.4,x,A2)
        WRITE(6,504)ALT_PRS,MAXPRS,MB_KM
 504    FORMAT("  Max. ",A8," : ",G9.4,x,A2)

        WRITE(6,*)'****************************************************'
        WRITE(6,*)'  '
 333    WRITE(6,*)'Enter northernmost latitude (degrees) :'
        READ(5,*,ERR=333) NLAT
        IF ((NLAT .LT. MINLAT).OR.(NLAT .GT. MAXLAT)) THEN
          WRITE(6,*)' Invalid latitude entered, try again. '
          GOTO 333
        ENDIF
        IF(MOD(NLAT,4) .NE. 0) THEN
         NLAT = NLAT + 1.0
         DO WHILE (MOD(INT(NLAT),4).NE.0)
          NLAT = NLAT + 1.0
         END DO
         NLAT = INT(NLAT)
         WRITE(6,334)NLAT
 334     FORMAT("Using next northernmost latitude at ",F5.1, " deg.")
        ENDIF

 444    WRITE(6,*)'Enter southernmost latitude (degrees) :'
        READ(5,*,ERR=444) SLAT
        IF((SLAT .LT. MINLAT).OR.(SLAT .GT. MAXLAT))THEN
          WRITE(6,*)' Invalid latitude entered, try again. '
          GOTO 444
        ELSEIF (SLAT .GT. NLAT)THEN
          WRITE(6,*)' Northernmost latitude is less than',
     %              ' southernmost latitude, try again.'
          GOTO 444
        ENDIF
        IF(MOD(SLAT,4) .NE. 0) THEN
         SLAT = SLAT - 1.0
         DO WHILE (MOD(INT(SLAT),4).NE.0)
          SLAT = SLAT - 1.0
         END DO
         SLAT = INT(SLAT)
         WRITE(6,445)SLAT
 445     FORMAT("Using next southernmost latitude at ",F5.1, " deg.")
        ENDIF

 555    WRITE(6,567)ALT_PRS,MB_KM
 567    FORMAT("Enter maximum ",A8," level (",A2,") :")
        READ(5,*,ERR=555) REQ_MAX_PRS
        IF (REQ_MAX_PRS .LT. MINPRS .OR. REQ_MAX_PRS .GT. MAXPRS) THEN
          WRITE(6,568)ALT_PRS
 568    FORMAT(" Invalid ",A8," entered, try again. ")
          GOTO 555
        ENDIF

        IF((REQ_MAX_PRS.NE.1000.0).AND.(REQ_MAX_PRS.NE.100.0).AND.
     %     (REQ_MAX_PRS.NE.10.0).AND.(REQ_MAX_PRS.NE.1.0)) THEN
          REQ_MIN_PLV =INT((-6)*LOG10(REQ_MAX_PRS/1000.0))
          REQ_MAX_PRS =1000.0*10.0**(REQ_MIN_PLV/(-6.0))
          WRITE(6,556)REQ_MAX_PRS
 556      FORMAT("Using next higher pressure level at ",G10.3, " mb.")
        ENDIF


 666    WRITE(6,667)ALT_PRS,MB_KM
 667    FORMAT("Enter minimum ",A8," level (",A2,") :")
        READ(5,*,ERR=666) REQ_MIN_PRS
        IF (REQ_MIN_PRS .LT. MINPRS .OR. REQ_MIN_PRS .GT. MAXPRS)THEN
          WRITE(6,668)ALT_PRS
 668      FORMAT(" Invalid ",A8," entered, try again.")
          GOTO 666
        ELSEIF (REQ_MIN_PRS .GT. REQ_MAX_PRS) THEN
          WRITE(6,669)ALT_PRS,ALT_PRS
 669      FORMAT("Minimum ",A8," should be less than ",
     %              "maximum ",A8,", try again.")
          GOTO 666
        ENDIF

        IF((REQ_MIN_PRS.NE.1000.0).AND.(REQ_MIN_PRS.NE.100.0).AND.
     %     (REQ_MIN_PRS.NE.10.0).AND.(REQ_MIN_PRS.NE.1.0)) THEN
          REQ_MAX_PLV =INT((-6)*LOG10(REQ_MIN_PRS/1000.0))+1
          REQ_MIN_PRS =1000.0*10.0**(REQ_MAX_PLV/(-6.0))
          WRITE(6,670)REQ_MIN_PRS
 670      FORMAT("Using next lower pressure level at ",G10.3, " mb.")
        ENDIF

        WRITE(6,*)'Data will start at -180.0 degrees longitude.'
 222    WRITE(6,*)'Enter longitudinal resolution between 1 and ',
     %  '5 degrees :'
        READ(5,'(I)',ERR=222) RESOL
        IF (RESOL .LT. 1 .OR. RESOL .GT. 5) THEN
           WRITE(6,*)' Invalid number, try again '
           GOTO 222
        ENDIF

        NTHETA= 360/RESOL

        DO 5 K = 1, NTHETA
          THETA(K) = RESOL*(K-1) - 180.0
    5   CONTINUE

        IF(SVFILE.EQ.'Y' .OR. SVFILE.EQ.'y') THEN
          OPEN (UNIT=21,FILE=OUTFILE,STATUS='UNKNOWN',ERR=98)
          WRITE(21,*) ''
          WRITE(21,*) 'File : ', FLNAME
        ELSE
          WRITE(6,*) ''
          WRITE(6,*) 'File : ', FLNAME
        END IF

C *************** Now Read it again ***********************

        REWIND LUN
        CLOSE (LUN)

        REC=0
        CALL OPN3B_TWO(LUN,FLNAME,OLDNEW,IOS)
        CALL RD3BCV_TWO_HDR(LUN,INSTR,PARAM,
     &  IORB_MODE,INSTR_MODE,
     &  INTIME,SFDU_BYTE,ITWO,IWR_TIME,VER_NUM,IOS)
   
        HITMISS = 'N'

       DO 234 WHILE (1 .EQ. 1)
         CALL RD3BCV_TWO(LUN,RLAT,IPR,NCOEF,CE,INTIME1,INTIME2,IOS)

        IF (IOS .EQ. 0) THEN
C         PRINT *, ''
c         PRINT *, '---------------------------'
c         PRINT *, 'DATA RECORDS',REC
c         PRINT *, '---------------------------'
c        WRITE(6,*) 'IPR       : ',IPR
         REC_PRS=1000.0*10.0**(IPR/(-6.0))
c        WRITE(6,*) 'NCOEF     : ',NCOEF
c        WRITE(6,*) 'RLAT      : ',RLAT
c        WRITE(6,*) 'INTIME1   : ',INTIME1(1),INTIME1(2)
c        WRITE(6,*) 'INTIME2   : ',INTIME2(1),INTIME2(2)
c         WRITE(6,*) 'CE(13 x 13): '
c         WRITE(6,*) (CE(N),N=1,169)
c         PRINT *, '---------------------------'
         IF(( RLAT .GE. SLAT) .AND. (RLAT.LE.NLAT) .AND.
C     %     (IPR .GE. REQ_MIN_PLV) .AND. (IPR .LE. REQ_MAX_PLV))THEN
     % (REC_PRS .GE. REQ_MIN_PRS) .AND. (REC_PRS .LE. REQ_MAX_PRS))THEN
         CALL ERROR_CV(CE,7,NTHETA,THETA,Z,ZM)  ! jj says 7 not 6
          HITMISS='Y'
          IF(SVFILE.EQ.'Y' .OR. SVFILE.EQ.'y') THEN
            WRITE(21,505) PARAM,UNITS,RLAT,ALT_PRS,REC_PRS,MB_KM
            WRITE(21,506) NTHETA,RESOL
              WRITE(21,507) (Z(N)*MULTI,N=1,NTHETA)

          ELSE
            WRITE(6,505) PARAM,UNITS,RLAT,ALT_PRS,REC_PRS,MB_KM
            WRITE(6,506) NTHETA,RESOL
 505       FORMAT(/,8('**********'),/,1X,A6,1X,A12,' at Latitude : ',
     %            F5.1,' and ',A8,' level : ',1PG10.3,' ',A2)
 506       FORMAT(/,I3,' data values beginning at longitude -180.0 ',
     %            'incremented every ',I1,' degrees :',/)

           WRITE(6,507) (Z(N)*MULTI,N=1,NTHETA)
 507   FORMAT(6(1X,1PE12.5))

             WRITE(6,*)

          ENDIF

         ENDIF

         REC=REC + 1
        ELSE

         GOTO 345 
        ENDIF

  234  CONTINUE
  345  CLOSE (LUN)
 
        IF(HITMISS.EQ.'N') THEN
          WRITE(6,777)ALT_PRS
 777      FORMAT('There are no data in the region and ',A8,
     %' level you selected.')
        ENDIF

        IF(SVFILE.EQ.'Y' .OR. SVFILE.EQ.'y') THEN
          CLOSE(21)                                                             
          OPEN (UNIT=21,FILE=OUTFILE,STATUS='OLD',ERR=99)
          IRETURN=FSTAT(21,STATS)
          WRITE(6,*) 'Total ',STATS(8),' bytes data are ',
     %       'saved in ',OUTFILE
          CLOSE (21)
        ENDIF
        GOTO 999
   98   WRITE(6,*)'Can not open ',OUTFILE
        GOTO 999
   99   WRITE(6,*)'Can not open existing ',OUTFILE
        GOTO 999

  999  STOP
       END

        SUBROUTINE OPN3B_TWO(LUN,FLNAME,OLDNEW,IOS)
C
C       THIS SUBROUTINE IS PART OF THE SOFTWARE TO GENERATE LEVEL
C       3B DATA FOR UARS
C       CONTACT: F.T. HUANG/SSAI 286-5908, 464-4751
C
C       THE ROUTINE OPENS THE FILE FOR WRITE WITH PROPER
C       ATTRIBUTES
C
C       BOTH ARGUMENTS ARE INPUT TO THE SUBROUTINE AND ARE THE
C       LOGICAL UNIT NUMBER (LUN, I*4) AND FILE NAME (FLNAME, CHAR *(*))
C
        CHARACTER *(*) FLNAME,OLDNEW
        INTEGER RECL_WORD/1000/
        IF(OLDNEW.EQ.'OLD') THEN
        OPEN(UNIT=LUN,NAME=FLNAME,TYPE=OLDNEW,READONLY,
     &  RECL=RECL_WORD,ORGANIZATION='SEQUENTIAL',ACCESS='SEQUENTIAL',
     &  RECORDTYPE='VARIABLE',FORM='UNFORMATTED',
     &  IOSTAT=IOS,ERR=200)
        ELSE
        OPEN(UNIT=LUN,NAME=FLNAME,TYPE=OLDNEW,
     &  RECL=RECL_WORD,ORGANIZATION='SEQUENTIAL',ACCESS='SEQUENTIAL',
     &  RECORDTYPE='VARIABLE',FORM='UNFORMATTED',
     &  IOSTAT=IOS,ERR=200)
        ENDIF
        RETURN
  200   CONTINUE
        CALL PERROR(IOS)
C       WRITE(6,*) 'ERR IN OPEN',IOS
  300   CONTINUE
        RETURN
        END


	SUBROUTINE RD3BCV_TWO_HDR(LUN,INSTR,PARAM,
     &	IORB_MODE,INSTR_MODE,
     &	INTIME,SFDU_BYTE,ITWO,IWR_TIME,VER_NUM,IOS)
C
C	THIS ROUTINE IS PART OF THE SOFTWARE TO GENERATE LEVEL 3B
C	FOR UARS.  IT WRITES THE TW_SIDED FOURIER COEEFICIENTS,
C	THE CORRESPONDING COVARIANCE MATRIX AND THE CONTROL
C	PARAMETERS WHICH WERE USED.
C
C	CONTACT F. T. HUANG/SSAI (301) 731-9300
C
C-----------------------------------------------------------------
C
C	ARGUMENT      TYPE        I/O      DESCRIPTION
C       ------------------------------------------------------------
C
C	LUN	      I*4	  I	   LOGICAL UNIT OF FILE
C	INSTR	      CHAR*6	  O	   INSTRMENT ID
C	PARAM	      CHAR*6      O	   MEASURED PARAMETER ID
C	IORB_MODE     I*4	  O	   1 ASCENDING MODE, 0 DESCEND
C	INSTR_MODE    I*4	  O        INSTRUMENT MODE AS DETERMINED
C					   FROM PARAMETER FILE.
C	INTIME	      I*4(2)	  O	   TIME OF OUTPUT OF COEFFS
C					   UDTF (DAY, 12:00GMT)
C	SFDU_BYTE     C(86)       O        SFDU INFORMATION
C	ITWO	      I*4         O        DENOTES IF FORWARD OR
C					   TWO-SIDED FILE (1:FRWRD,2:TWO)
C	IWR_TIME      I*4(2)	  O	   TIME IN UDTF WHEN FILE WAS CREATED.
C	VER_NUM	      R*4	  O	   L3B SOFTWARE VERSION NUMBER.
C	IOS           I*4	  O	   FORTRAN ERROR STATUS
C-----------------------------------------------------------------------	
	CHARACTER*(*) INSTR,PARAM
	INTEGER*4 INTIME(2),IWR_TIME(2),IREC_TYPE/0/
	BYTE SFDU_BYTE(86)
C
C	READ 86 BYTES OF SFDU INFORMATION
C
	ILTDY1 = 1
	READ(LUN,KEY=ILTDY1,KEYID=0,ERR=200,IOSTAT=IOS) ILTDY,SFDU_BYTE
C	WRITE(6,100) SFDU_BYTE
  100	FORMAT(1X,50A1)
C
C	READ 'HEADER'
C
	ILTDY1 = 2
C
C	NEW NEWER LONGER RECORDS (VERSION 3.06 OR LATER)
C	IF ERROR 67 THEN DATA OVERRUN. TRY SHORTER OLDER RECORDS
C
	IF(IREC_TYPE.EQ.0) THEN
	READ (LUN,KEY=ILTDY1,KEYID=0,ERR=200,IOSTAT=IOS)
     &	ILTDY,
     &	IORB_MODE,INSTR_MODE,ITWO,
     &	INTIME(1),INTIME(2),
     &	(INSTR(L:L),L=1,6),(PARAM(L:L),L=1,6),
     &	IWR_TIME(1),IWR_TIME(2),VER_NUM
  200	CONTINUE
	ENDIF
C
C	TRY READING OLDER RECORD TYPE IF ERROR IS 67
C
	IF(IOS.EQ.67.OR.IREC_TYPE.EQ.1) THEN
C
	IWR_TIME(1) = 0
	IWR_TIME(2) = 0 
	VER_NUM = 0.0
C
	IREC_TYPE = 1
	READ (LUN,KEY=ILTDY1,KEYID=0,ERR=300,IOSTAT=IOS)
     &	ILTDY,
     &	IORB_MODE,INSTR_MODE,ITWO,
     &	INTIME(1),INTIME(2),
     &	(INSTR(L:L),L=1,6),(PARAM(L:L),L=1,6)
	ENDIF
  300	CONTINUE
	RETURN
	END

	SUBROUTINE RD3BCV_TWO(LUN,RLAT,IPR,NCOEF,CE,
     &	INTIME1,INTIME2,IOS)
C
C	THIS ROUTINE IS PART OF THE SOFTWARE TO GENERATE LEVEL 3B
C	FOR UARS.  IT WRITES THE TW_SIDED FOURIER COEEFICIENTS,
C	THE CORRESPONDING COVARAINCE MATRIX AND OTHE CONTROL
C	PARAMETERS WHICH WERE USED.
C
C	CONTACT F. T. HUANG/SSAI (301) 464-4751
C
C-----------------------------------------------------------------
C
C	ARGUMENT      TYPE        I/O      DESCRIPTION
C       ------------------------------------------------------------
C
C	LUN	      I*4	  I	   LOGICAL UNIT OF FILE
C	RLAT	      R*4	  I	   LATITUDE OF COEFFS
C	IPR	      I*4	  I	   PRESSURE LEVEL INDEX
C	CE	      R*4 ARRAY	  I        COVARIANCE MATRIX
C		      (NCOEF BY NCOEF)
C	INTIME1	      I*4(2)	  I	   TIME OF LAST DATA POINT
C					   USED IN FORWARD RECURSION
C	INTIME2	      I*4(2)	  I	   TIME OF LAST DATA POINT USED
C				  	   FOR BACKWARD RECURSION.
C-----------------------------------------------------------------------	
	INTEGER*4 INTIME1(2),INTIME2(2)
	REAL*4 CE(13,13)
C
C	READ DATA. FIRST FIND KEY FROM RLAT.
C
	ILAT =  1 + (88.0 - RLAT)/4.0 
	ILTDY1 = ILAT * 50 + IPR
	READ (LUN,KEY=ILTDY1,KEYID=0,ERR=200,IOSTAT=IOS)
     &	ILTDY,IPR,
     &	NCOEF,RLAT,
     &	INTIME1(1),INTIME1(2),INTIME2(1),INTIME2(2),
     &	((CE(I,J),I=1,NCOEF),J=1,NCOEF)
  200	CONTINUE
	RETURN
	END


	SUBROUTINE ERROR_CV(COVAR,NCOF21,NTHETA,THETA,Z,ZMEAN)
C
C	THIS ROUTINE CALCULATES ERRORS (FROM VARIANCES) FROM THE INPUT
C	LEVEL 3B COVARIANCE MATRIX COVAR AT THE ANGLES THETA IN DEGREES.
C
C	NCOF21 IS SUCH  THAT 2*NCOF21 + 1 ARE THE TOTAL NUMBER
C	OF COEFFICIENTS. 
C	
C	THE FIRST NCOF21 ARE THE CONSTANT AND COSINE COFS
C	AND THE LAST NCOF21 ARE THE SINE COFS.
C
C	FOR THE UARS CASE, THERE NCOF21 IS 6, AND THE FIRST T COEFFS
C	ARE THE CONSTANT AND SIX COSINE COEFFS.  THE LAST 6 ARE THE
C	CORRESPONDING SINE COEFFS.
C
C	THETA(I) IS THE ARRAY OF ANGLES (FOR UARS THEY ARE 
C	LONGITUDES) IN DEGREES, NTHETA IS THE NUMBER OF
C
C	ARGUMENT	I/O	TYPE		DESCRIPTION
C	--------	---	----		-----------------------
C	COVAR 	        I       R*4 ARRAY	COVARIANCE ARRAY.
C				(2*NCOF21+1)x	FROM THE DIAGONAL TERMS
C                               (2*NCOF21+1)    ONE CAN GET THE MEAN ERROR.
C
C	NCOF21		I	I*4		DENOTES TOTAL NUMBER OF
C						COSINE TERMS, INCUDING
C						CONSTANT TERM.
C
C	NTHETA		I	I*4		NUMBER OF ANGLES 
C	
C	THETA		I	R*4(NTHETA)	ANGLES IN DEGREES WHERE
C						SERIES IS TO BE EVALUATED.
C
C	Z		O	R*4(NTHETA)	RESULTS OF ERRORS FOR EACH
C						NTHETA ANGLES.
C
C       BY JAMES JOHNSON, HUGHES STX (1996)
C       MODIFIED FROM SUBROUTINE FOURIER_SUM
C-------------------------------------------------------------------------
C
	DIMENSION Z(1),THETA(1),COVAR(13,13)
        
        REAL*4 PI, TWPI, CVK(13)

        PI = ACOS(-1.0)
        TWPI = 2.0*PI

C
C	DO OVER NTHETA ANGLES.
C
	DO 100 I=1,NTHETA
C
C	FIND ARGUMENTS FOR COS AND SINES
C

	ARG=TWPI * THETA(I)/360.0
	SUM = 0.0
C
C	FIND COSINES AND SINES AND FOURIER BASIS (EXPANSION) VECTOR.
C
        CVK(1) = 1.0
	DO 50 J=2,NCOF21
	  ARG1 = ARG * (J - 1)
	  COS1=COS(ARG1)
	  SIN1=SIN(ARG1)
          CVK(J) = COS1
          CVK(J+6) = SIN1
   50	CONTINUE

C
C       COMPUTE THE ERRORS = SQRT(VARIANCE)
C
        DO 70 N=1,13
          DO 60 M=1,13
	    SUM = SUM + CVK(M)*COVAR(M,N)*CVK(N)
   60     CONTINUE
   70   CONTINUE
	Z(I)=SQRT(SUM)
  100	CONTINUE

C
C       COMPUTE THE MEAN OF THE ERRORS
C
        ZMEAN=0.0
        DO 110 N=2,13
          ZMEAN = ZMEAN + COVAR(N,N) 
  110   CONTINUE
        ZMEAN = SQRT(COVAR(1,1) + 0.5*ZMEAN)
        
	RETURN
	END
C ****************************************************************************
C
C Name:         JUL
C
C Project:      PCDS Version 1.0
C
C Computer:     VAX/780 -- VAX/VMS Version 2.0
C
C Language:     FORTRAN 77
C
C Programmer:   L. SMIRNOVA
C
C Purpose:      Convert julian day to month and day.
C
C Calling arguments: YEAR,JULIAN,MONTH,DAY
C
C Calling routines: CONV
C
C Subroutines:  NONE
C
C
C Commons:      None
C
C Termination conditions:  None
C
C Inputs:       YEAR AND JULIAN DAY
C
C Outputs:      MONTH AND DAY
C
C
C METHOD:
c
C PDL:
C               PROC JUL:
C                       CHECK IF THE YEAR IS A LEAP YEAR
C                       DO UNTIL END OF MNTAB TABLE
C                          SUMMARIZE THE DAYS
C                          IF SUM = JULIAN DAY, THEN
C                             GET THE MONTH
C                             GET THE DAY
C                          END IF
C                       ENDDO
C               END
C
C******************************************************************************
        SUBROUTINE JUL(IYEAR,JULIAN,MONTH,IDAY)
        INTEGER MNTAB(12)
        DATA MNTAB/31,28,31,30,31,30,31,31,30,31,30,31/
        MNTAB(2)=28
        IF(MOD(IYEAR,4).EQ.0) MNTAB(2)=29
        MONTH=0
        IDAY=0
        ISUM=0
        DO I=1,12
           ISUM=ISUM+MNTAB(I)
           IF(JULIAN.LE.ISUM) THEN
              MONTH=I
              ISUM=ISUM-MNTAB(I)
              IDAY=JULIAN-ISUM
              GO TO 100
           END IF
        END DO
 100    CONTINUE
        RETURN
        END

