CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c
c   Program name: READ_L3B_COEF.F
c
c   Purpose:      read uars level 3b coefficient prod files.
c
c   Author:       Ke - Jun Sun   Hughes STX
c
c   E-mail:       sun@daac.gsfc.nasa.gov
c
c   Usage:        READ_L3B_COEF
c
c   Description:  This program allows a user to read the Fourier coefficients
c                 stored in the UARS level 3B data files, and transform these
c                 to geophysically meaningful values defined on 4 degree
c                 latitude intervals. The program allows a user to input a
c                 UARS level 3B coefficient data filename in the current
c                 directory, after which he has the option of dumping the
c                 results to the screen or saving the results in an ascii disk
c                 file. Options are provided for allowing selection of the
c                 longitudinal resolution, latitude range of interest, and
c                 minimum and maximum pressure levels.
c
c   To compile :  f77 -o READ_L3B_COEF READ_L3B_COEF.F 
c
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        CHARACTER     FLNAME*70
        CHARACTER*3   OLDNEW
        CHARACTER*6   INSTR,PARAM
        CHARACTER*12  PARAMS_20,UNITS
        CHARACTER*1   SVFILE,HITMISS
        BYTE  SFDU_BYTE(86)
        INTEGER*4 LUN,INTIME(2),IORB_MODE,INSTR_MODE,IPR_MIN,IPR_MAX
        INTEGER*4 ITWO,IDYLAG,IGAP_SIZE,NPTS_LSQ_MIN,NPTS_INIT,IVAR_SEL
        INTEGER*4 IHDR_FLAG,IWR_TIME(2),NPTS_BKWRD_MIN
        INTEGER*4 NPTS_BKWRD_DAY_MIN,IOS
        INTEGER*4 IPR,NCOEF,NPTSF,NPTSB,NPTS_DAY,NPTS_BAD
        INTEGER*4 INTIME1(2),INTIME2(2),NPTS_DEV,IFLAG_FW,IFLAG
        INTEGER*4 ISTD_FLAG,IVER,ICYC,ISTRT(2),ISTOP(2)
        REAL*4 RLAT_MIN,RLAT_MAX,DY_ALT,RNT_ALT,STD_FAC,VER_NUM
        REAL*4 RLAT,COEF(13),XLONG(90),TAUINV,TAVE,AVE_LST,AVE_SZA,RES_DEV1
        REAL*4 STDFAC,STD_CALC,SEC_FST,PARAM_VALUES(90),STD_OLD(90)
        REAL*4 THETA(360),Y(360),SLAT,NLAT
        REAL*4 MINLAT,MAXLAT,MINPRS,MAXPRS,REQ_MIN_PRS,REQ_MAX_PRS,REC_PRS
        INTEGER*4 RESOL,NTHETA,MIN_PLV,MAX_PLV,REQ_MIN_PLV,REQ_MAX_PLV
        CHARACTER*70 OUTFILE
        INTEGER*4 COUNT/1/,I,SELECT,NMLEN,IRETURN,STATS(13)
        INTEGER*4 MULTI
        INTEGER*4 AS_DS,YR,M,D,JDAY
        CHARACTER*3 MONTHS(12)/
     &  'JAN','FEB','MAR','APR','MAY',
     &  'JUN','JUL','AUG','SEP','OCT','NOV','DEC'/
        CHARACTER*8 ALT_PRS
        CHARACTER*2 MB_KM

        WRITE(6,*)'Enter UARS L3B Coefficient file name :'
        READ(5,'(A)')FLNAME

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
        CALL RD3B_TWO_HDR(LUN,INSTR,PARAM,INTIME,
     &  IORB_MODE,INSTR_MODE,RLAT_MIN,RLAT_MAX,IPR_MIN,IPR_MAX,
     &  SFDU_BYTE,ITWO,DY_ALT,RNT_ALT,
     &  IDYLAG,IGAP_SIZE,NPTS_LSQ_MIN,NPTS_INIT,
     &  PARAMS_20,IVAR_SEL,STDFAC,
     &  IHDR_FLAG,IWR_TIME,VER_NUM,
     &  NPTS_BKWRD_MIN,NPTS_BKWRD_DAY_MIN,ISTRT,ISTOP,IOS)

        IF(IOS.NE.0)THEN
         WRITE(6,*)'This is not a UARS L3B COEFFICIENT file.'
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
c       WRITE(6,*) 'INSTR      : ',INSTR
c       WRITE(6,*) 'PARAM      : ',PARAM
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

C SCAN THRU THE WHOLE FILE TO GET MIN/MAX LAT AND PRESSURE

 	DO 20 I=IFIX(RLAT_MAX),IFIX(RLAT_MIN),-4 
         DO 10 J=IPR_MIN,IPR_MIN+IPR_MAX-1
          CALL RD3B_TWO(LUN,RLAT,IPR,COEF,NCOEF,
     &    NPTSF,NPTSB,NPTS_DAY,NPTS_BAD,
     &    TAUINV,INTIME1,INTIME2,
     &    XLONG,TAVE,AVE_LST,AVE_SZA,
     &    RES_DEV1,STD_CALC,NPTS_DEV,
     &    SEC_FST,IFLAG_FW,IFLAG,PARAM_VALUES,ISTD_FLAG,
     &    IVER,ICYC,
     &    IGMT,STD_OLD,IOS)
 
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

        ENDIF

   10	CONTINUE
   20	CONTINUE

        MINPRS=1000.0*10.0**(MAX_PLV/(-6.0))
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
 501    FORMAT("  Min. latitude : ",F5.1," deg")
        WRITE(6,502)MAXLAT
 502    FORMAT("  Max. latitude : ",F5.1," deg")
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
          WRITE(6,556)ALT_PRS,REQ_MAX_PRS,MB_KM
 556      FORMAT("Using next higher ",A8," level at ",G10.3,x,A2)
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
          WRITE(6,670)ALT_PRS,REQ_MIN_PRS,MB_KM
 670      FORMAT("Using next lower ",A8," level at ",G10.3,x,A2,"." )
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

        CALL OPN3B_TWO(LUN,FLNAME,OLDNEW,IOS)

        CALL RD3B_TWO_HDR(LUN,INSTR,PARAM,INTIME,
     &  IORB_MODE,INSTR_MODE,RLAT_MIN,RLAT_MAX,IPR_MIN,IPR_MAX,
     &  SFDU_BYTE,ITWO,DY_ALT,RNT_ALT,
     &  IDYLAG,IGAP_SIZE,NPTS_LSQ_MIN,NPTS_INIT,
     &  PARAMS_20,IVAR_SEL,STDFAC,
     &  IHDR_FLAG,IWR_TIME,VER_NUM,
     &  NPTS_BKWRD_MIN,NPTS_BKWRD_DAY_MIN,ISTRT,ISTOP,IOS)

C ****** UNMASK THE FOLLOWING LINES IF YOU WANT TO SEE THE HEADER
C
c       WRITE(6,*) ' ***** FROM HEADER *****'
c       WRITE(6,*) 'Min. Latitude   : ',RLAT_MIN
c       WRITE(6,*) 'Max. Latitude   : ',RLAT_MAX
c       WRITE(6,*) 'Min. Plvl       : ',IPR_MIN
c       WRITE(6,*) 'Max. Plvl       : ',IPR_MAX
c       WRITE(6,*) 'IORB_MODE  : ',IORB_MODE
c       WRITE(6,*) 'INSTR_MODE : ',INSTR_MODE
c       WRITE(6,*) 'INTIME     : ',INTIME(1),INTIME(2)
c       WRITE(6,*) 'ITWO       : ',ITWO
c       WRITE(6,*) 'INSTR      : ',INSTR
c       WRITE(6,*) 'PARAM      : ',PARAM
c       WRITE(6,*) 'DY_ALT     : ',DY_ALT
c       WRITE(6,*) 'RNT_ALT    : ',RNT_ALT
c       WRITE(6,*) 'IDYLAG     : ',IDYLAG
c       WRITE(6,*) 'IGAP_SIZE  : ',IGAP_SIZE
c       WRITE(6,*) 'NPTS_LSQ_MIN : ',NPTS_LSQ_MIN
c       WRITE(6,*) 'NPTS_INIT  : ',NPTS_INIT
c       WRITE(6,*) 'PARAMS_20  : ',PARAMS_20
c       WRITE(6,*) 'IVAR_SEL   : ',IVAR_SEL
c       WRITE(6,*) 'STDFAC     : ',STDFAC
c       WRITE(6,*) 'IHDR_FLAG  : ',IHDR_FLAG
c       WRITE(6,*) 'IWR_TIME   : ',IWR_TIME(1),IWR_TIME(2)
c       WRITE(6,*) 'VER_NUM    : ',VER_NUM
c       WRITE(6,*) 'NPTS_BKWRD_MIN     : ',NPTS_BKWRD_MIN
c       WRITE(6,*) 'NPTS_BKWRD_DAY_MIN : ',NPTS_BKWRD_DAY_MIN
c       WRITE(6,*) 'ISTRT      : ',ISTRT(1),ISTRT(2)
c       WRITE(6,*) 'ISTOP      : ',ISTOP(1),ISTOP(2)

c        PRINT *, '---------------------------'
c        PRINT *, 'DATA RECORDS'
c        PRINT *, '---------------------------'

        HITMISS = 'N'

        DO 40 I=IFIX(RLAT_MAX),IFIX(RLAT_MIN),-4
         DO 25 J=IPR_MIN,IPR_MIN+IPR_MAX-1
          CALL RD3B_TWO(LUN,RLAT,IPR,COEF,NCOEF,
     &    NPTSF,NPTSB,NPTS_DAY,NPTS_BAD,
     &    TAUINV,INTIME1,INTIME2,
     &    XLONG,TAVE,AVE_LST,AVE_SZA,
     &    RES_DEV1,STD_CALC,NPTS_DEV,
     &    SEC_FST,IFLAG_FW,IFLAG,PARAM_VALUES,ISTD_FLAG,
     &    IVER,ICYC,
     &    IGMT,STD_OLD,IOS)
 
         IF (IOS .EQ. 0) THEN
C          PRINT *, ''
c          PRINT *, '---------------------------'
c          PRINT *, 'DATA RECORDS'
c          PRINT *, '---------------------------'

C  UNMASK THE FOLLOWING LINES IF YOU WANT TO SEE EVERY VALUE

c       WRITE(6,*) 'IPR       : ',IPR
        REC_PRS=1000.0*10.0**(IPR/(-6.0))
c       WRITE(6,*) 'NCOEF     : ',NCOEF
c       WRITE(6,*) 'RLAT      : ',RLAT
c       WRITE(6,*) 'NPTSF     : ',NPTSF
c       WRITE(6,*) 'NPTSB     : ',NPTSB
c       WRITE(6,*) 'NPTS_DAY  : ',NPTS_DAY
c       WRITE(6,*) 'NPTS_BAD  : ',NPTS_BAD
c       WRITE(6,*) 'TAUINV    : ',TAUINV
c       WRITE(6,*) 'INTIME1   : ',INTIME1(1),INTIME1(2)
c       WRITE(6,*) 'INTIME2   : ',INTIME2(1),INTIME2(2)
c       WRITE(6,*) 'TAVE      : ',TAVE
c       WRITE(6,*) 'AVE_LST   : ',AVE_LST
c       WRITE(6,*) 'AVE_SZA   : ',AVE_SZA
c       WRITE(6,*) 'RES_DEV1  : ',RES_DEV1
c       WRITE(6,*) 'STD_CALC  : ',STD_CALC
c       WRITE(6,*) 'NPTS_DEV  : ',NPTS_DEV
c       WRITE(6,*) 'SEC_FST   : ',SEC_FST
c       WRITE(6,*) 'IFLAG_FW  : ',IFLAG_FW
c       WRITE(6,*) 'IFLAG     : ',IFLAG
c       WRITE(6,*) ''
c       WRITE(6,*) 'COEF : '
c       WRITE(6,*) (COEF (K),K=1,NCOEF)
c       WRITE(6,*) ''
c       WRITE(6,*) 'XLONG : '
c       WRITE(6,*) (XLONG(K),K=1,NPTS_DAY)
c       WRITE(6,*) ''
c       WRITE(6,*) 'PARAM_VALUES : '
c       WRITE(6,*) (PARAM_VALUES(K),K=1,NPTS_DAY)
c       WRITE(6,*) ''
c       WRITE(6,*) 'ISTD_FLAG : ',ISTD_FLAG
c       WRITE(6,*) 'IVER      : ',IVER
c       WRITE(6,*) 'ICYC      : ',ICYC
c       WRITE(6,*) ''
c       WRITE(6,*) 'STD_OLD : '
c       WRITE(6,*) (STD_OLD(K),K=1,NPTS_DAY)
c         PRINT *, '---------------------------'

         IF(( RLAT .GE. SLAT) .AND. (RLAT.LE.NLAT) .AND.
     % (REC_PRS .GE. REQ_MIN_PRS) .AND. (REC_PRS .LE. REQ_MAX_PRS))THEN  
         CALL FOURIER_SUM(COEF,7,NTHETA,THETA,Y)  ! jj says 7 not 6
          HITMISS='Y'
          IF(SVFILE.EQ.'Y' .OR. SVFILE.EQ.'y') THEN
            WRITE(21,505) PARAM,UNITS,RLAT,ALT_PRS,REC_PRS,MB_KM
            WRITE(21,506) NTHETA,RESOL
            WRITE(21,507) (Y(N)*MULTI,N=1,NTHETA)

          ELSE
            WRITE(6,505) PARAM,UNITS,RLAT,ALT_PRS,REC_PRS,MB_KM
 505       FORMAT(/,8('**********'),/,1X,A6,1X,A12,' at Latitude : ',
     %            F5.1,' and ',A8,' level : ',1PG10.3,' ',A2)
            WRITE(6,506) NTHETA,RESOL
 506       FORMAT(/,I3,' data values beginning at longitude -180.0 ',
     %            'incremented every ',I1,' degrees :',/)
           WRITE(6,507) (Y(N)*MULTI,N=1,NTHETA)

 507   FORMAT(6(1X,1PE12.5))

             WRITE(6,*)
          ENDIF

         ENDIF

        ENDIF

   25   CONTINUE
   40   CONTINUE

        CLOSE (LUN)
        
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
  999   STOP
        END

	SUBROUTINE FOURIER_SUM(COF,NCOF21,NTHETA,THETA,Y)
C
C	THIS ROUTINE CALCULATES A FOURIER SERIES GIVEN THE INPUT
C	COEFFICIENTS COF AND THE ANGLES THETA IN DEGREES.
C
C	NCOF21 IS SUCH  THAT 2*NCOF21 + 1 ARE THE TOTAL NUMBER
C	OF COEFFICIENTS. 
C	
C	THE FIRST NCOF21 ARE THE CONSTANT AND COSINE COFS
C	AND THE LAST NCOF21 ARE THE SINE COFS.
C
C	FOR THE UARS CASE, THERE NCOF21 IS 7, AND THE FIRST T COEFFS
C	ARE THE CONSTANT AND SIX COSINE COEFFS.  THE LAST 6 ARE THE
C	CORRESPONDING SINE COEFFS.
C
C	THETA(I) IS THE ARRAY OF ANGLES (FOR UARS THEY ARE 
C	LONGITUDES) IN DEGREES, NTHETA IS THE NUMBER OF
C
C	ARGUMENT	I/O	TYPE		DESCRIPTION
C	--------	---	----		-----------------------
C	COF		I	R*4 ARRAY	ARRAY OF FOURIER COEFFS.
C				(2*NCOF21+1)	FIRST NCO21+1 NUMBER ARE
C						COSINE TERMS.  NEXT
C						NCOF21 ARE SINE TERMS.
C	NCOF21		I	I*4		DENOTES TOTAL NUMBER OF
C						COSINE TERMS, INCUDING
C						CONSTANT TERM.
C	NTHETA		I	I*4		NUMBER OF ANGLES 
C	
C	THETA		I	R*4(NTHETA)	ANGLES IN DEGREES WHERE
C						SERIES IS TO BE
C						EVALUATED.
C	Y		O	R*4(NTHETA)	ANSWERS. RESULTS OF 
C						FOURIER SERIES FOR EACH OF
C						NTHETA ANGLES.
C-------------------------------------------------------------------------
C
	DIMENSION Y(1),COF(1),THETA(1)
	DATA PI/3.141592/
	DATA TWPI/6.283192/
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
C	FIND COSINES AND SINES AND SUM FOR HIGHER TERMS.
C
	DO 50 J=2,NCOF21
	ARG1 = ARG * (J - 1)
	COS1=COS(ARG1)
	SIN1=SIN(ARG1)
	SUM = SUM + COF(J) * COS1 + COF(NCOF21 + J - 1) * SIN1
   50	CONTINUE
C
C	ADD CONSTANT (FIRST TERM)
C
	Y(I)=SUM+COF(1)
C
  100	CONTINUE
	RETURN
	END

	SUBROUTINE OPN3B_TWO(LUN,FLNAME,OLDNEW,IOS)
C
C	THIS SUBROUTINE IS PART OF THE SOFTWARE TO GENERATE LEVEL
C	3B DATA FOR UARS
C	CONTACT: F.T. HUANG/SSAI 286-5908, 464-4751
C
C	THE ROUTINE OPENS THE FILE FOR WRITE WITH PROPER
C	ATTRIBUTES
C
C	BOTH ARGUMENTS ARE INPUT TO THE SUBROUTINE AND ARE THE
C	LOGICAL UNIT NUMBER (LUN, I*4) AND FILE NAME (FLNAME, CHAR *(*))
C
	CHARACTER *(*) FLNAME,OLDNEW
	INTEGER RECL_WORD/1000/
	IF(OLDNEW.EQ.'OLD') THEN
     	OPEN(UNIT=LUN,NAME=FLNAME,TYPE=OLDNEW,READONLY,
     &	RECL=RECL_WORD,ORGANIZATION='SEQUENTIAL',ACCESS='SEQUENTIAL',
     &	RECORDTYPE='VARIABLE',FORM='UNFORMATTED',
     &	IOSTAT=IOS,ERR=200)
	ELSE
     	OPEN(UNIT=LUN,NAME=FLNAME,TYPE=OLDNEW,
     &	RECL=RECL_WORD,ORGANIZATION='SEQUENTIAL',ACCESS='SEQUENTIAL',
     &	RECORDTYPE='VARIABLE',FORM='UNFORMATTED',
     &	IOSTAT=IOS,ERR=200)
	ENDIF
	RETURN
  200	CONTINUE
        CALL PERROR(IOS)
C       WRITE(6,*) 'ERR IN OPEN',IOS
  300	CONTINUE
	RETURN
	END

	SUBROUTINE RD3B_TWO_HDR(LUN,INSTR,PARAM,INTIME,
     &	IORB_MODE,INSTR_MODE,RLAT_MIN,RLAT_MAX,IPR_MIN,IPR_MAX,
     &  SFDU_BYTE,ITWO,DY_ALT,RNT_ALT,
     &  IDYLAG,IGAP_SIZE,NPTS_LSQ_MIN,NPTS_INIT,
     &  PARAMS_20,IVAR_SEL,STDFAC,
     &  IHDR_FLAG,IWR_TIME,VER_NUM,
     &	NPTS_BKWRD_MIN,NPTS_BKWRD_DAY_MIN,ISTRT,ISTOP,IOS)
C
C	THIS ROUTINE IS PART OF THE SOFTWARE OF LEVEL 3B
C	DATA FOR UARS. IT READS THE SFDU AND 'HEADER' INFORMATION
C	FOR LEVLE 3B FILES.
C
C	CONTACT F.T. HUANG/SSAI  (301) 464-4751
C	FEB 1, 1991.
C
C-------------------------------------------------------------------
C
C	ARGUMENT      TYPE        I/O      DESCRIPTION
C       ------------------------------------------------------------
C
C	LUN	      I*4	  I	   LOGICAL UNIT OF FILE
C	INSTR	      CHAR*6	  O	   INSTRUMENT ID
C	PARAM	      CHAR*6      O	   MEASURED PARAMETER
C	INTIME        I*4(2)	  O	   UDTF (YRDAY,MSEC)
C					   CORRESPONDING
C					   TO COEFF WRITE TIME (12:00)
C	IORB_MODE     I*4	  O	   0 ASCENDING MODE, 1 DESCEND
C					   2 COMBINED (NOT USED)
C	INSTR_MODE    I*4	  O        INSTRUMENT MODE AS TAKEN FROM
C					   PARAMETER FILES 
C	RLAT_MIN      R*4	  O	   MINIMUM LATITUDE OF DATA AS
C					   READ FROM 3AL FILES. 3B LIN MAY
C					   BE DIFFERENT.
C	RLAT_MAX      R*4	  O        MAXIMUM LATITUDE OF DATA AS
C					   READ FROM 3AL FILES. 3B MAX 
C					   MAY BE DIFFERENT
C	IPR_MIN	      I*4	  O	   MINIMUM PRESSURE INDEX AS READ
C					   FROM 3AL FILE. 3B MIN MAY BE 
C					   DIFFERENT
C	IPR_MAX       I*4	  O        MAXIMUM PRESSURE INDEX AS READ
C					   FROM 3AL FILE. 3B MAX MAY BE
C					   DIFFERENT.
C	SFDU_BYTE     BYTE(86)    O  	   BYTE ARRAY CONTAINING SFDU
C					   INFORMATION. 
C	ITWO	      I*4	  O	   VALUE OF 2. DENOTES TWO_SIDED
C					   FILE.
C					   VALUE OF 0 MEANS NON-RECURSIVE
C					   TWO-SIDED FILE.
C	DY_ALT	      R*4         O        TBD
C	RNT_ALT       R*4	  O	   TBD
C	IDYLAG	      I*4	  O        NUMBER OF DAYS OF DATA
C					   USED IN BACKWARD RECURSION.
C					   IF NON-RECURSIVE (ITWO_ZERO
C					   = 0) THEN THIS  NUMBER/2 IS
C					   NUMBER OF
C					   DAYS FROM BEGINNING DAY TO
C					   DAY OF MIDDLE OF DATA.
C	IGAP_SIZE     I*4	  O	   MAXIMUM DATA DROP-OUT IN
C					   UNITS OF DAYS. IF GAP IS
C					   LARGER THAN THIS, THEN INITIAL
C					   ESITMATES ARE CALCULATED
C					   ANEW.
C	NPTS_LSQ_MIN  I*4	  O	   MINIMUM NUMBER OF POINTS
C				  	   NEEDED TO INITIALIZE
C					   BACKWARD ESTIMATION.
C	NPTS_INIT     I*4	  O	   MINIMUM NUMBER OF POINTS
C					   NEEDED TO MAKE INITIAL
C					   FORWARD ESTIMATION.
C	PARAMS_20     C*12	  O        DETERMINES IF (FIRST CHARACTER)
C				           3A PARAMETER FILES ARE TO BE
C					   USED, AND IF SO IF ADDITONAL
C					   INFORMATION USED.
C	IVAR_SEL      I*4	  O	   0:ASSUME L3AL DATA GIVES
C					   STANDARD DEVIATION. DATA
C					   WITH NEGATIVE VALUES OF
C					   STANDARD DEVIATION ARE
C					   REJECTED.
C					   -1:ASSUME L3AL DATA GIVES
C					   STANDARD DEVIATION. DATA
C					   WITH NEGATIVE VALUES OF
C					   STANDARD DEVIATION ARE
C					   ACCEPTED, AND THE POSITIVE
C					   VALUE OF THE STANDARD DEVIATION
C					   ARE USED.
C				           1:ASSUME L3AL DATA GIVES
C					   VARIANCE. DATA WITH NEGATIVE
C					   VALUES OF VARIANCE
C					   ARE REJECTED.
C				           -2:ASSUME L3AL DATA GIVES
C					   VARIANCE. DATA WITH NEGATIVE
C					   VALUES OF VARIANCE ARE ACCEPTED
C					   AND THE POSITIVE VALUE OF THE
C					   VARIANCE IS USED.
C	STDFAC	      R*4	  O	   IF GT 0, THEN L3AL DATA IS FILTERED
C					   BASED ON STDFAC*DATA_SIGMA 
C					   DIFFERENCE. DATA_SIGMA IS THE
C					   STANDARD DEVIATION AS FOUND FROM
C					   THE DATA ITSELF, ASSUMING A 
C					   STAIONARY STATE.
C					   ONLY THE CHARACTERISTIC OF THE
C					   VALUE IS USED FOR THIS PURPOSE.
C					   THE MANTISSA IS USED TO SET A
C					   MINIMUM FOR THE VALUE OF DATA_SIGMA.
C					   IF THE VALUE IS NEGATIVE, THEN
C					   THE POSITIVE VALUE IS USED AS 
C					   DESCRIBED ABOVE, BUT NEGATIVE
C					   VALUES OF DATA ARE REJECTED.
C	IHDR_FLAG     I*4	  O	   TBD 
C	IWR_TIME      I*4(2)	  O	   TIME IN UDTF WHEN FILE
C					   WAS CREATED.
C	VER_NUM	      R*4	  O	   L3B SOFTWARE VERSION NUMBER.
C	NPTS_BKWRD_MIN
C		      I*4	  O	   THE TWO RIGHT-MOST DIGITS DENOTE THE
C					   MINIMUM NUMBER OF DATA POINTS
C					   REQUIRED IN BACKWARD RECURSION.
C					   IF THE NUMBER OF DATA POINTS
C					   WAS LESS THAN THIS, NO TWO-SIDED
C					   COEFFICIENTS WERE GENERATED.
C					   THE SECOND AND THIRD (COUNTING FROM
C					   ZERO) RIGHT-MOST DIGITS DENOTE
C					   THE POWER (BASE 10) OF THE MINIMUM
C					   VALUE OF THE STANDARD DEVIATION
C					   OF THE 3AL DATA.
C					   THE FOURTH AND FIFTH (COUNTING FROM
C					   ZERO) RIGHT-MOST DIGITS DENOTE
C					   THE POWER (BASE 10) OF THE MINIMUM
C					   VALUE OF THE 3AL DATA. IF
C					   NPTS_BKWRD_MIN BEGINS WITH A MINUS,
C					   THE EXPONENTS ARE NEGATIVE.
C					   NOTE THAT IF THE EXPONENT IS ZERO,
C					   THEN A MINIMUM VALUE TEST IS NOT
C					   MADE, AND ALL VALUES ARE ACCEPTED.  
C					   FOR EXAMPLE, A CHARACTER STRING
C					   OF -000131800 MEANS THAT THE
C					   MINIMUM NUMBER OF POINTS REQUIRED
C					   IS 0, THE MINIMUM VALUE OF STANDARD
C					   DEVIATION ACCEPTED IS 
C					   10 ** (-18), AND THE
C					   MINIMUM DATA VALUE ACCEPTED IS
C					   10 ** (-13).
C	NPTS_BKWRD_DAY_MIN
C		      I*4	  O	   MINIMUM NUMBER OF DATA POINTS
C					   IN LAST HALF OF TWO-SIDED DAY
C					   REQUIRED. IF THE NUMBER OF DATA
C					   POINTS WAS LESS THAN THIS, THEN
C					   THE DATA FOR THE TW-SIDED DAY
C					   WAS DISCARDED, AND PROCESSING
C					   CONTINUED AS IF THERE WAS NO DATA
C					   ON THAT DAY.
C	IOS	      I*4	  O        FORTRAN STATUS CODE
C-----------------------------------------------------------------------	
	INTEGER*4 INTIME(2),IWR_TIME(2),IREC_TYPE/0/
        INTEGER*4 ISTRT(2),ISTOP(2)
	CHARACTER*(*) INSTR,PARAM
	BYTE SFDU_BYTE(86)
	CHARACTER*12 PARAMS_20
C	CHARACTER*86 SFDU_CHAR
C	EQUIVALENCE (SFDU_CHAR,SFDU_BYTE)
C
C	READ SFDU INFORMATION
C
 	ILTDY1 = 1
	READ(LUN,ERR=200,IOSTAT=IOS) ILTDY,SFDU_BYTE
C	PRINT *, '---------------------------'
C	PRINT *, 'SFDU RECORD'
C	PRINT *, '---------------------------'
C       WRITE(6,*) 'ILTDY : ',ILTDY
C	WRITE(6,100) SFDU_BYTE
C 100	FORMAT(50A1)
C
 	ILTDY1 = 2
C
C	NEW NEWER LONGER RECORDS (VERSION 3.06 OR LATER)
C	IF ERROR 67 THEN DATA OVERRUN. TRY SHORTER OLDER RECORDS
C
	IF(IREC_TYPE.EQ.0) THEN
	READ(LUN,ERR=200,IOSTAT=IOS)
     &	ILTDY,RLAT_MIN,RLAT_MAX,IPR_MIN,IPR_MAX,
     &	IORB_MODE,INSTR_MODE,INTIME(1),INTIME(2),ITWO,
     &	(INSTR(L:L),L=1,6),(PARAM(L:L),L=1,6),
     &	DY_ALT,RNT_ALT,
     &	IDYLAG,IGAP_SIZE,NPTS_LSQ_MIN,NPTS_INIT,
     &	PARAMS_20,IVAR_SEL,STDFAC,IHDR_FLAG,
     &	IWR_TIME(1),IWR_TIME(2),VER_NUM,
     &	NPTS_BKWRD_MIN,NPTS_BKWRD_DAY_MIN, 
     &	ISTRT(1),ISTRT(2),ISTOP(1),ISTOP(2)
C	PRINT *, ''
C	PRINT *, '---------------------------'
C	PRINT *, 'LABEL RECORD'
C	PRINT *, '---------------------------'
C       WRITE(6,*) 'ILTDY : ',ILTDY
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
	READ(LUN,ERR=300,IOSTAT=IOS)
     &	ILTDY,RLAT_MIN,RLAT_MAX,IPR_MIN,IPR_MAX,
     &	IORB_MODE,INSTR_MODE,INTIME(1),INTIME(2),ITWO,
     &	(INSTR(L:L),L=1,6),(PARAM(L:L),L=1,6),
     &	DY_ALT,RNT_ALT,
     &	IDYLAG,IGAP_SIZE,NPTS_LSQ_MIN,NPTS_INIT,
     &	PARAMS_20,IVAR_SEL,STDFAC,IHDR_FLAG
	ENDIF
  300	CONTINUE
	RETURN
	END

	SUBROUTINE RD3B_TWO(LUN,RLAT,IPR,COEF,NCOEF,
     &	NPTSF,NPTSB,NPTS_DAY,NPTS_BAD,
     &	TAUINV,INTIME1,INTIME2,
     &	XLONG,TAVE,AVE_LST,AVE_SZA,
     &  RES_DEV1,STD_CALC,NPTS_DEV,
     &	SEC_FST,IFLAG_FW,IFLAG,PARAM_VALUES,ISTD_FLAG,
     &	IVER,ICYC,
     &	IGMT,STD_OLD,IOS)
C
C
C
C	THIS ROUTINE IS PART OF THE SOFTWARE TO GENERATE LEVEL 3B
C	FOR UARS.  IT READS THE TWO_SIDED FOURIER COEFFICIENTS,
C	THE CORRESPONDING COVARIANCE MATRIX AND OTHE CONTROL
C	PARAMETERS WHICH WERE USED.
C
C	CONTACT F. T. HUANG/SSAI (301) 731-9300
C	AUGUST, 1992
C-----------------------------------------------------------------
C	ARGUMENT      TYPE        I/O      DESCRIPTION
C       ------------------------------------------------------------
C
C	LUN	      I*4	  I	   LOGICAL UNIT OF FILE
C	RLAT	      R*4	  I/O	   REQUESTED LATITUDE/RETURNED 
C					   LATITUDE (DEGREES).
C	IPR	      I*4	  I/O	   REQUESTED PRESSURE LEVEL/
C					   RETURNED PRESSURE INDEX
C	COEF	      R*4 VECTOR  O        FOURIER COEFFICIENTS.
C		      (NCOEF)
C	NCOEF         I*4	  O	   NUMBER OF FOURIER COEFFICIENTS.
C	NPTSF	      I*4	  O	   NUMBER OF DATA POINTS
C					   FOR NEXT 24 HRS FROM 12GMT.
C	NPTSB	      I*4	  O	   NUMBER OF DATA POINTS
C					   FOR NEXT 24 HRS FROM 12GMT.
C	NPTS_DAY      I*4	  O	   NUMBER OF DATA POINTS
C					   FOR + - 12 HRS FROM 12GMT.
C	NPTS_BAD      I*4	  O        NUMBER OF DATA POINTS FILTERED
C					   (IN BACKWARD RECURSION)
C					   OUT WHEN STDFAC IS NOT ZERO.
C	TAUINV	      R*4	  O	   PARAMETER DENOTING
C					   FILTER MEMORY LENGTH.
C	INTIME1	      I*4(2)	  O	   UDTF TIME OF LAST DATA POINT USED
C					   IN FORWARD RECURSION.
C	INTIME2	      I*4(2)	  O	   UDTF TIME OF EARLIEST (TIMEWISE)
C					   DATA USED IN BACKWARD RECURSION.
C	XLONG	      R*4         O	   LONGITUDE VALUES OF 3AL DATA
C		      (NPTS_DAY) 	   IN CURRENT DAY.
C	TAVE	      R*4	  O	   AVERAGE OF DATA POINTS FOR THAT
C					   DAY
C	AVE_LST	      R*4	  O	   AVERAGE LOCAL SOLAR TIME AT
C					   LATITUDE
C	AVE_SZA	      R*4	  O	   AVERAGE SOLAR ZENITH ANGLE AT
C					   LATITUDE
C
C	RES_DEV1      R*4	  O	   CUMMULATIVE RESIDUAL DEVIATION
C					   AT PREVIOUS DATA POINT (FROM 12GMT) 
C					   BETWEEN FIT AND DATA, NORMALIZED
C					   TO 3AL DATA STANDARD DEVIATION.
C	STD_CALC      R*4	  O	   STANDARD DEVIATION CALCULATED FROM
C					   3AL DATA ASSUMING STATIONARY
C					   TIME SERIES, OR FROM A PRIORI
C					   FOURIER COEFFICIENTS (IF AVAILABLE).
C	NPTS_DEV      I*4	  O        CUMMULATIVE NUMBER POINTS USED
C					   IN FINDING RES_DEV1.
C	SEC_FST	      R*4	  N/A	   TIME IN SECONDS OF FIRST DATA POINT.
C	IFLAG_FW      I*4	  O        STATUS FLAG FROM FORWARD FILE.
C	IFLAG	      I*4	  O        STATUS FLAG: 0:NOMINAL
C					   12: NUMBER OF DATA IN NEXT 36HRS
C					       FROM TWO_SIDED OUTPUT TIME
C					       LESS THAN 10.
C					   13: TOTAL NUMBER OF POINTS IN
C					       NEXT IDYLAG DAYS LESS THAN 30.
C					   14: BOTH 12 AND 13 ARE TRUE.
C					   100 IS ADDED TO IFLAG WHEN
C					   THERE ALSO ARE NOT ENOUGH
C					   DATA POINTS TO BEGIN BACKWARD
C					   RECURSION AND FORWARD COEFFS
C					   ARE USED TO START BACKWARD
C					   RECURSION.
C	PARAM_VALUES  R*4(NPTS_DAY)   O	   3AL 3AL DATA VALUES
C	ISTD_FLAG     I*4	  O	   A VALUE OF 1 MEANS THAT ONE
C					   OR MORE NEGATIVE STANDARD
C					   DEVIATIONS OCCURED IN THE
C					   3AL DATA. IN THE MLS CASE,
C					   THIS MEANS THAT THE DATA VALUE
C					   DEPENDS TO A SIGNIFICANT
C					   DEGREE ON THE INITIAL OR
C					   A PRIORI VALUE.
C	IVER	      I*4         O        3AL VERSION NUMBER
C	ICYC	      I*4	  O	   3AL CYCLE NUMBER
C	IGMT	      I*4         O	   TBD
C	STD_OLD	      R*4	  O	   TBD
C	IOS	      I*4	  O	   FORTRAN ERROR STATUS.
C-----------------------------------------------------------------------	
	DIMENSION COEF(1),XLONG(1),PARAM_VALUES(100),STD_OLD(200)
	CHARACTER*6 INSTR,PARAM
	CHARACTER*86 SFDU_CHAR
	INTEGER*4 INTIME1(2),INTIME2(2),IGMT(200)
	EQUIVALENCE (SFDU_CHAR,SFDU_BYTE)
	INTEGER*4 IREC_TYPE/0/
C
C	READ DATA
C
	ILAT =  1 + (88.0 - RLAT)/4.0 
	ILTDY1 = ILAT * 50 + IPR
C
C	NEW NEWER LONGER RECORDS (VERSION 3.06 OR LATER)
C	IF ERROR 67 THEN DATA OVERRUN. TRY SHORTER OLDER RECORDS
C
	IF(IREC_TYPE.EQ.0) THEN
	READ (LUN,ERR=200,IOSTAT=IOS)
     &	ILTDY,IPR,NCOEF,
     &	RLAT,NPTSF,NPTSB,NPTS_DAY,NPTS_BAD,
     &	TAUINV,INTIME1(1),INTIME1(2),INTIME2(1),INTIME2(2),
     &	TAVE,AVE_LST,AVE_SZA,
     &	RES_DEV1,STD_CALC,NPTS_DEV,
     &	SEC_FST,IFLAG_FW,IFLAG,
     &	(COEF(I),I=1,NCOEF),
     &	(XLONG(I),I=1,NPTS_DAY), 
     &	(PARAM_VALUES(I),I=1,NPTS_DAY),
     &	ISTD_FLAG,IVER,ICYC, 
C    &	(IGMT(I),I=1,NPTS_DAY),
     &	(STD_OLD(I),I=1,NPTS_DAY)
C	WRITE(6,*) 'ILTDY : ',ILTDY
  200	CONTINUE
	ENDIF
C
C	TRY READING OLDER RECORD TYPE IF ERROR IS 67
C
	IF(IOS.EQ.67.OR.IREC_TYPE.EQ.1) THEN
C
	IVER = 0
	ICYC = 0
C
	IREC_TYPE = 1
	READ (LUN,ERR=300,IOSTAT=IOS)
     &	ILTDY,IPR,
     &	NCOEF,
     &	RLAT,NPTSF,NPTSB,NPTS_DAY,IDYLAG,
     &	TAUINV,INTIME1(1),INTIME1(2),INTIME2(1),INTIME2(2),
     &	TAVE,AVE_LST,AVE_SZA,
     &	RES_DEV1,RES_DEV2,NPTS_DEV,
     &	SIGMA,IFLAG_FW,IFLAG,
     &	(COEF(I),I=1,NCOEF),
     &	(XLONG(I),I=1,NPTS_DAY),
     &	(PARAM_VALUES(I),I=1,NPTS_DAY)
	ENDIF
  300	CONTINUE
C       WRITE(6,*) 'ERR IN OPEN',IOS
C
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

