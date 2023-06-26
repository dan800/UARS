/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   Program name: read_l3b_cov.c

   Purpose:      read uars level 3b corvariance prod files.

   Author:       Ke - Jun Sun   Hughes STX

   E-mail:       sun@daac.gsfc.nasa.gov

   Usage:        read_l3b_cov

   Description:  This program allows a user to read the Fourier covariances
                 stored in the UARS level 3B data files, and transform these
                 to geophysically meaningful values defined on 4 degree
                 latitude intervals. The program allows a user to input a
                 UARS level 3B covariance data filename in the current
                 directory, after which he has the option of dumping the
                 results to the screen or saving the results in an ascii disk
                 file. Options are provided for allowing selection of the
                 latitude range of interest, and minimum and maximum pressure
                 levels.

   to compile : cc -o read_l3b_cov read_l3b_cov.c -lm 

 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*   4 bytes recl (say N) | N bytes | 4 bytes recl  */

#include <stdio.h>
#include <string.h>
#include <math.h>
#define MAX 20

void check_int_flt(char numbstr[MAX], int *outint, float *outfloat, int flag);
void error_cv(float covar[13][13],int ncof21,int ntheta,float theta[360],
              float z[360], float zmin);
void jul(int yr, int juln, int *month, int *day);
main()
{
 char line[1000],sfdu[87];
 char temp[30],parameter[9];
 char instmt[7],*unit,*alt_prs,*mb_km;
 int *int4,recl,coef,pts_day,ntheta,addr,count,recs,ret,i,j,k,resolution;
 float *real4;
 float coefs[13][13];
 float theta[360],z[360];
 float maxlat,minlat,req_maxlat,req_minlat,rec_lat;
 float rec_prs,maxprs,minprs,req_maxprs,req_minprs,dummyflt,zmin;
 int req_p_level1,req_p_level2;
 int rec_p_level,min_p_level,max_p_level;
 char outfile[80],defaultfile[80],infile[80]; 
 FILE *ifp,*ptr,*outfp;
 char tmpstr[50],numbstr[20];
 char valid,svfile,hitmiss;
 int  dummyint,m,n,multi,yr,d,jday,ds_as;
 static char *months[12] = {"JAN","FEB","MAR","APR","MAY",
                  "JUN","JUL","AUG","SEP","OCT","NOV","DEC"};


 printf("Enter UARS L3B Covariance file name :\n");
 scanf("%s",infile);getchar();
 ifp=fopen(infile,"r");
 if (ifp==NULL)
  {
   printf(" Can not open data file %s\n",infile);
   exit(0);
  }

/* *********** Scan thru whole file to get min/max lat and prs ************ */

 fread(line,1,4,ifp);
 int4=(int *)&line[0];
 recl=*int4;
 fread(line,1,recl,ifp);
 fread(line,1,4,ifp);
 int4=(int *)&line[0];
 recl=*int4;

 fread(line,1,4,ifp);
 int4=(int *)&line[0];
 recl=*int4;

 if (recl != 48)
  {
   printf("This is not a UARS L3B COVARIANCE file.\n");
   printf("Please check file name and try again.\n");
   exit(1);
  }

 fread(line,1,recl,ifp);
 fread(line,1,4,ifp);
 int4=(int *)&line[0];
 recl=*int4;

 max_p_level=0; min_p_level=100;
 maxlat=-90.0;  minlat=90.0;

 do
 {
  ret =  fread(line,1,4,ifp);
  if (ret > 0)
   { 
    recs++;
    int4=(int *)&line[0];
    recl=*int4;
/*    printf("\n\n ******** DATA RECORD %d  *******\n",recs); */
    fread(line,1,recl,ifp);
    int4=(int *)&line[4];
    rec_p_level=*int4;
    real4=(float *)&line[12];
    rec_lat=*real4;
    if(rec_p_level > max_p_level)
      max_p_level=rec_p_level;
     else
      if (rec_p_level < min_p_level)min_p_level=rec_p_level;

    if(rec_lat > maxlat)maxlat=rec_lat;
     else
      if (rec_lat < minlat)minlat=rec_lat;

    fread(line,1,4,ifp);
    int4=(int *)&line[0];
    recl=*int4;
    }     /*  endif ret > 0 */
   } while (ret > 0);

 rewind(ifp);

 maxprs = 1000.0 * pow(10.0,(min_p_level/(-6.0)));
 minprs = 1000.0 * pow(10.0,(max_p_level/(-6.0)));

 printf("Do you wish to save result to an ASCII file (y/n) ?\n");
 svfile=getchar();getchar();
 if(svfile=='Y' || svfile=='y')
  {
   strncpy(defaultfile,infile,strlen(infile));
   defaultfile[strlen(infile)]='\0';
   strcat(defaultfile,".DUMP");
   printf("\nEnter output file name, or <CR> for default name :\n");
   printf(" %s\n",defaultfile);
   gets(outfile);
   if (strlen(outfile)==0)
     strcpy(outfile,defaultfile);
   if ((outfp = fopen(outfile,"w")) ==NULL)
    {
     printf(" cannot open file %s \n",outfile);
     exit(1);
    }
   else
    printf("The result will be stored in %s\n",outfile);
  }

/* ********** Now read it 2nd time **************** */

 fread(line,1,4,ifp);
 int4=(int *)&line[0];
 recl=*int4;
 fread(line,1,recl,ifp);

/*
 int4=(int *)&line[0];
 printf(" LTDY : %d\n",*int4);
 strncpy(sfdu,&line[4],54);
 printf(" %s\n",sfdu);
 strncpy(sfdu,&line[58],34);
 printf(" %s\n",sfdu);
*/

 fread(line,1,4,ifp);
 int4=(int *)&line[0];
 recl=*int4;
 fread(line,1,4,ifp);
 int4=(int *)&line[0];
 recl=*int4;
/*  printf(" ******** DATA HDR RECORD ********\n");  */

 fread(line,1,recl,ifp);

/* <- remove this pair if you want to see what's in header

 int4=(int *)&line[0];
 printf(" LTDY      : %d\n",*int4);
*/
 int4=(int *)&line[4];
 ds_as=*int4;
/*
 printf(" ORB_MODE  : %d\n",*int4);
 int4=(int *)&line[8];
 printf(" INSTR_MODE: %d\n",*int4);
 int4=(int *)&line[12];
 printf(" TWO       : %d\n",*int4);
*/
 int4=(int *)&line[16];
 jday=(*int4)%1000;
 yr=(*int4)/1000;
 yr += 1900;
 jul(yr,jday,&m,&d);
/*
 printf(" INTIME(1) : %d\n",*int4);
 int4=(int *)&line[20];
 printf(" INTIME(2) : %d\n",*int4);
 strncpy(temp,&line[24],6);
 printf(" INSTR     : %s\n",temp);
 strncpy(temp,&line[30],6);
 printf(" PARAM     : %s\n",temp);
 int4=(int *)&line[36];
 printf(" WRTIME(1) : %d\n",*int4);
 int4=(int *)&line[40];
 printf(" WRTIME(2) : %d\n",*int4);
 real4=(float *)&line[44];
 printf(" VER_NUM   : %f\n",*real4);

 ->  */

 strncpy(temp,&line[24],6);
 sscanf(temp,"%s",instmt);

 strncpy(temp,&line[30],6);
 sscanf(temp,"%s",parameter);

 if(isdigit(parameter[0]))
 {
  strcpy(temp,"AEROSOL");
  strcat(temp,parameter);
  sscanf(temp,"%s",parameter);
 }
 
 multi = 1;
 if (strcmp(instmt,"HRDI")==0||
    strcmp(instmt,"WINDII")==0)
    unit = "(m/s)";
  else if (strcmp(parameter,"TEMP")==0)
    unit = "(degrees K)";
  else if (strcmp(parameter,"CF2CL2")==0||strcmp(parameter,"CFCL3")==0
         ||strcmp(parameter,"CH4")==0||strcmp(parameter,"CO")==0
         ||strcmp(parameter,"CLONO2")==0||strcmp(parameter,"H2O")==0
         ||strcmp(parameter,"HNO3")==0||strcmp(parameter,"N2O5")==0
         ||strcmp(parameter,"N2O")==0||strcmp(parameter,"NO2")==0
         ||strcmp(parameter,"NO")==0||strcmp(parameter,"HCL")==0
         ||strcmp(parameter,"HF")==0||strcmp(parameter,"CLO")==0
         ||strcmp(parameter,"O3_183")==0||strcmp(parameter,"O3_205")==0
         ||strcmp(parameter,"O3B8")==0||strcmp(parameter,"O3B9")==0
         ||strcmp(parameter,"O3")==0)
   {
    unit = "(ppmv)";
    multi = 1000000;
   }
  else
    unit = "(1/km)";

 fread(line,1,4,ifp);
 int4=(int *)&line[0];
 recl=*int4;

 if((strcmp(instmt,"HRDI")==0 || strcmp(instmt,"WINDII")==0) &&
    (strcmp(parameter,"MWIN_A")==0 || strcmp(parameter,"ZWIN_A")==0)||
    (strcmp(instmt,"WINDII")==0 && strcmp(parameter,"TEMP")==0))
   {alt_prs="altitude";mb_km="km";}
 else
   {alt_prs="pressure";mb_km="mb";}

 printf("****************************************************\n");
 printf(" Instrument    : %s\n",instmt);
 printf(" Parameter     : %s\n",parameter);
 printf(" Date          : %02d-%3s-%ld \n",d,months[m-1],yr);
 printf(" Time          : 12:00:00.00 GMT\n");
 if(ds_as == 0)
 printf(" Orbit Node    : Ascending\n");
 else
 printf(" Orbit Node    : Descending\n");
 printf(" Min. latitude : %6.2f deg\n",minlat);
 printf(" Max. latitude : %6.2f deg\n",maxlat);
 printf(" Min. %8s : %-9.4#G %2s\n",alt_prs,minprs,mb_km);
 printf(" Max. %8s : %-9.4#G %2s\n",alt_prs,maxprs,mb_km);
 printf("****************************************************\n");

 hitmiss='n';

   do
   {
    valid = 'y';
    printf("\nEnter Northern most latitude (degrees) :\n");
/*    scanf("%f",&req_maxlat);  */
    fgets(numbstr,MAX,stdin);
    check_int_flt(numbstr,&dummyint,&req_maxlat,2);
    if( req_maxlat > maxlat || req_maxlat < minlat)
    {
     printf(" Invalid latitude entered, try again. \n");
     valid = 'n';
    }
   } while(valid=='n');

   if((int)(ceil(req_maxlat)) %4 !=0)
    {
     while((int)(ceil(req_maxlat)) % 4 != 0 && req_maxlat <= maxlat - 1)
      req_maxlat++;
     req_maxlat = ceil(req_maxlat);
     printf("Using next northernmost latitude at %5.1f deg.\n",req_maxlat);
    }

   do
   {
    valid = 'y';
    printf("Enter Southern most latitude (degrees) :\n");
/*    scanf("%f",&req_minlat);  */
    fgets(numbstr,MAX,stdin);
    check_int_flt(numbstr,&dummyint,&req_minlat,2);
    if ( req_minlat > maxlat || req_minlat < minlat)
    {
     printf(" Invalid latitude entered, try again. \n");
     valid = 'n';
    }
    if ( req_maxlat < req_minlat )
    {
     printf(" Northern most latitude is less than Southern most latitude \n");
     printf(" Try again. \n");
     valid = 'n';
    }

   } while(valid=='n');

   if((int)(floor(req_minlat)) %4 != 0)
    {
     while((int)(floor(req_minlat)) % 4 != 0 && req_minlat >= minlat + 1)
      req_minlat--;
     req_minlat = floor(req_minlat);
     printf("Using next southernmost latitude at %5.1f deg.\n",req_minlat);
    }

   do
   {
    valid = 'y';
   printf("Enter maximum %8s level (%2s) : \n",alt_prs,mb_km);
    fgets(numbstr,MAX,stdin);
    check_int_flt(numbstr,&dummyint,&req_maxprs,2);
    if( req_maxprs > maxprs || req_maxprs < minprs)
    {
     printf(" Invalid %8s entered, try again. \n",alt_prs);
     valid = 'n';
    }
   } while(valid=='n');
   if(req_maxprs!=1000.0 && req_maxprs!=100.0 &&
      req_maxprs!=10.0 && req_maxprs!=1.0)
    {
    req_p_level1 = (int)(floor((-6)*log10(req_maxprs/1000.0)));
    req_maxprs = 1000.0*pow(10.0,(req_p_level1)/(-6.0));
    printf("Using next higher %8s level at %10.3#G %2s.\n",alt_prs,req_maxprs,
           mb_km);
    }

   do
   {
    valid = 'y';
    printf("Enter minimum %8s level (%2s) : \n",alt_prs,mb_km);
    fgets(numbstr,MAX,stdin);
    check_int_flt(numbstr,&dummyint,&req_minprs,2);
    if ( req_minprs > maxprs || req_minprs < minprs)
    {
     printf(" Invalid %8s entered, try again. \n",alt_prs);
     valid = 'n';
    }
    if ( req_maxprs < req_minprs )
    {
     printf(" Maximum %8s less than Minimum %8s,\n",alt_prs,alt_prs);
     printf(" Try again. \n");
     valid = 'n';
    }

   } while(valid=='n');

   if(req_minprs!=1000.0 && req_minprs!=100.0 &&
      req_minprs!=10.0 && req_minprs!=1.0)
    {
      req_p_level2 = (int)(ceil((-6)*log10(req_minprs/1000.0)));
      req_minprs = 1000.0*pow(10.0,(req_p_level2)/(-6.0));
      printf("Using next lower pressure level at %10.3#G mb.\n",req_minprs);
    }

 printf(" Data will start at -180.0 degrees longitude.\n");
 do
 { valid = 'y';
   printf(" Enter longitudinal resolution between 1 and 5 degrees.\n");
   fgets(numbstr,MAX,stdin);
   check_int_flt(numbstr,&resolution,&dummyflt,1);
   if (resolution < 1 || resolution > 5)
    {
     printf(" Invalid number, try again \n");
     valid = 'n';
    }
 }while (valid == 'n');

 ntheta = 360/resolution;
 for(k=0;k<ntheta;k++)
   theta[k] = resolution*(k) - 180.0;

 if(svfile=='Y' || svfile=='y')
   fprintf(outfp,"\nFile : %s",infile);
 else
   printf("\nFile : %s",infile);

/*          Read in record after record        */

/*  remove those comment characters if you want to see all vars */

 do
 {
  ret =  fread(line,1,4,ifp);
  if (ret > 0)
   {
    recs++;
    int4=(int *)&line[0];
    recl=*int4;
/*    printf("\n\n ******** DATA RECORD %d  *******\n",recs);  */

    fread(line,1,recl,ifp);
/*
    int4=(int *)&line[0];
    printf(" LTDY        : %d\n",*int4);
*/
    int4=(int *)&line[4];
    rec_p_level=*int4;
    rec_prs=1000.0*pow(10.0,(rec_p_level)/(-6.0));

/*    printf(" IPR         : %d\n",rec_p_level);  */

    int4=(int *)&line[8];
    coef = *int4;
/*    printf(" COEF        : %d\n",coef);  */  

    real4=(float *)&line[12];
    rec_lat=*real4;
/*
    printf(" RLAT        : %f\n",*real4);
    int4=(int *)&line[16];
    printf(" INTIME1(1)  : %d\n",*int4);
    int4=(int *)&line[20];
    printf(" INTIME1(2)  : %d\n",*int4);
    int4=(int *)&line[24];
    printf(" INTIME2(1)  : %d\n",*int4);
    int4=(int *)&line[28];
    printf(" INTIME2(2)  : %d\n",*int4);
*/
   if(rec_lat >= req_minlat && rec_lat <= req_maxlat &&
      rec_prs >= req_minprs && rec_prs <= req_maxprs)
    {
      hitmiss='y';
      addr=32;
      for(n=0;n<13;n++)
       for(m=0;m<13;m++)
       {
        real4=(float *)&line[addr];
        coefs[n][m]=*real4;
        addr += 4;
       }
      error_cv(coefs,7,ntheta,theta,z,zmin);

      if (svfile=='y')
       {
        fprintf(outfp,"\n\n********************************************************************************\n");
        fprintf(outfp," %-6s %-12s at Latitude : %5.1f and %8s level : %10.3#G %2s\n",parameter,unit,rec_lat,alt_prs,rec_prs,mb_km);
        fprintf(outfp,"\n%3d data values beginning at longitude -180.0 incremented every %1d degrees :\n\n",ntheta,resolution);

      for (count=0;count<ntheta;count++)
       {
         fprintf(outfp," %12.5e",z[count]*multi);
         if ((count+1)%6==0) fprintf(outfp,"\n");
       }
       fprintf(outfp,"\n");
       }
      else
      {
        printf("\n\n********************************************************************************\n");
        printf(" %-6s %-12s at Latitude : %5.1f and %8s level : %10.3#G %2s\n" ,
        parameter,unit,rec_lat,alt_prs,rec_prs,mb_km);
        printf("\n%3d data values beginning at longitude -180.0 incremented every %1d degrees :\n\n",ntheta,resolution);

      for (count=0;count<ntheta;count++)
       {
         printf(" %12.5e",z[count]*multi);
         if ((count+1)%6==0) printf("\n");
       }
       printf("\n");
      }
    }
    fread(line,1,4,ifp);
    int4=(int *)&line[0];
    recl=*int4;

    }     /*  endif ret > 0 */
   } while (ret > 0);

 if(hitmiss == 'n' && svfile=='n')
  printf("There is no data in the region and %8s level you selected\n",alt_prs);


 if(svfile == 'y')
  {
   if(hitmiss == 'y')
    printf("Total %ld bytes data are saved in file %s\n",ftell(outfp),outfile);
   else
    printf("No data is being saved in ascii file %s\n",outfile);

   fclose(outfp);
  }

 fclose(ifp);

}
/* ************************************************************************

 Function check_int_flt is to convert a input string to a integer or
 a float number.

 ************************************************************************** */

void check_int_flt(char numbstr[MAX], int *outint, float *outfloat, int flag)

{
 int ret;
 do
 {
  if(flag == 1)
   {
    ret=sscanf(numbstr,"%d",outint);
    if (ret != 1)
     {
      printf("Try again, please enter an integer ==>");
      fgets(numbstr,MAX,stdin);
      printf("\n");
     }
   }
  else
   {
    ret=sscanf(numbstr,"%f",outfloat);
    if (ret != 1)
     {printf("Try again, please enter a real number ==>");
      fgets(numbstr,MAX,stdin);
      printf("\n");
     }
   }
 }while(ret != 1);
}
/*
        SUBROUTINE ERROR_CV(COVAR,NCOF21,NTHETA,THETA,Z,ZMEAN)
C
C       THIS ROUTINE CALCULATES ERRORS (FROM VARIANCES) FROM THE INPUT
C       LEVEL 3B COVARIANCE MATRIX COVAR AT THE ANGLES THETA IN DEGREES.
C
C       NCOF21 IS SUCH  THAT 2*NCOF21 + 1 ARE THE TOTAL NUMBER
C       OF COEFFICIENTS.
C
C       THE FIRST NCOF21 ARE THE CONSTANT AND COSINE COFS
C       AND THE LAST NCOF21 ARE THE SINE COFS.
C
C       FOR THE UARS CASE, THERE NCOF21 IS 6, AND THE FIRST T COEFFS
C       ARE THE CONSTANT AND SIX COSINE COEFFS.  THE LAST 6 ARE THE
C       CORRESPONDING SINE COEFFS.
C
C       THETA(I) IS THE ARRAY OF ANGLES (FOR UARS THEY ARE
C       LONGITUDES) IN DEGREES, NTHETA IS THE NUMBER OF
C
C       ARGUMENT        I/O     TYPE            DESCRIPTION
C       --------        ---     ----            -----------------------
C       COVAR           I       R*4 ARRAY       COVARIANCE ARRAY.
C                               (2*NCOF21+1)x   FROM THE DIAGONAL TERMS
C                               (2*NCOF21+1)    ONE CAN GET THE MEAN ERROR.
C
C       NCOF21          I       I*4             DENOTES TOTAL NUMBER OF
C                                               COSINE TERMS, INCUDING
C                                               CONSTANT TERM.
C
C       NTHETA          I       I*4             NUMBER OF ANGLES
C
C       THETA           I       R*4(NTHETA)     ANGLES IN DEGREES WHERE
C                                               SERIES IS TO BE EVALUATED.
C
C       Z               O       R*4(NTHETA)     RESULTS OF ERRORS FOR EACH
C                                               NTHETA ANGLES.
C
C       BY JAMES JOHNSON, HUGHES STX (1996)
C       MODIFIED FROM SUBROUTINE FOURIER_SUM
C-------------------------------------------------------------------------
 */
void error_cv(float covar[13][13], int ncof21, int ntheta, float theta[360],
              float z[360], float zmean)
{
 float pi=acos(-1.0);
 float twpi=2.0*pi;
 int i,j,m,n;
 float arg,arg1,cos1,sin1,sum,cvk[13];

 for(i=0;i<ntheta;i++)
  {
   arg=twpi * theta[i]/360.0;
   sum=0.0;
   cvk[0] = 1.0;
   for(j=1;j<ncof21;j++)
    {
     arg1 = arg * (j);
     cos1 = cos(arg1);
     sin1 = sin(arg1);
     cvk[j] = cos1;
     cvk[j+6] = sin1;
    }
   for(n=0;n<13;n++)
    for(m=0;m<13;m++)
     sum = sum + cvk[m]*covar[n][m]*cvk[n];
   z[i] = sqrt(sum);

   zmean = 0.0;
   for(n=1;n<13;n++)
    zmean = zmean + covar[n][n];
   zmean = sqrt(covar[0][0] + 0.5 * zmean);
  }
}
void jul(yr,juln,month,day)

int yr,juln,*month,*day;
{
 int mntab[12]={31,28,31,30,31,30,31,31,30,31,30,31};
 int sum,mn;
 if ((yr % 4)== 0)
   mntab[1]=29;

  *month = 0;
  *day = 0;
  sum = 0;
  for (mn = 0; mn <= 11; mn++)
   {
    sum += mntab[mn];
    if(juln <= sum)
     {
      * month = mn+1;
      sum -= mntab[mn];
      *day = juln - sum;
      break;
     }
   }
}

