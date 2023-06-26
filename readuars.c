/* **********************************************************************
 *
 * Name: readuars
 *
 * Synopsis:
 *     readuars [-s <slat> <nlat> [<wlon> <elon>]] [-v <zmin> <zmax>] \
 *      	[-u <uparam> ] [-w <start> <end>] <filename>
 *
 * Description:
 *
 *     This program reads and dumps the contents of all UARS level 3A
 *     data files (3AL, 3AT, 3LP, 3TP, 3AS/3BS, and UARS Correlative).
 *
 *     The input file may be either the actual data file (PROD extension),
 *     or the metadata file (META extension).  Output can be redirected to
 *     an ASCII file.  Subsetting by latitude and vertical height can also
 *     be done.
 *
 * Notes:
 *
 *     UARS files from the DAAC are native binary with big-endian ordered data.
 *     Users reading these files on a PC or other litte-endian machines must
 *     use the byteswap option -b, or data values will come out wrong.
 *
 * Compile:
 *
 *     To compile this progam : cc -o readuars readuars.c -lm
 *
 * ^AUTHOR:  James Johnson
 *
 * ^ADDRESS: jjohnson@daac.gsfc.nasa.gov
 *           NASA GSFC/SSAI
 *           Code 902
 *           Greenbelt, MD  20771
 *
 * ^CREATED: Nov 27, 1996 (original)
 *           Jul 17, 2003 (modified)
 *
 * ************************************************************************ */


#ifndef UARS3ATDEF_H
#define UARS3ATDEF_H

#define MIN3ATSIZE 148
#define MIN3ALSIZE 176
#define MAXRECSIZE 65535
#define FILL -9999.0
#define NaN 0x7fffffff


/*
 * Meta -- structure for the UARS metadata (*META file) or CDHF catalog dump
 */
struct Meta {
    char    type[12+1];		/* instrument type */
    char    subtype[12+1];	/* parameter or subtype */
    char    level[3+1];		/* data process level */
    char    day[4+1];		/* UARS day number */
    char    start_time[20+1];	/* start time of data */
    char    stop_time[20+1];	/* stop time of data */
    char    calib_id[12+1];	/* calibration id */
    char    version[4+1];	/* CCB version number */
    char    cycle[2+1];		/* file cycle number */
    char    file_name[80+1];	/* name of data file */
    char    file_size[5+1];	/* size of data file (original VMS blocks) */
    char    rec_size[5+1];	/* size of records in data file */
    char    uars_pi[24+1];	/* UARS PI name (correlative) */
    char    l3_base_idx[2+1];	/* base index of data */
    char    l3_num_pts[3+1];	/* number of points in data */
    char    base_wavelen[13+1];	/* units of altitude */
    char    min_alt[13+1];	/* units of altitude */
    char    max_alt[13+1];	/* units of altitude */
    char    units_alt[2+1];	/* units of altitude */
    char    min_lat[13+1];	/* min (south) latitude */
    char    max_lat[13+1];	/* max (north) latitude */
    char    min_lon[13+1];	/* min (west) longitude */
    char    max_lon[13+1];	/* max (east) longitude */
    char    source[12+1];	/* source of correlative */
    char    corr_pi[24+1];	/* correlative PI name */
    char    station_id[12+1];	/* station id */
    char    inst_id[12+1];	/* instrument id */
    char    param[5][12+1];	/* parameter names */
    char    comments[80+1];	/* comments (if any) */
    char    qual_uars[3+1];	/* UARS assigned quality for data */
    char    qual_pi[3+1];	/* PI assigned quality for data */
    char    create_job[21+1];	/* create job id */
};


/*
 * Sfdu -- structure for the UARS SFDU record
 */
struct Sfdu {
    char key[20];	/* VMS record key (for 3AL/3LP files only) */
    char Tz[12];	/* SFDU type (Tz) field */
    char Lz[8];		/* SFDU length (Lz) field */
    char Ti[12];	/* SFDU type (Ti) field */
    char Li[8];		/* SFDU length (Li) field */
};


/*
 * HeadL3A -- structure for the UARS level 3AL/3AT header (or label) record
 */
struct HeadL3A {
    char key[20];	/* VMS record key (for 3AL/3LP files only) */
    char Sat_ID[4];	/* Satellite identifier (always set to UARS) */
    char Rec_type[2];	/* Record type (always set to '1') */
    char Inst_ID[12];	/* Instrument identifier (CLAES, HALOE, etc.) */
    char Subtype[12];	/* Data subtype or species (H2O, O3, TEMP, etc.) */
    char Format_ver[4];	/* Format version number (always set to '1') */
    char Rec_count[8];	/* Physical record count (always set to '1', 1st rec) */
    char N_cont[4];	/* Number of continuation label records (0 for UARS) */
    char N_recs[8];	/* Number of physical records in file */
    char Create[23];	/* File creation time in VAX/VMS ASCII format */
    char Y_1st[3];	/* Year for first data record */
    char D_1st[3];	/* Day of year for first data record */
    char Ms_1st[8];	/* Milliseconds of day for first data record */
    char Y_last[3];	/* Year for last data record */
    char D_last[3];	/* Day of year for last data record */
    char Ms_last[8];	/* Milliseconds of day for last data record */
    char Level[3];	/* Data process level */
    char Day[4];	/* UARS day number */
    char N_pts[4];	/* Number of data points or 32-bit words per record */
    char Bs_index[4];	/* Base index of data points (not used in 3LP/3TP/3S) */
    char Rec_len[5];	/* Record length in bytes */
    char Min_lat[3];	/* Minimum Latitude (3AL/3LP only) */
    char Max_lat[3];	/* Maximum Latitude (3AL/3LP only) */
    char Version[9];	/* CCB version number */
    char Cycle[5];	/* File cycle number (wrong, see META file instead) */
    char Virt_flag[1];	/* Virtual file flag (always <blank> for UARS) */
    char TV_file[4];	/* Total number of time/version entries in file (= 0)*/
    char TV_rec[4];	/* Number of time/version entries in record (= 0)*/
};


/*
 * HeadL3S -- structure for the UARS level 3AS/3BS header (or label) record
 */
struct HeadL3S {
    char key[20];	/* VMS record key (for 3AL/3LP files only) */
    char Sat_ID[4];	/* Satellite identifier (always set to UARS) */
    char Rec_type[2];	/* Record type (always set to '1') */
    char Inst_ID[12];	/* Instrument identifier (CLAES, HALOE, etc.) */
    char Subtype[12];	/* Data subtype or species (H2O, O3, TEMP, etc.) */
    char Format_ver[4];	/* Format version number (always set to '1') */
    char Rec_count[8];	/* Physical record count (always set to '1', 1st rec) */
    char N_cont[4];	/* Number of continuation label records (0 for UARS) */
    char N_recs[8];	/* Number of physical records in file */
    char Create[23];	/* File creation time in VAX/VMS ASCII format */
    char Y_1st[3];	/* Year for first data record */
    char D_1st[3];	/* Day of year for first data record */
    char Ms_1st[8];	/* Milliseconds of day for first data record */
    char Y_last[3];	/* Year for last data record */
    char D_last[3];	/* Day of year for last data record */
    char Ms_last[8];	/* Milliseconds of day for last data record */
    char Level[3];	/* Data process level */
    char Day[4];	/* UARS day number */
    char N_pts[4];	/* Number of data points or 32-bit words per record */
    char Bs_wavelen[6];	/* Base wavelength of data point values */
    char Rec_len[5];	/* Record length in bytes */
    char Version[9];	/* CCB version number */
    char Cycle[5];	/* File cycle number (wrong, see META file instead) */
    char Virt_flag[1];	/* Virtual file flag (always <blank> for UARS) */
    char Ntv_file[4];	/* Total number of time/version entries in file (= 0)*/
    char Ntv_rec[4];	/* Number of time/version entries in record (= 0)*/
};


/*
 * HeadCORR -- structure for the UARS Correlative header (or label) record
 */
struct HeadCORR {
    char SFDU_Tz[12];	/* SFDU type (Tz) field */
    char SFDU_Lz[8];	/* SFDU length (Lz) field */
    char SFDU_Ti[12];	/* SFDU type (Ti) field */
    char SFDU_Li[8];	/* SFDU length (Li) field */
    char SFDU_Vi[48];	/* SFDU variable (Vi) field */
    char Project[4];	/* Project name (always set to UARS) */
    char UARS_PI[20];	/* UARS Principal Investigators name */
    char UARS_CMI[20];	/* UARS Correlative Measurement Investigators name */
    char Class[8];	/* Correlative data class */
    char Inst_ID[12];	/* Instrument identifier (CLAES, HALOE, etc.) */
    char Station_ID[12];/* Observing station identifier */
    char File_type[12];	/* Correlative data file type */
    char Start_time[23];/* Start time of data in file */
    char Stop_time[23];	/* Start time of data in file */
    char Max_lat[7];	/* Maximum latitude of data in file */
    char Min_lat[7];	/* Minimum latitude of data in file */
    char Max_lon[7];	/* Maximum longitude of data in file */
    char Min_lon[7];	/* Minimum longitude of data in file */
    char Max_alt_km[8];	/* Maximum altitude in kilometers of data in file */
    char Min_alt_km[8];	/* Minimum altitude in kilometers of data in file */
    char Max_alt_mb[8];	/* Maximum altitude in millibars of data in file */
    char Min_alt_mb[8];	/* Minimum altitude in millibars of data in file */
    char Rec_size[6];	/* Size of each record in the file */
    char N_recs[6];	/* The total number of records in the file */
    char Quality1[3];	/* Data quality word #1 */
    char Quality2[3];	/* Data quality word #2 */
    char Comments[80];	/* Comments supplied by the PI on the assmilation data*/
    char Param1[12];	/* Correlative data parameter #1 */
    char Param2[12];	/* Correlative data parameter #2 */
    char Param3[12];	/* Correlative data parameter #3 */
    char Param4[12];	/* Correlative data parameter #4 */
    char Param5[12];	/* Correlative data parameter #5 */
};


/*
 * DataL3A -- structure for the UARS level 3AL/3AT data record
 */
struct DataL3A {
    char key[20];	/* VMS record key (for 3AL/3LP files only) */
    char Sat_ID[4];	/* Satellite Identifier (always set to UARS) */
    char Rec_type[2];	/* Record type (always set to '1') */
    char Inst_ID[12];	/* Instrument Identifier (HALOE, ISAMS, etc.) */
    char Rec_count[8];	/* Physical record count (>= 1) */
    char Spare[2];	/* Spare - unused */
    int Total;		/* Total number of data points in the record */
    int Actual;		/* Number of actual data points in the record */
    int St_index;	/* Starting index of first actual point */
    int Time[2];	/* Record time in UDTF format [YYDDD,MMMMMMMM] */
    float Lat;		/* Geodetic latitude (-90 to +90) */
    float Lon;		/* Geodetic longitude (0 to 360) */
    float LST;		/* Local solar time (0 to 24) */
    float SZA;		/* Solar zenith angle (0 to 180) */
    float *data;	/* pointer to array of data values */
    float *quality;	/* pointer to array of standard deviations */
};


/*
 * ParmL3P -- structure for the UARS level 3LP/3TP parameter record
 */
struct ParmL3P {
    char key[20];	/* VMS record key (for 3AL/3LP files only) */
    char Sat_ID[4];	/* Satellite Identifier (always set to UARS) */
    char Rec_type[2];	/* Record type (always set to '1') */
    char Inst_ID[12];	/* Instrument Identifier (HALOE, ISAMS, etc.) */
    char Rec_count[8];	/* Physical record count (>= 2) */
    char Spare[2];	/* Spare - unused */
    int Total;		/* Total number of 32-bit words in the record */
    int Actual;		/* Number of actual 32-bit words in the record */
    char Spare2[4];	/* Spare - unused */
    int Time[2];	/* Record time in UDTF format [YYDDD,MMMMMMMM] */
    float Lat;		/* Geodetic latitude (-90 to +90) */
    float Lon;		/* Geodetic longitude (0 to 360) */
    char Spare3[8];	/* Spare - unused */
    int Num_words;	/* Number of 32-bit words */
    char *parm;		/* pointer to array of parameter values */
};


/*
 * DataL3S -- structure for the UARS level 3AS/3BS data record
 */
struct DataL3S {
    char key[20];	/* VMS record key (for 3AL/3LP files only) */
    char Sat_ID[4];	/* Satellite Identifier (always set to UARS) */
    char Rec_type[2];	/* Record type (always set to '1') */
    char Inst_ID[12];	/* Instrument Identifier (HALOE, ISAMS, etc.) */
    char Rec_count[8];	/* Physical record count (>= 1) */
    char Spare[2];	/* Spare - unused */
    int Total;		/* Total number of data points in the record */
    int Actual;		/* Number of actual data points in the record */
    float St_wavelen;	/* Starting wavelength of first actual point */
    int Time[2];	/* Record time in UDTF format [YYDDD,MMMMMMMM] */
    char Spare2[16];	/* Spare - unused */
    float *data;	/* pointer to array of data values */
    float *quality;	/* pointer to array of standard deviations */
    int Nparams;	/* number of solar parameters */
    char *param;	/* pointer to array of solar parameter indices */
};


/*
 * DataNMC -- structure for the UARS Correlative NMC data record
 */
struct DataNMC {
    int data_type;	/* The type of data contained in the file */
    int surf_type1;	/* First Surface type code */
    float surf_value1;	/* First Surface values */
    int surf_type2;	/* Second Surface type code */
    float surf_value2;	/* Second Surface values */
    int year;		/* Year of the data */
    int month;		/* Month of year of the data */
    int day;		/* Day of month of the data */
    int hour;		/* Hour of day of the data */
    int hemisphere;	/* Hemisphere or area type of the data */
    char time_value[23];/* Processing time of the data in the file */
    char extra_byte[1];	/* Extra byte to fill the record (spare) */
    float data[65][65];	/* Array of data values */
};


/*
 * DataUKMO -- structure for the UARS Correlative UKMO data record
 */
struct DataUKMO {
    int year_valid;	/* Year for which grid data are valid */
    int month_valid;	/* Month of year for which grid data are valid */
    int day_valid;	/* Day of month for which grid data are valid */
    int hour_valid;	/* Hour of day for which grid data are valid */
    int minutes_valid;	/* Minutes past hour for which grid data are valid */
    int yearday_valid;	/* Day of year for which grid data are valid */
    int year_data;	/* Year for data */
    int month_data;	/* Month of year for data */
    int day_data;	/* Day of month for data */
    int hour_data;	/* Hour of day for data */
    int minutes_data;	/* Minutes past hour for data */
    int yearday_data;	/* Day of year for data */
    int time_indicator;	/* Time indicator for data */
    int forecast_p;	/* Forecast period for data */
    int field_rec_size;	/* Size of data array */
    int grid_code;	/* Code for type of grid used */
    int hemisphere;	/* Code for hemisphere used */
    int num_rows;	/* Number of rows in data array */
    int num_cols;	/* Number of rows in data array */
    int extra_dsize;	/* Extra data size */
    int packing_method;	/* Packing method used to store data */
    int hdr_release;	/* Header release number */
    int field_code;	/* Field code for data */
    int field_2_code;	/* Second Field code for data */
    int process_code;	/* Processing code for data */
    int level_type;	/* Level type code */
    int ref_level_type;	/* Reference level type */
    int experiment_num;	/* Experiment number */
    int max_chunk_size;	/* Maximum chunk size */
    int extra_dsize_fp;	/* Extra data size Fpformat */
    int meto2_proj_num;	/* The fields file projection number (obsolete) */
    int meto2_fld_type;	/* The fields file field type code (obsolete) */
    int meto2_lvl_code;	/* The fields file level code (obsolete) */
    int reserved1[4];	/* Reserved for future use */
    int spare1;		/* Spare #1 for users use */
    int spare2[7];	/* Spare #2 for users use */
    float reserved2[4];	/* Reserved for future use */
    float datum;	/* Constant subtracted from each value in field */
    float pack_accuracy;/* Indicator of packing accuracy for packed fields */
    float level;	/* Level at which the grid point field is valid */
    float ref_level;	/* Reference level value (not used) */
    float level_a;	/* 'A' value of level of data (not used) */
    float ref_level_a;	/* 'A' value of reference level of data (not used) */
    float lat_ps_npol;	/* Real latitude of 'pseudo' N pole of projection */
    float lon_ps_npol;	/* Real longitude of 'pseudo' N pole of projection */
    float grid_orient;	/* Grid orientation (not used) */
    float lat_0_row;	/* Latitude of 'zeroth' row */
    float lat_interval;	/* Latitude intervals between rows */
    float lon_0_row;	/* Longitude of 'zeroth' row */
    float lon_interval;	/* Longitude intervals between rows */
    float missing_val;	/* Value used to indicate points with missing data */
    float scale;	/* Scaling factor of data */
    float *data;	/* Gridded Data Array */
};

union {
    struct DataL3A *l3a;
    struct DataL3S *l3s;
    struct ParmL3P *l3p;
} data;

struct Meta *ReadMETA(const char *);

void PrintSFDU(struct Sfdu *);
void PrintHeadL3A(struct HeadL3A *);
void PrintDataL3A(char *, char *, int, struct DataL3A *,
		  float, float, float, float, float, float);
void PrintDataL3AHeading(char *, char *);
void PrintDataL3S(char *, char *, int, struct DataL3S *);
void PrintHeadCorr(struct HeadCORR *);
void PrintDataUKMO(int, struct DataUKMO *,
		   float, float, float, float, float, float, char);
void PrintDataNMC(int, struct DataNMC *,
		  float, float, float, float, float, float);
void SubsetNMCGrid(float, float, float, float, int *, int **, int *, int **);
void w3fb05(float, float, float, float, float *, float *);
void w3fb04(float, float, float, float, float *, float *);
void Usage(void);			/* function prints usage of program */

#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

main(int argc, char *argv[])
{
    char file_name[1024];		/* name of the input file */
    char data_file[1024];		/* name of the data file */
    char type[12+1];			/* type of instrument or correlative */
    char subtype[12+1];			/* parameter or subtype of data */
    char lvlsrc[3+1];			/* processing level or source */
    int rec_size;			/* record size */
    int nrecs;				/* number of records */
    int nan=NaN;			/* not a number value */
    float slat=-90.0, nlat=90.0;	/* south and north latitude */
    float wlon=0.0, elon=360.0;		/* east and west longitude */
    float wavmin=0.0, wavmax=500.0;	/* start and end wavelength */
    float zmin, zmax;			/* min and max altitude */
    char uparam='a';			/* UKMO parameter */
    char swap;				/* flag to swap bytes */
    struct Meta *md;			/* structure for metadata record */
    struct Sfdu *sr;			/* structure for L3A SFDU record */
    struct HeadL3A *hr;			/* structure for L3A header record */
    struct DataL3A *dr;			/* structure for L3A data record */
    struct HeadL3S *hs;			/* structure for L3(A|B)S header rec */
    struct DataL3S *ds;			/* structure for L3(A|B)S data record */
    struct HeadCORR *hc;		/* structure for Correlative header */
    struct DataNMC *dn;			/* structure for NMC data record */
    struct DataUKMO *du;		/* structure for UKMO data record */

    zmin = *(float *) &nan;
    zmax = *(float *) &nan;

    if (ParseArgs(argc, argv, file_name, &slat, &nlat, &wlon, &elon,
		  &zmin, &zmax, &uparam, &wavmin, &wavmax, &swap) == -1)
    {
	fprintf(stderr, "Usage:  readuars [-?] [-b] [-s <slat> <nlat> "
			"[<wlon> <elon>]] [-v <zmin> <zmax>] [-u <uparam>] "
			"[-w <wavmin> <wavmax>] <filename>\n");
	exit(-1);
    }

    if (slat > nlat) {
	fprintf(stderr, "Error: south latitude must be less than north latitude\n");
	exit(-1);
    }

    if (fmod(elon - wlon, 360.0) != 0.0)
    {
	wlon = fmod(wlon, 360.0);
	if (wlon < 0) wlon = wlon + 360.0;
	elon = fmod(elon, 360.0);
	if (elon < 0) elon = elon + 360.0;
    }

    if (wavmin < 0 || wavmax < 0) {
	fprintf(stderr, "Error: wavelength must be a positive number\n");
	exit(-1);
    }

    if (wavmin > wavmax) {
	fprintf(stderr, "Error: min wavelength greater than max wavlength\n");
	exit(-1);
    }

    if (zmin < 0 || zmax < 0) {
	fprintf(stderr, "Error: altitude must be a positive number\n");
	exit(-1);
    }

    /* See if argument supplied file is a META file */
    if ((md = ReadMETA(file_name)) != NULL) {
	strcpy(data_file, md->file_name);
	rec_size = atoi(md->rec_size);
    }
    else {
	strcpy(data_file, file_name);
	if ((rec_size = GetRecSize(data_file, type, subtype, lvlsrc)) < 0)
	   exit(-1);
    }

    /* Check min and max altitude and set to defaults if necessary */

    if ((strncmp(type, "HRDI", 4) == 0 && strstr(subtype, "_A") != NULL) ||
	 strncmp(type, "WINDII", 6) == 0 || strncmp(type, "PEM", 3) == 0)
    {
	if (*(int *) &zmin == nan) zmin = 0.0;
	if (*(int *) &zmax == nan) zmax = 400.0;

	if (zmin > zmax) {	/* When using km zmin < zmax ! */
	    fprintf(stderr,"Error: min altitude not lower than max altitude\n");
	    exit(-1);
	}
    }
    else
    {
	if (*(int *) &zmin == nan) zmin = 1000.0;
	if (*(int *) &zmax == nan) zmax = 1.0e-6;

	if (zmin < zmax) {	/* When using hPa or mbar zmin > zmax ! */
	    fprintf(stderr,"Error: min altitude not lower than max altitude\n");
	    exit(-1);
	}
    }

    if (strncmp(lvlsrc, "NMC", 3) == 0) {

        if ((nrecs = ReadNMC(data_file, rec_size, swap, &hc, &dn)) < 0) {
            exit(-1);
        }

        PrintHeadCorr(hc);
        PrintDataNMC(nrecs, dn, slat, nlat, wlon, elon, zmin, zmax);
    }
    else if (strncmp(lvlsrc, "UKMO", 4) == 0) {

	if (uparam != 'a' && uparam != 'h' && uparam != 't' &&
            uparam != 'm' && uparam != 'z' && uparam != 'o')
        {
	    fprintf(stderr,"Error: not a valid choice for UKMO parameter\n");
	    exit(-1);
        }

        if ((nrecs = ReadUKMO(data_file, rec_size, swap, &hc, &du)) < 0) {
            exit(-1);
        }

        if (uparam == 'o' && strncmp(hc->Param5, "OMEGA", 5) != 0) {
	    fprintf(stderr,"Error: the OMEGA parameter is not in this file\n");
	    exit(-1);
        }

        PrintHeadCorr(hc);
        PrintDataUKMO(nrecs, du, slat, nlat, wlon, elon, zmin, zmax, uparam);
    }
    else if (strncmp(lvlsrc,"3AS",3) == 0 || strncmp(lvlsrc,"3BS",3) == 0) {

        if ((nrecs = ReadL3S(data_file, rec_size, swap, &sr, &hs, &ds)) < 0) {
            exit(-1);
        }

        PrintSFDU(sr);
/*      PrintHeadL3S(hs); */
	PrintDataL3S(type, subtype, nrecs, ds);
    }
    else {

        if ((nrecs = ReadL3A(data_file, rec_size, swap, &sr, &hr, &dr)) < 0) {
            exit(-1);
        }

        PrintSFDU(sr);
        PrintHeadL3A(hr);
	PrintDataL3A(type, subtype, nrecs, dr, slat, nlat, wlon, elon,
		     zmin, zmax);
    }

    exit(0);
}


/*
 * ReadMETA
 */

struct Meta *ReadMETA(const char *metafile)
{
    char   line[255+1];			/* string for lines in *META file */
    char   *ptr;			/* pointer for manipulation line */
    int    nlines;			/* number of lines */
    struct Meta meta, *md;		/* structure which contains metadata */
    FILE   *fp;				/* *META file pointer */

    if ((fp = fopen(metafile, "r")) == NULL) return NULL;

    memset(&meta, 0, sizeof(struct Meta));

    nlines = 0;
    do
    {
	fgets(line, 100, fp); nlines++;
        ptr = strchr(line, ':'); ptr+=2;

        if (strstr(line, "! TYPE : ") != NULL) {
            strncpy(meta.type, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! SUBTYPE : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.subtype, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! LEVEL : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.level, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! DAY : ") != NULL) {
            strncpy(meta.day, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! START_TIME : ") != NULL) {
            strncpy(meta.start_time, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! STOP_TIME : ") != NULL) {
            strncpy(meta.stop_time, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! CALIBRATION_ID : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.calib_id, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! VERSION : ") != NULL) {
            strncpy(meta.version, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! CYCLE : ") != NULL) {
            strncpy(meta.cycle, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! FILE_NAME : ") != NULL) {
            strncpy(meta.file_name, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! FILE_SIZE : ") != NULL) {
            strncpy(meta.file_size, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! RECORD_SIZE : ") != NULL) {
            strncpy(meta.rec_size, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! UARS_PI : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.uars_pi, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! L3_BASE_INDEX : ") != NULL) {
            strncpy(meta.l3_base_idx, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! L3_NBR_POINTS : ") != NULL) {
            strncpy(meta.l3_num_pts, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! BASE_WAVELENGTH : ") != NULL) {
            strncpy(meta.base_wavelen, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! MIN_ALTITUDE : ") != NULL) {
            strncpy(meta.min_alt, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! MAX_ALTITUDE : ") != NULL) {
            strncpy(meta.max_alt, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! UNITS_ALTITUDE : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.units_alt, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! MIN_LATITUDE : ") != NULL) {
            strncpy(meta.min_lat, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! MAX_LATITUDE : ") != NULL) {
            strncpy(meta.max_lat, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! MIN_LONGITUDE : ") != NULL) {
            strncpy(meta.min_lon, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! MAX_LONGITUDE : ") != NULL) {
            strncpy(meta.max_lon, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! SOURCE : ") != NULL) {
            strncpy(meta.source, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! CORRELATIVE_PI : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.corr_pi, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! STATION_ID : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.station_id, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! INSTRUMENT_ID : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.inst_id, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! PARAMETER(1) : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.param[0], ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! PARAMETER(2) : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.param[1], ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! PARAMETER(3) : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.param[2], ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! PARAMETER(4) : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.param[3], ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! PARAMETER(5) : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.param[4], ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! COMMENTS : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.comments, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! DATA_QUALITY_UARS : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.qual_uars, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! DATA_QUALITY_PI : ") != NULL) {
            if (ptr[0] != ' ')
            strncpy(meta.qual_pi, ptr, strlen(ptr)-1);
        }
        else if (strstr(line, "! CREATE_JOB : ") != NULL) {
            strncpy(meta.create_job, ptr, strlen(ptr)-1);
        }
        else {
	    if (nlines == 1) return NULL;
            fprintf(stderr, "Error: unknown META field %s\n", line);
	    break;
        }

    } while (! feof(fp));

    fclose(fp);

    if (strlen(meta.type) == 0) return NULL;

    md = &meta;
    return md;
}


/*
 * ReadL3A
 */

int ReadL3A(const char *file_name, int rec_size, char *swap,
	struct Sfdu **out_sr, struct HeadL3A **out_hr, struct DataL3A **out_dr)
{
    char rec[MAXRECSIZE];		/* array to hold records */
    char temp[MAXRECSIZE];		/* temporary array for manipulation */
    char sfdu[60+1];			/* array to hold SFDU record */
    char header[176+1];			/* array to hold header/label record */
    char *ptr;				/* string pointer */
    int i, j;				/* increment counters */
    int key;				/* size of index key (L3AL only) */
    int nrecs;				/* total number of data records */
    int npts;				/* number of data array elements */
    int *i32ptr;			/* pointer for 32-bit integers */
    int year, month, day, msec;		/* year, month, day, msec */
    float *f32ptr;			/* pointer for 32-bit floats */
    float *data, *dptr;			/* pointer for data value array */
    float *qual, *qptr;			/* pointer for data quality array */
    struct Sfdu *sr;			/* structure for SFDU record */
    struct HeadL3A *hr;			/* structure for header record */
    struct DataL3A *dr, *dtmp;		/* structure for data record */
    FILE *fp;				/* *PROD file pointer */


    if ((fp = fopen(file_name, "rb")) == NULL) {
        fprintf(stderr, "Error: unable to open %s\n", file_name);
        return -1;
    }

    /* the 0th record contains the SFDU information */
    fread(rec, 1, rec_size, fp);
    if ((ptr = strstr(rec, "CCSD1Z")) == NULL) {
	fclose(fp);
	fprintf(stderr, "Error: %s is not a valid UARS data file\n", file_name);
	return -1;
    }
    key = ptr - rec;
    sr = (struct Sfdu *) sfdu;
    memcpy(&sfdu[20-key], rec, 40+key); sfdu[40+20] = 0;

    /* the 1st record is the header or label */
    fread(rec, 1, rec_size, fp);
    memset(header, 0, sizeof(header));
    hr = (struct HeadL3A *) header;
    memcpy(&header[20-key], rec, 125+key); header[125+20] = 0;
    if (strncmp(hr->Level,"3AL",3) == 0 || strncmp(hr->Level,"3LP",3) == 0) {
	strncpy(&header[125+20], &rec[key+125], 29);
    }
    else {
  	strncpy(&header[131+20], &rec[key+125], 23);
    }

    /* get the number of data records in the file */
    strncpy(temp, hr->N_recs, 8); temp[8] = 0;
    nrecs = atoi(temp) - 1;  /* minus 1 because 1st record is the header */

    /* get the data array size (number of levels) in the file */
    strncpy(temp, hr->N_pts, 4); temp[4] = 0;
    npts = atoi(temp);

    /* the rest of the records are the data records, allocate memory */
    if ((dr = calloc(nrecs, sizeof(struct DataL3A))) == NULL) {
	fprintf(stderr, "Error: can't allocate memory for data records.\n");
	return -1;
    }
    if ((data = calloc(npts*nrecs, sizeof(float))) == NULL) {
	fprintf(stderr, "Error: can't allocate memory for data arrays.\n");
	return -1;
    }
    if ((qual = calloc(npts*nrecs, sizeof(float))) == NULL) {
	fprintf(stderr, "Error: can't allocate memory for quality arrays.\n");
	return -1;
    }

    memset(temp, 0, MAXRECSIZE);
    for (i=0, dptr=data, qptr=qual; i<nrecs; i++, dptr+=npts, qptr+=npts)
    {
	fread(rec, 1, rec_size, fp);

/* ---------------------------------------------------------------- *
 * This section is for byte swapping on little-endian machines.     *
 * NOTE: For level 3A files swapping is needed from byte 28+key on. *
 *       9 + 2*npts = number of 4 byte words to swap.               *
 * ---------------------------------------------------------------- */
	if (swap) {
	    ptr = rec + 28+key;
	    for (j = 0; j < 9 + 2*npts; j++, ptr+=4) byteswap(ptr, 4);
	}

	dtmp = (struct DataL3A *) temp;
	memcpy(&temp[20-key], rec, 64+key);

	f32ptr = (float *) &rec[64+key];
	i32ptr = (int *) &rec[64+key];

	for (j = 0; j < npts; j++, f32ptr++, i32ptr++) {
	    if (*i32ptr == 0x7fffffff) data[npts*i + j] = FILL;
	    else data[npts*i + j] = *f32ptr;
	}

	for (j = 0; j < npts; j++, f32ptr++, i32ptr++) {
	    if (*i32ptr == 0x7fffffff) qual[npts*i + j] = FILL;
	    else qual[npts*i + j] = *f32ptr;
	}

	dr[i] = *dtmp;
	dr[i].data = dptr;
	dr[i].quality = qptr;
    }
	
    fclose(fp);

    *out_sr = sr;
    *out_hr = hr;
    *out_dr = dr;

    return nrecs;
}


/*
 * ReadL3S
 */

int ReadL3S(const char *file_name, int rec_size, char *swap,
	struct Sfdu **out_sr, struct HeadL3S **out_hr, struct DataL3S **out_dr)
{
    char rec[MAXRECSIZE];		/* array to hold records */
    char temp[MAXRECSIZE];		/* temporary array for manipulation */
    char sfdu[60+1];			/* array to hold SFDU record */
    char header[176+1];			/* array to hold header/label record */
    char *ptr;				/* string pointer */
    int i, j;				/* increment counters */
    int key;				/* size of index key (L3AL only) */
    int nrecs;				/* total number of data records */
    int npts;				/* number of data array elements */
    int *i32ptr;			/* pointer for 32-bit integers */
    int year, month, day, msec;		/* year, month, day, msec */
    float *f32ptr;			/* pointer for 32-bit floats */
    float *data, *dptr;			/* pointer for data value array */
    float *qual, *qptr;			/* pointer for data quality array */
    struct Sfdu *sr;			/* structure for SFDU record */
    struct HeadL3S *hr;			/* structure for header record */
    struct DataL3S *dr, *dtmp;		/* structure for data record */
    FILE *fp;				/* *PROD file pointer */


    if ((fp = fopen(file_name, "rb")) == NULL) {
        fprintf(stderr, "Error: unable to open %s\n", file_name);
        return -1;
    }

    /* the 0th record contains the SFDU information */
    fread(rec, 1, rec_size, fp);
    if ((ptr = strstr(rec, "CCSD1Z")) == NULL) {
	fclose(fp);
	fprintf(stderr, "Error: %s is not a valid UARS data file\n", file_name);
	return -1;
    }
    key = ptr - rec;
    sr = (struct Sfdu *) sfdu;
    memcpy(&sfdu[20-key], rec, 40+key); sfdu[40+20] = 0;

    /* the 1st record is the header or label */
    fread(rec, 1, rec_size, fp);
    memset(header, 0, sizeof(header));
    hr = (struct HeadL3S *) header;
    memcpy(&header[20-key], rec, 150+key); header[150+20] = 0;

    /* get the number of data records in the file */
    strncpy(temp, hr->N_recs, 8); temp[8] = 0;
    nrecs = atoi(temp) - 1;  /* minus 1 because 1st record is the header */

    /* get the data array size (number of levels) in the file */
    strncpy(temp, hr->N_pts, 4); temp[4] = 0;
    npts = atoi(temp);

    /* the rest of the records are the data records, allocate memory */
    if ((dr = calloc(nrecs, sizeof(struct DataL3S))) == NULL) {
	fprintf(stderr, "Error: can't allocate memory for data records.\n");
	return -1;
    }
    if ((data = calloc(npts*nrecs, sizeof(float))) == NULL) {
	fprintf(stderr, "Error: can't allocate memory for data arrays.\n");
	return -1;
    }
    if ((qual = calloc(npts*nrecs, sizeof(float))) == NULL) {
	fprintf(stderr, "Error: can't allocate memory for quality arrays.\n");
	return -1;
    }

    memset(temp, 0, MAXRECSIZE);
    for (i=0, dptr=data, qptr=qual; i<nrecs; i++, dptr+=npts, qptr+=npts)
    {
	fread(rec, 1, rec_size, fp);

/* ---------------------------------------------------------------- *
 * This section is for byte swapping on little-endian machines.     *
 * NOTE: For level 3S files swapping is needed from byte 28+key on. *
 *       9 + 2*npts = number of 4 byte words to swap.               *
 * ---------------------------------------------------------------- */
	if (swap) {
	    ptr = rec + 28+key;
	    for (j = 0; j < 9 + 2*npts; j++, ptr+=4) byteswap(ptr, 4);
	}

	dtmp = (struct DataL3S *) temp;
	memcpy(&temp[20-key], rec, 64+key);

	f32ptr = (float *) &rec[64+key];
	i32ptr = (int *) &rec[64+key];

	for (j = 0; j < npts; j++, f32ptr++, i32ptr++) {
	    if (*i32ptr == 0x7fffffff) data[npts*i + j] = FILL;
	    else data[npts*i + j] = *f32ptr;
	}

	for (j = 0; j < npts; j++, f32ptr++, i32ptr++) {
	    if (*i32ptr == 0x7fffffff) qual[npts*i + j] = FILL;
	    else qual[npts*i + j] = *f32ptr;
	}

	if (swap) {
	    ptr = rec + 64 + 4*2*npts + key;
	    byteswap(ptr, 4);
        }

        i32ptr = (int *) &rec[64 + 4*2*npts + key];
        dtmp->Nparams = *i32ptr;

        dtmp->param = calloc(2*dtmp->Nparams, 20);
	memcpy(dtmp->param, &rec[64 + 4*2*npts + 4 + key], 2*20*dtmp->Nparams);

	dr[i] = *dtmp;
	dr[i].data = dptr;
	dr[i].quality = qptr;
    }
	

    fclose(fp);

    *out_sr = sr;
    *out_hr = hr;
    *out_dr = dr;

    return nrecs;
}


/*
 * ReadNMC
 */

int ReadNMC(const char *file_name, int rec_size, char *swap,
	struct HeadCORR **out_hc, struct DataNMC **out_dn)
{
    char rec[MAXRECSIZE];		/* array to hold records */
    char temp[MAXRECSIZE];		/* temporary array for manipulation */
    char sfdu[60+1];			/* array to hold SFDU record */
    char header[440+1];			/* array to hold header/label record */
    char *ptr;				/* string pointer */
    int i, j;				/* increment counters */
    int key;				/* size of index key (L3AL only) */
    int nrecs;				/* total number of data records */
    int npts;				/* number of data array elements */
    int offset;	
    int *i32ptr;			/* pointer for 32-bit integers */
    int year, month, day, msec;		/* year, month, day, msec */
    float *f32ptr;			/* pointer for 32-bit floats */
    float *data, *dptr;			/* pointer for data value array */
    struct HeadCORR *hc;		/* structure for header record */
    struct DataNMC *dn, *dtmp;		/* structure for data record */
    FILE *fp;				/* *PROD file pointer */


    if ((fp = fopen(file_name, "rb")) == NULL) {
        fprintf(stderr, "Error: unable to open %s\n", file_name);
        return -1;
    }

    /* the 1st record contains the header information */
    fread(rec, 1, rec_size, fp);
    if ((ptr = strstr(rec, "CCSD1Z")) == NULL) {
	fclose(fp);
	fprintf(stderr, "Error: %s is not a valid UARS data file\n", file_name);
	return -1;
    }
    memset(header, 0, sizeof(header));
    hc = (struct HeadCORR *) header;
    strncpy(header, rec, 380); header[380] = 0;

    /* get the number of data records in the file */
    strncpy(temp, hc->N_recs, 6); temp[6] = 0;
    nrecs = atoi(temp)-1;	/* Doesn't include the header record */

    /* the rest of the records in the file are the data records */
    if ((dn = calloc(nrecs, sizeof(struct DataNMC))) == NULL) {
	fprintf(stderr, "Error: can't allocate memory for data records.\n");
	return -1;
    }

    for (i=0; i<nrecs; i++)
    {
	fread(rec, 1, rec_size, fp);

/* ----------------------------------------------------- *
 * First swap the ten data header values                 *
 * ----------------------------------------------------- */
	if (swap) {
	    for (j=0, ptr=rec; j < 10; j++, ptr+=4) byteswap(ptr, 4);
	}

/* ----------------------------------------------------- *
 * Skip the ASCII time stamp, and continue swapping      *
 * the actual data values                                *
 * ----------------------------------------------------- */
  	if (swap) {
	    for (j=0, ptr+=24; j < 65*65; j++, ptr+=4) byteswap(ptr, 4);
	}

	dtmp = (struct DataNMC *) temp;
	memcpy(temp, rec, rec_size);
	dn[i] = *dtmp;
    }

    fclose(fp);

    *out_hc = hc;
    *out_dn = dn;

    return nrecs;
}


/*
 * ReadUKMO
 */

int ReadUKMO(const char *file_name, int rec_size, char *swap,
	struct HeadCORR **out_hc, struct DataUKMO **out_du)
{
    char rec[MAXRECSIZE];		/* array to hold records */
    char temp[MAXRECSIZE];		/* temporary array for manipulation */
    char sfdu[60];			/* array to hold SFDU record */
    char header[440+1];			/* array to hold header/label record */
    char *ptr;				/* string pointer */
    int i, j;				/* increment counters */
    int key;				/* size of index key (L3AL only) */
    int nrecs;				/* total number of data records */
    int npts;				/* number of data array elements */
    int offset;	
    int *i32ptr;			/* pointer for 32-bit integers */
    int year, month, day, msec;		/* year, month, day, msec */
    float *f32ptr;			/* pointer for 32-bit floats */
    float *data, *dptr;			/* pointer for data value array */
    struct HeadCORR *hc;		/* structure for header record */
    struct DataUKMO *du, *dtmp;		/* structure for data record */
    FILE *fp;				/* *PROD file pointer */


    if ((fp = fopen(file_name, "rb")) == NULL) {
        fprintf(stderr, "Error: unable to open %s\n", file_name);
        return -1;
    }

    /* the 1st record contains the header information */
    fread(rec, 1, rec_size, fp);
    if ((ptr = strstr(rec, "CCSD1Z")) == NULL) {
	fclose(fp);
	fprintf(stderr, "Error: %s is not a valid UARS data file\n", file_name);
	return -1;
    }
    hc = (struct HeadCORR *) header;
    memset(header, 0, sizeof(header));
    memcpy(header, rec, 86); header[86] = 0; header[87] = 0;
    memcpy(&header[88], &rec[86], 438-86); header[440] = 0;

    /* get the number of data records in the file */
    strncpy(temp, hc->N_recs, 6); temp[6] = 0;
    nrecs = atoi(temp)-1;	/* Doesn't include the header record */

    /* the rest of the records are data records, allocate memory */
    if ((du = calloc(nrecs, sizeof(struct DataUKMO))) == NULL) {
	fprintf(stderr, "Error: can't allocate memory for data records.\n");
	return -1;
    }
    if ((data = calloc(96*73*nrecs, sizeof(float))) == NULL) {
	fprintf(stderr, "Error: can't allocate memory for data array.\n");
	return -1;
    }

    for (i=0, dptr=data, offset=0; i<nrecs; i++)
    {
	fread(rec, 1, rec_size, fp);

/* ----------------------------------------------------- *
 * Swap all of the data values                           *
 * ----------------------------------------------------- */
	if (swap) {
	    for (j=0, ptr=rec; j < rec_size/4; j++, ptr+=4) byteswap(ptr, 4);
	}

	dtmp = (struct DataUKMO *) temp;
	memcpy(temp, rec, rec_size);
	du[i] = *dtmp;

        memcpy(&data[offset], &temp[256], du[i].field_rec_size*sizeof(float));
        du[i].data = dptr;

        dptr += du[i].field_rec_size;
        offset = dptr - data;
    }

    fclose(fp);

    *out_hc = hc;
    *out_du = du;

    return nrecs;
}


void PrintSFDU(struct Sfdu *sr)
{
    printf("-------------------------------\n");
    printf("    SFDU Descriptor\n");
    printf("-------------------------------\n");

/* Uncomment if you want to see this value. */
    if (strlen(sr->key) != 0) printf("Record key : %.20s\n", sr->key);

    printf("T[z] field : %.12s\n", sr->Tz);
    printf("L[z] field : %.8s\n", sr->Lz);
    printf("T[i] field : %.12s\n", sr->Ti);
    printf("L[i] field : %.8s\n\n", sr->Li);

    return;
}


void PrintHeadL3A(struct HeadL3A *hr)
{
    printf("-------------------------------\n");
    printf("    Header Record\n");
    printf("-------------------------------\n");

/* Uncomment if you want to see this value. */
    if (strlen(hr->key) != 0) printf("Record key : %.20s\n", hr->key);

    printf("Satellite ID  : %.4s\n", hr->Sat_ID);
/* Uncomment if you want to see this value. */
    printf("Record type   : %.2s\n", hr->Rec_type);
    
    printf("Instrument ID : %.12s\n", hr->Inst_ID);
    printf("Data subtype or species : %.12s\n", hr->Subtype);
/* Uncomment if you want to see these values. */
    printf("Format version number   : %.4s\n", hr->Format_ver);
    printf("Physical record count   : %.8s\n", hr->Rec_count);
    printf("Number of continuation records for file label : %.4s\n",hr->N_cont);
    printf("Number of physical records in file : %.8s\n", hr->N_recs);
    
    printf("File Creation Time : %.23s\n", hr->Create);
    printf("Year of 1st data record  : %.3s\n", hr->Y_1st);
    printf("Day of 1st data record   : %.3s\n", hr->D_1st);
    printf("Milliseconds of 1st data record  : %.8s\n", hr->Ms_1st);
    printf("Year of last data record : %.3s\n", hr->Y_last);
    printf("Day of last data record  : %.3s\n", hr->D_last);
    printf("Milliseconds of last data record : %.8s\n", hr->Ms_last);
    printf("Data level : %.3s\n", hr->Level);
    printf("UARS day number : %.4s\n", hr->Day);
/* Uncomment if you want to see these values. */
    printf("Number of data points per record : %.4s\n", hr->N_pts);
    printf("Base index of data point value : %.4s\n", hr->Bs_index);
    printf("Record length in bytes : %.5s\n", hr->Rec_len);
    
    if (strncmp(hr->Level,"3AL",3) == 0 || strncmp(hr->Level,"3LP",3) == 0) {
        printf("Minimum latitude for records in file : %.3s\n", hr->Min_lat);
        printf("Maximum latitude for records in file : %.3s\n", hr->Max_lat);
    }
    printf("CCB version number : %.9s\n", hr->Version);
    printf("File cycle number      : %.5s\n", hr->Cycle);
/* Uncomment if you want to see these values. */
    printf("Virtual file flag : %.1s\n", hr->Virt_flag);
    printf("Total number of time/version entries in file : %.4s\n",hr->TV_file);
    printf("Number of time/version entries in record : %.4s\n", hr->TV_rec);

    return;
}


void PrintDataL3A(char *instrmt, char *subtype, int nrecs, struct DataL3A *dr,
	 float slat, float nlat, float wlon, float elon, float zmin, float zmax)
{
    int i, j;			/* increment counters */
    int izmin, izmax;		/* UARS altitude index for zmin and zmax */
    float pr_alt[90];		/* pressure/altitude array */
    struct DataL3A *d;		/* structure for data record */

    /* First calculate the pressure or altitude level values */
    if ((strncmp(instrmt,"HRDI",4) == 0 && strstr(subtype,"_A") != NULL) ||
	 strncmp(instrmt,"WINDII",6) == 0 ||
         strncmp(instrmt,"PEM",3) == 0)
    {
	for (i=0; i<90; i++)
	{
	    if (i < 13) pr_alt[i] = 5.0*(float)i;
	    else if (i < 33) pr_alt[i] = 60.0 + 3.0*(float)(i - 12);
	    else pr_alt[i] = 120.0 + 5.0*(float)(i - 32);
	}

	if (zmin < 60) izmin = (int)(zmin/5.0 + 0.5);
	else if (zmin < 120) izmin = (int)((zmin-60.0)/3.0 + 12 + 0.5);
	else izmin = (int)((zmin-120.0)/5.0 + 32 + 0.5);

	if (zmax < 60) izmax = (int)(zmax/5.0 + 0.5);
	else if (zmax < 120) izmax = (int)((zmax-60.0)/3.0 + 12 + 0.5);
	else izmax = (int)((zmax-120.0)/5.0 + 32 + 0.5);
    }
    else {
	for (i=0; i<45; i++) pr_alt[i] = 1000.0*pow(10.0, -((float)i/6.0));

        izmin = (int)(6*(3 - log10(zmin)) + 0.5);
        izmax = (int)(6*(3 - log10(zmax)) + 0.5);
    }


    for (i=0, d=dr; i<nrecs; i++, d++)
    {
	if (d->Lat < slat || nlat < d->Lat) continue;

	if (wlon > elon) {	/* Crosses Greenwich */
	    if (d->Lon < wlon && elon < d->Lon) continue;
        }
	else {
	    if (d->Lon < wlon || elon < d->Lon) continue;
        }

	printf("\n---------------------------------------------\n");
	printf("    Data Record Number %.8s\n", d->Rec_count);
	printf("---------------------------------------------\n");

/* These are always the same for each record so they are commented out. */
	if (strlen(d->key) != 0) printf("Record key : %.20s\n", d->key);
	printf("Satellite ID  : %.4s\n", d->Sat_ID);
	printf("Record type   : %.2s\n", d->Rec_type);
	printf("Instrument ID : %.12s\n", d->Inst_ID);

/* If you want to see these values just uncomment this section. */
	printf("Physical record count : %.8s\n", d->Rec_count);
	printf("Total number points in record : %d\n", d->Total);
	printf("Number of Actual points : %d\n", d->Actual);
	printf("Starting index of first actual point : %d\n", d->St_index);

	printf("Time in UDTF format : %6d %8d (yyddd, millisec)\n",
		d->Time[0], d->Time[1]);
	printf("Latitude  : %f (degrees)\n", d->Lat);
	printf("Longitude : %f (degrees)\n", d->Lon);
	printf("Local solar time   : %f (hours)\n", d->LST);
	printf("Solar zenith angle : %f (degrees)\n\n", d->SZA);

        PrintDataL3AHeading(instrmt, subtype);

	for (j = d->St_index; j < d->Actual + d->St_index; j++)
	{
	    if (j < izmin || j > izmax) continue;

/*	    int4 = (int *) &line[64+(j-d->St_index)*4]; */ /* Check for NaN fills */
/*	    if (*int4 != 0x7fffffff) { */
		printf("   %2d  %13g  %#13g       %#13g\n", j, pr_alt[j],
			d->data[j-1], d->quality[j-1]);
/*			d->data[j - d->St_index], d->quality[j - d->St_index]); */
/*	    }
	    else {
		printf("   %2d  %13g        NaN                 NaN\n",
			j, pr_alt[j]);
	    } */
	}
    }

    return;
}


void PrintDataL3S(char *instrmt, char *subtype, int nrecs, struct DataL3S *dr)
	 
{
    int i, j;			/* increment counters */
    int izmin, izmax;		/* UARS altitude index for zmin and zmax */
    float pr_alt[90];		/* pressure/altitude array */
    struct DataL3S *d;		/* structure for data record */


    for (i=0, d=dr; i<nrecs; i++, d++)
    {
	printf("\n---------------------------------------------\n");
	printf("    Data Record Number %.8s\n", d->Rec_count);
	printf("---------------------------------------------\n");

/* These are always the same for each record so they are commented out. */
	printf("Satellite ID  : %.4s\n", d->Sat_ID);
	printf("Record type   : %.2s\n", d->Rec_type);
	printf("Instrument ID : %.12s\n", d->Inst_ID);

/* If you want to see these values just uncomment this section. */
	printf("Physical record count : %.8s\n", d->Rec_count);
	printf("Total number points in record : %d\n", d->Total);
	printf("Number of Actual points : %d\n", d->Actual);
	printf("Starting wavelength of first actual point : %f\n", d->St_wavelen);

	printf("Time in UDTF format : %6d %8d (yyddd, millisec)\n\n",
		d->Time[0], d->Time[1]);

	printf("Index  Wavelength(nm) Data (W/m^2)        Quality (W/m^2)\n");
        printf("-----  -------------  ------------------  ---------------------\n");

	for (j = 0; j < d->Actual; j++)
	{
/*	    if (j < wavmin || j > wavmax) continue; */

/*	    int4 = (int *) &line[64+(j-d->St_index)*4]; */ /* Check for NaN fills */
/*	    if (*int4 != 0x7fffffff) { */
		printf("  %3d  %13g  %#13g       %#13g\n", j, d->St_wavelen+j,
			d->data[j], d->quality[j]);
/*	    }
	    else {
		printf("  %3d  %13g        NaN                 NaN\n",
			j, d->St_wavelen+j);
	    } */
	}

	printf("\nAdditional Parameters : %d\n", d->Nparams);

        for (j = 0; j < d->Nparams; j++)
        {
           printf("  %2d) %.20s = %.20s\n", j,
                   &d->param[j*40], &d->param[j*40 + 20]);
        }
    }

    return;
}


void PrintDataL3AHeading(char *instrmt, char *subtype)
{
    if (strncmp(instrmt, "PEM", 3) == 0)
    {
	if (strncmp(subtype, "EDEP3AT_ELEC", 12) == 0 ||
	    strncmp(subtype, "EDEP3AT_PROT", 12) == 0)
	    printf("Index  Altitude (km)  "
		   "Data (erg/cm**3/s)  Quality (erg/cm**3/s)\n");
	else
	    printf("Index  Altitude (km)  "
		   "Data (keV/g/s)      Quality (keV/g/s)\n");
    }
    else if (strncmp(instrmt, "HRDI", 4) == 0 ||
	     strncmp(instrmt, "WINDII", 6) == 0)
    {

	if      (strncmp(subtype, "MERWIN_P", 8) == 0 ||
	         strncmp(subtype, "ZONWIN_P", 8) == 0)
	    printf("Index  Pressure (mb)  "
		   "Data (m/s)          Quality (m/s)**2\n");
	else if (strncmp(subtype, "MERWIN_A", 8) == 0 ||
                 strncmp(subtype, "ZONWIN_A", 8) == 0)
	    printf("Index  Altitude (km)  "
		   "Data (m/s)          Quality (m/s)**2\n");
	else if (strncmp(subtype, "VOLER_P", 7) == 0)
	    printf("Index  Pressure (mb)  "
		   "Data (1/cm**3/s)    Quality (1/cm**3/s)**2\n");
	else if (strncmp(subtype, "VOLER_A", 7) == 0)
	    printf("Index  Altitude (km)  "
		   "Data (1/cm**3/s)    Quality (1/cm**3/s)**2\n");
	else if (strncmp(subtype, "TEMP_P", 6) == 0)
	    printf("Index  Pressure (mb)  "
		   "Data (Kelvin)       Quality (Kelvin)**2\n");
	else
	    printf("Index  Altitude (km)  "
		   "Data (Kelvin)       Quality (Kelvin)**2\n");
    }
    else
    {
	if      (strncmp(subtype, "AEXT", 4) == 0 ||
	         strncmp(subtype, "AERO", 4) == 0)
	    printf("Index  Pressure (mb)  "
		   "Data (1/km)         Quality (1/km)\n");
	else if (strncmp(subtype, "GPH", 3) == 0 ||
		 strncmp(subtype, "ALT", 3) == 0)
	    printf("Index  Pressure (mb)  "
		   "Data (m)            Quality (m)\n");
	else if (strncmp(subtype, "UTH", 3) == 0)
	    printf("Index  Pressure (mb)  "
		   "Data (%%)            Quality (%%)\n");
	else if (strncmp(subtype, "TEMP", 4) == 0)
	    printf("Index  Pressure (mb)  "
		   "Data (Kelvin)       Quality (Kelvin)\n");
	else
	    printf("Index  Pressure (mb)  "
		   "Data (vmr)          Quality (vmr)\n");
    }

    printf("-----  -------------  ------------------  ---------------------\n");

    return;
}


void PrintHeadCorr(struct HeadCORR *hc)
{
    printf("---------------------------------------------\n");
    printf("    Header Record\n");
    printf("---------------------------------------------\n");

    printf("SFDU T[z] field               : %.12s\n", hc->SFDU_Tz);
    printf("SFDU L[z] field               : %.8s\n", hc->SFDU_Lz);
    printf("SFDU T[i] field               : %.12s\n", hc->SFDU_Ti);
    printf("SFDU L[i] field               : %.8s\n", hc->SFDU_Li);
    printf("SFDU V[i] field               : %.48s\n", hc->SFDU_Vi);
    printf("Project Name                  : %.4s\n", hc->Project);
    printf("UARS PI                       : %.20s\n", hc->UARS_PI);
    printf("UARS CMI                      : %.20s\n", hc->UARS_CMI);
    printf("Correlative Data Class        : %.8s\n", hc->Class);

/* Not used */
/*  printf("Instrument Type               : %.12s\n", hc->Inst_ID); */
/*  printf("Observing Station ID          : %.12s\n", hc->Station_ID); */

    printf("Correlative Data File Type ID : %.12s\n", hc->File_type);
    printf("Start Time of File            : %.23s\n", hc->Start_time);
    printf("Stop Time of File             : %.23s\n", hc->Stop_time);
    printf("Maximum Latitude              : %.7s\n", hc->Max_lat);
    printf("Minimum Latitude              : %.7s\n", hc->Min_lat);
    printf("Maximum Longitude             : %.7s\n", hc->Max_lon);
    printf("Minimum Longitude             : %.7s\n", hc->Min_lon);

/* Not used */
/*  printf("Maximum Altitude Kilometers   : %.8s\n", hc->Max_alt_km); */
/*  printf("Minimum Altitude Kilometers   : %.8s\n", hc->Min_alt_km); */

    printf("Maximum Altitude Millibars    : %.8s\n", hc->Max_alt_mb);
    printf("Minimum Altitude Millibars    : %.8s\n", hc->Min_alt_mb);
    printf("Record Size                   : %.6s\n", hc->Rec_size);
    printf("Number Records in File        : %.6s\n", hc->N_recs);

/* Not used */
/*  printf("Data Quality Word #1          : %.3s\n", hc->Quality1); */
/*  printf("Data Quality Word #2          : %.3s\n", hc->Quality2); */

    printf("User Comments                 : %.80s\n", hc->Comments);
    if (strncmp(hc->Param1, "            ", 12) != 0 && strlen(hc->Param1) > 0)
	printf("Correlative Data Parameter #1 : %.12s\n", hc->Param1);
    if (strncmp(hc->Param2, "            ", 12) != 0 && strlen(hc->Param2) > 0)
	printf("Correlative Data Parameter #2 : %.12s\n", hc->Param2);
    if (strncmp(hc->Param3, "            ", 12) != 0 && strlen(hc->Param3) > 0)
	printf("Correlative Data Parameter #3 : %.12s\n", hc->Param3);
    if (strncmp(hc->Param4, "            ", 12) != 0 && strlen(hc->Param4) > 0)
	printf("Correlative Data Parameter #4 : %.12s\n", hc->Param4);
    if (strncmp(hc->Param5, "            ", 12) != 0 && strlen(hc->Param5) > 0)
	printf("Correlative Data Parameter #5 : %.12s\n", hc->Param5);

    return;
}


void PrintDataNMC(int nrecs, struct DataNMC *drec,
	float slat, float nlat, float wlon, float elon, float zmin, float zmax)
{
    int i, j, k, n;		/* increment counters */
    int lat1st, latend;		/* first latitude and end latitude point */
    int lon1st, lonend;		/* first longitude and end longitude point */
    int izmin, izmax, izval;	/* altitude index min and max and value */
int num_n, *idx_n;
int num_s, *idx_s;
    char prname[15];		/* parameter name and units */
    struct DataNMC *d;		/* structure for data record */
    float latgrd[65][65][2];	/* north and south lat cell center points */
    float latedg[66][66][2];	/* north and south lat cell corner points */
    float longrd[65][65][2];	/* north and south lon cell center points */
    float lonedg[66][66][2];	/* north and south lon cell corner points */
				/* NMC pressure levels */
    float nmclvl[18] = {1000, 850, 700, 500, 400, 300, 250, 200,
			150, 100, 70, 50, 30, 10, 5, 2, 1, 0.4};

    izmin = 18;
    izmax = 1;
    for (i=0; i<18-1; i++) {
        j = 18 - i;
        if (zmin >= nmclvl[j]-(nmclvl[j]-nmclvl[j+1])/2.0) izmin = j;
        if (zmax <= nmclvl[i]-(nmclvl[i]-nmclvl[i+1])/2.0) izmax = i+1;
    }

    /* First set up latitude and longitude grid cell center points. */
    for (j=0; j<65; j++)
    {
	for (i=0; i<65; i++)
	{
	    /* Because the data were originally made on a VAX, reverse the row
             * and columns in the output grid arrays for latitude and longitude
             * to match the array order of the data array.
             */
	    w3fb05(i-32.0, j-32.0, 381, 80, &latgrd[i][j][0], &longrd[i][j][0]);
	    w3fb05(i-32.0, j-32.0,-381,260, &latgrd[i][j][1], &longrd[i][j][1]);

	    /* w3fb05 returns west longitude, so change to east longitude! */
	    longrd[i][j][0] = 360 - longrd[i][j][0];
	    longrd[i][j][1] = 360 - longrd[i][j][1];
	}
    }

    /* Next set up latitude and longitude grid cell corner points. */
    for (j=0; j<66; j++)
    {
	for (i=0; i<66; i++)
	{
	    w3fb05(i-32.5, j-32.5, 381, 80, &latedg[i][j][0], &lonedg[i][j][0]);
	    w3fb05(i-32.5, j-32.5,-381,260, &latedg[i][j][1], &lonedg[i][j][1]);

	    /* w3fb05 returns west longitude, so change to east longitude! */
	    lonedg[i][j][0] = 360 - lonedg[i][j][0];
	    lonedg[i][j][1] = 360 - lonedg[i][j][1];
	}
    }

    SubsetNMCGrid(slat, nlat, wlon, elon, &num_n, &idx_n, &num_s, &idx_s);

    for (n=0, d=drec; n<nrecs; n++, d++)
    {
        if (d->hemisphere == 27 && num_n <= 0) continue;
        else if (d->hemisphere == 28 && num_s <= 0) continue;

	if (d->surf_value1 > nmclvl[izmin] || d->surf_value1 < nmclvl[izmax]) {
	    continue;
        }

	printf("\n---------------------------------------------\n");
	printf("    Data Record Number %2d\n", n+1);
	printf("---------------------------------------------\n");

  	printf("Data Type         : %d\n", d->data_type);
	printf("Surface Type 1    : %d\n", d->surf_type1);
	printf("Surface Value 1   : %.1f\n", d->surf_value1);
/* Not used */
/*	printf("Surface Type 2    : %d\n", d->surf_type2);
	printf("Surface Value 2   : %.1f\n", d->surf_value2);
*/
	printf("Year Data         : %d\n", d->year);
	printf("Month Data        : %02d\n", d->month);
	printf("Day Data          : %02d\n", d->day);
	printf("Hour Data         : %02d\n", d->hour);
	printf("Hemisphere        : %d\n", d->hemisphere);
	printf("Time Indicator    : %.23s\n\n", d->time_value);

	/* Write out values for the data grid */
	printf("       NMC 65 by 65 Data Grid\n");

        if (d->hemisphere == 27)
        {
	    printf("        Northern Hemisphere:\n\n");
	    printf("Col Row Longitude  Latitude  Data Value\n");
	    printf("--- --- --------- --------- -----------\n");

	    for (k=0; k<num_n; k++)
	    {
		i = idx_n[2*k+0];
		j = idx_n[2*k+1];

		if((slat < latedg[i][j][0] && nlat > latedg[i][j][0] &&
		    wlon < lonedg[i][j][0] && elon > lonedg[i][j][0]) ||
		   (slat < latedg[i+1][j][0] && nlat > latedg[i+1][j][0] &&
		    wlon < lonedg[i+1][j][0] && elon > lonedg[i+1][j][0]) ||
		   (slat < latedg[i+1][j+1][0] && nlat > latedg[i+1][j+1][0] &&
		    wlon < lonedg[i+1][j+1][0] && elon > lonedg[i+1][j+1][0]) ||
		   (slat < latedg[i][j+1][0] && nlat > latedg[i][j+1][0] &&
		    wlon < lonedg[i][j+1][0] && elon > lonedg[i][j+1][0]))
		{
		    printf("%3d %3d %#9.5f %#9.5f %#11.5g\n", i-32, j-32,
			   longrd[i][j][0], latgrd[i][j][0], d->data[j][i]);
		}
	    }
	}
	else if (d->hemisphere == 28)
	{
	    printf("        Southern Hemisphere:\n\n");
	    printf("Col Row Longitude  Latitude  Data Value\n");
	    printf("--- --- --------- --------- -----------\n");

	    for (k=0; k<num_s; k++)
	    {
		i = idx_s[2*k+0];
		j = idx_s[2*k+1];

		if((slat < latedg[i][j][1] && nlat > latedg[i][j][1] &&
		    wlon < lonedg[i][j][1] && elon > lonedg[i][j][1]) ||
		   (slat < latedg[i+1][j][1] && nlat > latedg[i+1][j][1] &&
		    wlon < lonedg[i+1][j][1] && elon > lonedg[i+1][j][1]) ||
		   (slat < latedg[i+1][j+1][1] && nlat > latedg[i+1][j+1][1] &&
		    wlon < lonedg[i+1][j+1][1] && elon > lonedg[i+1][j+1][1]) ||
		   (slat < latedg[i][j+1][1] && nlat > latedg[i][j+1][1] &&
		    wlon < lonedg[i][j+1][1] && elon > lonedg[i][j+1][1]))
		{
		    printf("%3d %3d %#9.5f %#9.5f %#11.5g\n", i-32, j-32,
			   longrd[i][j][1], latgrd[i][j][1], d->data[j][i]);
		}
	    }
	}
	else
	{
	    printf("        Unknown Hemisphere!\n\n");
	}
    }

    return;
}

void SubsetNMCGrid( float slat, float nlat, float wlon, float elon,
		    int *out_n, int **out_idxn, int *out_s, int **out_idxs )
{
    int i, j;
    int num_n, num_s;
    int *idx_n, *idx_s, *n, *s;
    float lonedg[66][66][2];	/* east and west lon cell corner points */
    float latedg[66][66][2];	/* north and south lat cell corner points */


    /* Set up latitude and longitude grid cell corner points to
     * find out which grid cell is in the spatial search box.
     */
    for (j=0; j<66; j++)
    {
	for (i=0; i<66; i++)
	{
	    w3fb05(i-32.5, j-32.5, 381, 80, &latedg[i][j][0], &lonedg[i][j][0]);
	    w3fb05(i-32.5, j-32.5,-381,260, &latedg[i][j][1], &lonedg[i][j][1]);

	    /* w3fb05 returns west longitude, so change to east longitude! */
	    lonedg[i][j][0] = 360 - lonedg[i][j][0];
	    lonedg[i][j][1] = 360 - lonedg[i][j][1];
	}
    }

    idx_n = calloc(65*65*2, sizeof(int));
    idx_s = calloc(65*65*2, sizeof(int));

    for (j=0, n=idx_n, s=idx_s, num_n=0, num_s=0; j<65; j++)
    {
	for (i=0; i<65; i++)
	{
	    if ( ( slat < latedg[i][j][0] && nlat > latedg[i][j][0] &&
		   wlon < lonedg[i][j][0] && elon > lonedg[i][j][0] ) ||
		 ( slat < latedg[i+1][j][0] && nlat > latedg[i+1][j][0] &&
		   wlon < lonedg[i+1][j][0] && elon > lonedg[i+1][j][0] ) ||
		 ( slat < latedg[i+1][j+1][0] && nlat > latedg[i+1][j+1][0] &&
		   wlon < lonedg[i+1][j+1][0] && elon > lonedg[i+1][j+1][0] ) ||
		 ( slat < latedg[i][j+1][0] && nlat > latedg[i][j+1][0] &&
		   wlon < lonedg[i][j+1][0] && elon > lonedg[i][j+1][0] ) )
	    {
		num_n++;
		*n=i; n++;
		*n=j; n++;
	    }
	    if ( ( slat < latedg[i][j][1] && nlat > latedg[i][j][1] &&
		   wlon < lonedg[i][j][1] && elon > lonedg[i][j][1] ) ||
		 ( slat < latedg[i+1][j][1] && nlat > latedg[i+1][j][1] &&
		   wlon < lonedg[i+1][j][1] && elon > lonedg[i+1][j][1] ) ||
		 ( slat < latedg[i+1][j+1][1] && nlat > latedg[i+1][j+1][1] &&
		   wlon < lonedg[i+1][j+1][1] && elon > lonedg[i+1][j+1][1] ) ||
		 ( slat < latedg[i][j+1][1] && nlat > latedg[i][j+1][1] &&
		   wlon < lonedg[i][j+1][1] && elon > lonedg[i][j+1][1] ) )
	    {
		num_s++;
		*s=i; s++;
		*s=j; s++;
	    }
	}
    }

    /* realloc? let's not waste unneeded memory! */
    if (num_n > 0) {
	idx_n = realloc(idx_n, 2*num_n*sizeof(int));
    }
    else {
	free(idx_n);
	idx_n = NULL;
    }

    if (num_s > 0) {
	idx_s = realloc(idx_s, 2*num_s*sizeof(int));
    }
    else {
	free(idx_s);
	idx_s = NULL;
    }

    *out_idxn = idx_n;
    *out_n = num_n;
    *out_idxs = idx_s;
    *out_s = num_s;


    return;
}


void PrintDataUKMO(int nrecs, struct DataUKMO *drec,
	float slat, float nlat, float wlon, float elon,
	float zmin, float zmax, char uparam)
{
    int i, j, k, n;		/* increment counters */
    int lat1st, latend;		/* first latitude and end latitude point */
    int lon1st, lonend;		/* first longitude and end longitude point */
    int izmin, izmax, izval;	/* altitude index min and max and value */
    char prname[15];		/* parameter name and units */
    struct DataUKMO *d;		/* structure for data record */

    izmin = (int)(6*(3 - log10(zmin)) + 0.5);
    izmax = (int)(6*(3 - log10(zmax)) + 0.5);

    for (n=0, d=drec; n<nrecs; n++, d++)
    {
	if ( d->field_code !=  1 && d->field_code != 16 &&
             d->field_code != 48 && d->field_code != 49 &&
             d->field_code != 56 && d->field_code != 57 &&
             d->field_code != 40 ) continue;

        if ( d->field_code ==  1 && (uparam != 'a' && uparam != 'h') ) continue;
        if ( d->field_code == 16 && (uparam != 'a' && uparam != 't') ) continue;
        if ( d->field_code == 48 || d->field_code == 56 &&
	    (uparam != 'a' && uparam != 'z') ) continue;
        if ( d->field_code == 49 || d->field_code == 57 &&
	    (uparam != 'a' && uparam != 'm') ) continue;
        if ( d->field_code == 40 &&
	    (uparam != 'a' && uparam != 'o') ) continue;

	izval = (int)(6*(3-log10(d->level)) + 0.5);
	if (izval < izmin || izval > izmax) continue;

	printf("\n---------------------------------------------\n");
	printf("    Data Record Number %2d\n", n+1);
	printf("---------------------------------------------\n");

/* The valid times and data times are the same */
/*
	printf("Year Valid               : %4d\n", d->year_valid);
	printf("Month Valid              : %02d\n", d->month_valid);
	printf("Day Valid                : %02d\n", d->day_valid);
	printf("Hour Valid               : %02d\n", d->hour_valid);
	printf("Minutes Valid            : %02d\n", d->minutes_valid);
	printf("Yearday Valid            : %3d\n", d->yearday_valid);
*/
	printf("Year Data                : %4d\n", d->year_data);
	printf("Month Data               : %02d\n", d->month_data);
	printf("Day Data                 : %02d\n", d->day_data);
	printf("Hour Data                : %02d\n", d->hour_data);
	printf("Minutes Data             : %02d\n", d->minutes_data);
	printf("Yearday Data             : %3d\n", d->yearday_data);
	printf("Time Indicator           : %d\n", d->time_indicator);
	printf("Forecast Period          : %d\n", d->forecast_p);
	printf("Field Record Size        : %d\n", d->field_rec_size);
	printf("Grid Code                : %d\n", d->grid_code);
	printf("Hemisphere               : %d\n", d->hemisphere);
	printf("Number of Rows           : %d\n", d->num_rows);
	printf("Number of Columns        : %d\n", d->num_cols);
	printf("Extra Data Size          : %d\n", d->extra_dsize);
	printf("Packing Method           : %d\n", d->packing_method);
	printf("Header Release Number    : %d\n", d->hdr_release);
	printf("Field Code               : %d\n", d->field_code);
	printf("Second Field Code        : %d\n", d->field_2_code);
	printf("Processing Code          : %d\n", d->process_code);
	printf("Level Type               : %d\n", d->level_type);
	printf("Reference Level Type     : %d\n", d->ref_level_type);
	printf("Experiment Number        : %d\n", d->experiment_num);
	printf("Maximum Chunk Size       : %d\n", d->max_chunk_size);
	printf("Extra Data Size Fpformat : %d\n", d->extra_dsize_fp);

/* The following are not used */
/*
	printf("Meto2 Project Number     : %d\n", d->meto2_proj_num);
	printf("Meto2 Field Type Code    : %d\n", d->meto2_fld_type);
	printf("Meto2 Level Code         : %d\n", d->meto2_lvl_code);
	printf("Reserved #1              : %d %d %d %d\n", d->reserved1[0],
							   d->reserved1[1],
							   d->reserved1[2],
							   d->reserved1[3]);
	printf("Spare #1                 : %d\n", d->spare1);
	printf("Spare #2                 : %d %d %d %d %d %d %d\n",
							   d->spare2[0],
							   d->spare2[1],
							   d->spare2[2],
							   d->spare2[3],
							   d->spare2[4],
							   d->spare2[5],
							   d->spare2[6]);
	printf("Reserved #2              : %f %f %f %f\n", d->reserved2[0],
							   d->reserved2[1],
							   d->reserved2[2],
							   d->reserved2[3]);
*/
	printf("Datum Value              : %.1f\n", d->datum);
	printf("Packing Accuracy         : %.1f\n", d->pack_accuracy);
	printf("Level Value              : %.3f\n", d->level);

/* The following are not used */
/*
	printf("Reference Level Value    : %f\n", d->ref_level);
	printf("Level A-value            : %f\n", d->level_a);
	printf("Reference Level A-value  : %f\n", d->ref_level_a);
	printf("Latitude Pseudo N-pole   : %f\n", d->lat_ps_npol);
	printf("Longitude Pseudo N-pole  : %f\n", d->lon_ps_npol);
	printf("Grid Orientation         : %f\n", d->grid_orient);
*/
	printf("Latitude Zeroth Point    : %.2f\n", d->lat_0_row);
	printf("Latitude Interval        : %.2f\n", d->lat_interval);
	printf("Longitude Zeroth Point   : %.3f\n", d->lon_0_row);
	printf("Longitude Interval       : %.3f\n", d->lon_interval);
	printf("Missing Data Indicator   : %g\n", d->missing_val);
	printf("MKS Scaling Factor       : %.1f\n", d->scale);

        if ( d->field_code ==  1 )
		strcpy(prname, "HEIGHT (m)");
        else if ( d->field_code == 16 )
		strcpy(prname, "TEMP (K)");
        else if ( d->field_code == 48 || d->field_code == 56 )
		strcpy(prname, "ZONWIN_P (m/s)");
        else if ( d->field_code == 49 || d->field_code == 57 )
		strcpy(prname, "MERWIN_P (m/s)");
        else if ( d->field_code == 40 )
		strcpy(prname, "OMEGA (Pa/s)");
	else
		strcpy(prname, "Parameter ?");


	lat1st = (int)((nlat - d->lat_0_row)/d->lat_interval + 0.49999);
	if (lat1st < 1) lat1st = 1;

	latend = (int)((slat - d->lat_0_row)/d->lat_interval + 0.5);
	if (latend > d->num_rows) latend = d->num_rows;

	lon1st = (int)((wlon - d->lon_0_row)/d->lon_interval + 0.49999);
	if (fmod(lon1st, d->num_cols) != lonend)
	    lon1st = fmod(lon1st, d->num_cols);
	if (lon1st < 1) lon1st = 1;

	lonend = (int)((elon - d->lon_0_row)/d->lon_interval + 0.5);
	if (fmod(lonend, d->num_cols) != lon1st)
	    lonend = fmod(lonend, d->num_cols);
	if (lonend > d->num_cols) lonend = d->num_cols;


	printf("\n                      Gridded Data Values: %s\n", prname);
	printf("                    ----------------------------------------");
	printf("\n\nLongitude Array :\n");
	if (lon1st <= lonend) {
            for (i=lon1st-1, k=0; i<lonend; i++, k++) {
        	printf("%13.3f", d->lon_interval*(i+1) + d->lon_0_row);
		if ((k+1)%6 == 0) printf("\n");
            }
	}
	else {
            for (i=lon1st-1, k=0; i<d->num_cols; i++, k++) {
        	printf("%13.3f", d->lon_interval*(i+1) + d->lon_0_row);
		if ((k+1)%6 == 0) printf("\n");
            }
            for (i=0; i<lonend; i++, k++) {
        	printf("%13.3f", d->lon_interval*(i+1) + d->lon_0_row);
		if ((k+1)%6 == 0) printf("\n");
            }
	}
	if (k%6 != 0) printf("\n");
	
        for (j=lat1st-1; j<latend; j++)
        {
	    printf("\nLatitude : %7.3f\n", d->lat_interval*(j+1)+d->lat_0_row);

	    if (lon1st <= lonend) {
                for (i=lon1st-1, k=0; i<lonend; i++, k++) {
        	    printf("%13.5e",d->scale*d->data[j*d->num_cols+i]+d->datum);
		    if ((k+1)%6 == 0) printf("\n");
                }
	    }
	    else {
                for (i=lon1st-1, k=0; i<d->num_cols; i++, k++) {
        	    printf("%13.3e",d->scale*d->data[j*d->num_cols+i]+d->datum);
		    if ((k+1)%6 == 0) printf("\n");
                }
                for (i=0; i<lonend; i++, k++) {
        	    printf("%13.5e",d->scale*d->data[j*d->num_cols+i]+d->datum);
		    if ((k+1)%6 == 0) printf("\n");
                }
	    }
	    if (k%6 != 0) printf("\n");

        }
    }

    return;
}


/************************************************************************
 * ParseArgs -- parse arguments from the command line in main program
 ************************************************************************/

int ParseArgs(int argc,                 /* input number of arguments */
              char *argv[],             /* input list of arguments */
              char *filename,           /* name of input HDF file */
              float *slat,	        /* southernmost latitude */
              float *nlat,	        /* northernmost latitude */
              float *wlon,	        /* westernmost longitude */
              float *elon,	        /* easternmost longitude */
              float *zmin,	        /* minimum altitude */
              float *zmax,	        /* maximum altitude */
              char *uparam,	        /* UKMO parameter */
              float *wavmin,	        /* minimum wavelength */
              float *wavmax,	        /* maximum wavelength */
              char *swap		/* flag to swap bytes */
)
{
    int i, n;                   /* increment counters */
    char bflag=0;               /* flag indicates if one should byteswap */
    char fflag=0;               /* flag indicates if filename exists */
    char sflag=0;               /* flag indicates if spatial values set */
    char vflag=0;               /* flag indicates if vertical values set */
    char wflag=0;               /* flag indicates if wavelength values set */
    char uflag=0;               /* flag indicates if UKMO parameter set */


    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-')
        {
            switch(argv[i][1])
            {
                case 'b':       /* byteswap option */
		    bflag = 1;
                    break;
                case 's':       /* spatial subset option */
		    if (sflag) return -1;
		    i++;
		    n = GetLatLon(argc-i, &argv[i], slat, nlat, wlon, elon);
		    if (n != 2 && n != 4) return -1;
		    i+=(n-1);
		    sflag = 1;
                    break;
                case 'v':       /* vertical subset option */
		    if (vflag) return -1;
		    i++;
		    n = GetVertical(argc-i, &argv[i], zmin, zmax);
		    if (n != 2) return -1;
		    i+=(n-1);
		    vflag = 1;
                    break;
                case 'w':       /* wavelength subset option */
		    if (wflag) return -1;
		    i++;
		    n = GetWavelength(argc-i, &argv[i], wavmin, wavmax);
		    if (n != 2) return -1;
		    i+=(n-1);
		    wflag = 1;
                    break;
                case 'u':       /* UKMO parameter option */
		    if (uflag) return -1;
		    i++;
		    *uparam = *argv[i];
		    uflag = 1;
                    break;
                case '?':       /* show usage */
		    Usage();
		    exit (0);
                    break;
                default:
		    return -1;
            }
        }
        else if (! fflag)
        {
            fflag = 1;
            strcpy(filename, argv[i]);
        }
        else return -1;
    }

    *swap = bflag;
    if (! fflag) return -1;
    else return 0;
}


int GetRecSize(const char *file_name, char *type, char *subtype, char *lvlsrc)
{
    char rec[MAXRECSIZE];		/* array to hold records */
    char *ptr;				/* string pointer */
    int rec_size;			/* size of data record */
    int key;				/* level 3AL index key size */
    int i;				/* increment counter */
    FILE *fp;				/* input file pointer */

    if ((fp = fopen(file_name, "rb")) == NULL) {
        fprintf(stderr, "Error: unable to open %s\n", file_name);
        return -1;
    }

    fread(rec, 1, MAXRECSIZE, fp);
    if ((ptr = strstr(rec, "CCSD1Z")) == NULL) {
	fclose(fp);
	fprintf(stderr, "Error: %s is not a valid UARS data file\n", file_name);
	return -1;
    }
    key = ptr - rec;

    if ( (ptr = strstr(rec, "nmc")) - rec == 132 ||
         (ptr = strstr(rec, "NMC")) - rec == 132 )
    {
	strcpy(type, "CORR");
	strcpy(subtype, "?");
        strcpy(lvlsrc, "NMC");

        rec_size = 16964;	/* Fixed */
	ptr = &rec[rec_size];
    }
    else if ( (ptr = strstr(rec, "ukmo")) - rec == 130 ||
              (ptr = strstr(rec, "UKMO")) - rec == 130 )
    {
	strcpy(type, "CORR");
	strcpy(subtype, "ASSIM");
	strcpy(lvlsrc, "UKMO");

        rec_size = 28288;	/* Fixed */
    }
    else
    {
        for (i=0, rec_size=0; i<MAXRECSIZE; i+=4)
        {
            if ((ptr = strstr(&rec[i], "UARS")) != NULL) {
	        rec_size = ptr - rec - key;
	        break;
            }
        }

        if (rec_size == 0) {
	    fprintf(stderr, "Error: unable to determine record size\n");
	    return -1;
        }

	ptr += 6;  strncpy(type, ptr, 12);
	ptr += 12; strncpy(subtype, ptr, 12);
        ptr += 87; strncpy(lvlsrc, ptr, 3);
    }

    return rec_size;
}


int GetVertical(int argc, char *argv[], float *zmin, float *zmax)
{
    int n=0;
    char *ptr;
    double dval;

    if (argc < 2) return -1;

    dval = strtod(argv[0], &ptr);
    if (ptr == argv[0]) return -1;
    *zmin = (float) dval; n++;

    dval = strtod(argv[1], &ptr);
    if (ptr == argv[1]) return -1;
    *zmax = (float) dval; n++;

    return n;
}


int GetWavelength(int argc, char *argv[], float *wavmin, float *wavmax)
{
    int n=0;
    char *ptr;
    double dval;

    if (argc < 2) return -1;

    dval = strtod(argv[0], &ptr);
    if (ptr == argv[0]) return -1;
    *wavmin = (float) dval; n++;

    dval = strtod(argv[1], &ptr);
    if (ptr == argv[1]) return -1;
    *wavmax = (float) dval; n++;

    return n;
}


int GetLatLon(int argc, char *argv[], float *slat, float *nlat,
	      float *wlon, float *elon)
{
    int n=0;
    char *ptr;
    double dval;

    if (argc < 2) return -1;

    dval = strtod(argv[0], &ptr);
    if (ptr == argv[0]) return -1;
    *slat = (float) dval; n++;

    dval = strtod(argv[1], &ptr);
    if (ptr == argv[1]) return -1;
    *nlat = (float) dval; n++;

    if (argc >= 4)
    {
	dval = strtod(argv[2], &ptr);
	if (ptr != argv[2])
	{
	    *wlon = (float) dval; n++;

	    dval = strtod(argv[3], &ptr);
	    if (ptr == argv[3]) return -1;
	    *elon = (float) dval; n++;
	}
    }

    return n;
}


void w3fb04 (float alat, float along, float xmeshl, float orient, 
     float *xi, float *xj) {

/*
C
C SUBPROGRAM: W3FB04         LATITUDE, LONGITUDE TO GRID COORDINATES
C   AUTHOR: MCDONELL,J.      ORG: W345       DATE: 90-06-04
C
C ABSTRACT: CONVERTS THE COORDINATES OF A LOCATION ON EARTH FROM THE
C   NATURAL COORDINATE SYSTEM OF LATITUDE/LONGITUDE TO THE GRID (I,J)
C   COORDINATE SYSTEM OVERLAID ON A POLAR STEREOGRAPHIC MAP PRO-
C   JECTION TRUE AT 60 DEGREES N OR S LATITUDE. W3FB04 IS THE REVERSE
C   OF W3FB05.
C
C PROGRAM HISTORY LOG:
C   77-05-01  J. MCDONELL 
C   89-01-10  R.E.JONES   CONVERT TO MICROSOFT FORTRAN 4.1
C   90-06-04  R.E.JONES   CONVERT TO SUN FORTRAN 1.3
C   93-01-26  B. Doty     converted to C
C
C USAGE:  CALL W3FB04 (ALAT, ALONG, XMESHL, ORIENT, XI, XJ)
C
C   INPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     ALAT   ARG LIST  LATITUDE IN DEGREES (<0 IF SH)
C     ALONG  ARG LIST  WEST LONGITUDE IN DEGREES
C     XMESHL ARG LIST  MESH LENGTH OF GRID IN KM AT 60 DEG LAT(<0 IF SH)
C                   (190.5 LFM GRID, 381.0 NH PE GRID,-381.0 SH PE GRID)
C     ORIENT ARG LIST  ORIENTATION WEST LONGITUDE OF THE GRID
C                   (105.0 LFM GRID, 80.0 NH PE GRID, 260.0 SH PE GRID)
C
C   OUTPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     XI     ARG LIST  I OF THE POINT RELATIVE TO NORTH OR SOUTH POLE
C     XJ     ARG LIST  J OF THE POINT RELATIVE TO NORTH OR SOUTH POLE
C
C   SUBPROGRAMS CALLED:
C     NAMES                                                   LIBRARY
C     ------------------------------------------------------- --------
C     COS SIN                                                 SYSLIB
C
C   REMARKS: ALL PARAMETERS IN THE CALLING STATEMENT MUST BE
C     REAL. THE RANGE OF ALLOWABLE LATITUDES IS FROM A POLE TO
C     30 DEGREES INTO THE OPPOSITE HEMISPHERE.
C     THE GRID USED IN THIS SUBROUTINE HAS ITS ORIGIN (I=0,J=0)
C     AT THE POLE IN EITHER HEMISPHERE, SO IF THE USER'S GRID HAS ITS
C     ORIGIN AT A POINT OTHER THAN THE POLE, A TRANSLATION IS NEEDED
C     TO GET I AND J. THE GRIDLINES OF I=CONSTANT ARE PARALLEL TO A
C     LONGITUDE DESIGNATED BY THE USER. THE EARTH'S RADIUS IS TAKEN
C     TO BE 6371.2 KM.
C
C ATTRIBUTES:
C   LANGUAGE: SUN FORTRAN 1.4
C   MACHINE:  SUN SPARCSTATION 1+
C*/

static float radpd = 0.01745329;
static float earthr = 6371.2;

float re,xlat,wlong,r;

      re    = (earthr * 1.86603) / xmeshl;
      xlat =  alat * radpd;
 
      if (xmeshl>0.0) { 
        wlong = (along + 180.0 - orient) * radpd;
        r     = (re * cos(xlat)) / (1.0 + sin(xlat));
        *xi    = r * sin(wlong);
        *xj    = r * cos(wlong);
 
      } else {
 
        re    = -re;
        xlat =  -xlat;
        wlong = (along - orient) * radpd; 
        r     = (re * cos(xlat)) / (1.0+ sin(xlat));
        *xi   =  r * sin(wlong);
        *xj   = -r * cos(wlong);
      }
}



void w3fb05 (float xi, float xj, float xmeshl, float orient, 
     float *alat, float *along) {
/*
C
C SUBPROGRAM: W3FB05         GRID COORDINATES TO LATITUDE, LONGITUDE
C   AUTHOR: JONES,R.E.       ORG: W345       DATE: 86-07-17
C
C ABSTRACT: CONVERTS THE COORDINATES OF A LOCATION FROM THE GRID(I,J]
C   COORDINATE SYSTEM OVERLAID ON THE POLAR STEREOGRAPHIC MAP PROJEC-
C   TION TRUE AT 60 DEGREES N OR S LATITUDE TO THE NATURAL COORDINATE
C   SYSTEM OF LATITUDE/LONGITUDE ON THE EARTH. W3FB05 IS THE REVERSE
C   OF W3FB04.
C
C PROGRAM HISTORY LOG:
C   86-07-17  R.E.JONES
C   88-06-20  R.E.JONES   CHANGE TO MICROSOFT FORTRAN 4.10
C   89-03-29  R.E.JONES   CHANGE TO VAX-11 FORTRAN 
C
C USAGE:  CALL W3FB05 (XI, XJ, XMESHL, ORIENT, ALAT, ALONG]
C
C   INPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     XI     ARG LIST  I OF THE POINT RELATIVE TO THE NORTH OR S. POLE
C     XJ     ARG LIST  J OF THE POINT RELATIVE TO THE NORTH OR S. POLE
C     XMESHL ARG LIST  MESH LENGTH OF GRID IN KM AT 60 DEGREES(<0 IF SH]
C                   (190.5 LFM GRID, 381.0 NH PE GRID,-381.0 SH PE GRID]
C     ORIENT ARG LIST  ORIENTATION WEST LONGITUDE OF THE GRID
C                    (105.0 LFM GRID, 80.0 NH PE GRID, 260.0 SH PE GRID]
C
C   OUTPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     ALAT   ARG LIST  LATITUDE IN DEGREES  (<0 IF SH]
C     ALONG  ARG LIST  WEST LONGITUDE IN DEGREES
C
C   SUBPROGRAMS CALLED:
C     NAMES                                                   LIBRARY
C     ------------------------------------------------------- --------
C     ASIN   ATAN2                                            SYSLIB
C
C   REMARKS: ALL PARAMETERS IN THE CALLING STATEMENT MUST BE
C     REAL. THE RANGE OF ALLOWABLE LATITUDES IS FROM A POLE TO
C     30 DEGREES INTO THE OPPOSITE HEMISPHERE.
C     THE GRID USED IN THIS SUBROUTINE HAS ITS ORIGIN (I=0,J=0]
C     AT THE POLE, SO IF THE USER'S GRID HAS ITS ORIGIN AT A POINT
C     OTHER THAN A POLE, A TRANSLATION IS REQUIRED TO GET I AND J FOR
C     INPUT INTO W3FB05. THE SUBROUTINE GRID IS ORIENTED SO THAT
C     GRIDLINES OF I=CONSTANT ARE PARALLEL TO A WEST LONGITUDE SUP-
C     PLIED BY THE USER. THE EARTH'S RADIUS IS TAKEN TO BE 6371.2 KM.
C
C   LANGUAGE: VAX-11 FORTRAN
C   MACHINE:  VAX-11 730, 750, 780, 8600, ETC.
C
C   FOR DETAILS OBTAIN WRITEUP FROM NMC/AD/SAB
C */

  float degprd = 57.2957795;
  float earthr = 6371.2;
  float gi2,r2,angle;

  gi2 = pow(((1.86603 * earthr)/(xmeshl)), 2);
  r2 = xi*xi + xj*xj;
  if (r2 == 0.0)
  {
    *along = 0.0;
    *alat = 90.0;
    if (xmeshl < 0.0) *alat = -1 * *alat;
  } 
  else
  {
   *alat = asin((gi2 - r2) / (gi2 + r2)) * degprd;
   angle = degprd * atan2(xj,xi);
   if (angle < 0.0) angle += 360.0;

   if (xmeshl >= 0.0)
     *along = 270.0 + orient - angle;
   else
   {
     *along = angle + orient - 270.0;
     *alat = -1 * *alat;
   }
   if (*along < 0.0) *along += 360.0;
   if (*along >= 360.0) *along -= 360.0;
  }

  return;
}


/************************************************************************
 * Usage -- tell the user what to do.
 ************************************************************************/

void Usage(void)
{
    fprintf(stderr, "Usage:  readuars [-?] [-b] "
		    "[-s <slat> <nlat> [<wlon> <elon>]] "
		    "[-v <zmin> <zmax>] [-u <uparam>] <filename>\n\n");
    fprintf(stderr, "Options/Arguments:\n\n");
    fprintf(stderr, "\t-? -- help: prints this page\n\n");
    fprintf(stderr, "\t-b -- swap bytes for integers and floats\n\n");
    fprintf(stderr, "\t-s <slat> <nlat> <wlon> <elon> "
		    "-- spatial region to subset by where\n"
		    "\t\t<slat> = southernmost latitude (-90 to +90)\n"
		    "\t\t<nlat> = northernmost latitude (-90 to +90)\n"
		    "\t\t<wlon> = westernmost longitude (-180 to +180 or 0 to 360)\n"
		    "\t\t<elon> = easternmost longitude (-180 to +180 or 0 to 360)\n\n");
    fprintf(stderr, "\t-v <zmin> <zmax> "
		    "-- vertical range to subset by where\n"
		    "\t\t<zmin> = minimum altitude (1000 to 1.0e-6 hPa or 0 to 400 km)\n"
		    "\t\t<zmax> = maximum altitude (1000 to 1.0e-6 hPa or 0 to 400 km)\n\n");
    fprintf(stderr, "\t-u <uparam> "
		    "-- UKMO parameter where options include h=HEIGHT, t=TEMP,\n"
		    "\t\t\tm=MERWIN_P, z=ZONWIN_P, o=OMEGA, or a=ALL (default)\n\n");
    fprintf(stderr, "\t<filename> -- name of the input file (META or PROD)\n");

  return;
}


/* ****************************************************************************
  Function : byteswap

  Description : This function performs byte swapping for 2 and 4 byte integers, 
		as well as 4 and 8 byte floats.
**************************************************************************** */

int byteswap(char *ptr, int len)
{
	char tmp;
	
	switch(len)
	{
	default:
		break;
	case 2:
		tmp=*ptr; *ptr=*(ptr+1);*(ptr+1)=tmp;
		break;
	case 4:
		tmp=*ptr; *ptr=*(ptr+3);*(ptr+3)=tmp;
		tmp=*(ptr+1); *(ptr+1)=*(ptr+2);*(ptr+2)=tmp;
		break;
	case 8:
		tmp=*ptr; *ptr=*(ptr+7);*(ptr+7)=tmp;
		tmp=*(ptr+1); *(ptr+1)=*(ptr+6);*(ptr+6)=tmp;
		tmp=*(ptr+2); *(ptr+2)=*(ptr+5);*(ptr+5)=tmp;
		tmp=*(ptr+3); *(ptr+3)=*(ptr+4);*(ptr+4)=tmp;
		break;
	}

	return(0);
}

