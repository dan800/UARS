# output from chatgpt on providing readuars C classes and defs.

class Sfdu:
    def __init__(self):
        self.key = ""  # VMS record key (for 3AL/3LP files only)
        self.Tz = ""  # SFDU type (Tz) field
        self.Lz = ""  # SFDU length (Lz) field
        self.Ti = ""  # SFDU type (Ti) field
        self.Li = ""  # SFDU length (Li) field


class HeadL3A:
    def __init__(self):
        self.key = ""  # VMS record key (for 3AL/3LP files only)
        self.Sat_ID = ""  # Satellite identifier (always set to UARS)
        self.Rec_type = ""  # Record type (always set to '1')
        self.Inst_ID = ""  # Instrument identifier (CLAES, HALOE, etc.)
        self.Subtype = ""  # Data subtype or species (H2O, O3, TEMP, etc.)
        self.Format_ver = ""  # Format version number (always set to '1')
        self.Rec_count = ""  # Physical record count (always set to '1', 1st rec)
        self.N_cont = ""  # Number of continuation label records (0 for UARS)
        self.N_recs = ""  # Number of physical records in file
        self.Create = ""  # File creation time in VAX/VMS ASCII format
        self.Y_1st = ""  # Year for first data record
        self.D_1st = ""  # Day of year for first data record
        self.Ms_1st = ""  # Milliseconds of day for first data record
        self.Y_last = ""  # Year for last data record
        self.D_last = ""  # Day of year for last data record
        self.Ms_last = ""  # Milliseconds of day for last data record
        self.Level = ""  # Data process level
        self.Day = ""  # UARS day number
        self.N_pts = ""  # Number of data points or 32-bit words per record
        self.Bs_index = ""  # Base index of data points (not used in 3LP/3TP/3S)
        self.Rec_len = ""  # Record length in bytes
        self.Min_lat = ""  # Minimum Latitude (3AL/3LP only)
        self.Max_lat = ""  # Maximum Latitude (3AL/3LP only)
        self.Version = ""  # CCB version number
        self.Cycle = ""  # File cycle number (wrong, see META file instead)
        self.Virt_flag = ""  # Virtual file flag (always <blank> for UARS)
        self.TV_file = ""  # Total number of time/version entries in file (= 0)
        self.TV_rec = ""  # Number of time/version entries in record (= 0)


class HeadCORR:
    def __init__(self):
        self.SFDU_Tz = ""  # SFDU type (Tz) field
        self.SFDU_Lz = ""  # SFDU length (Lz) field
        self.SFDU_Ti = ""  # SFDU type (Ti) field
        self.SFDU_Li = ""  # SFDU length (Li) field
        self.SFDU_Vi = ""  # SFDU variable (Vi) field
        self.Project = ""  # Project name (always set to UARS)
        self.UARS_PI = ""  # UARS Principal Investigators name
        self.UARS_CMI = ""  # UARS Correlative Measurement Investigators name
        self.Class = ""  # Correlative data class
        self.Inst_ID = ""  # Instrument identifier (CLAES, HALOE, etc.)
        self.Station_ID = ""  # Observing station identifier
        self.File_type = ""  # Correlative data file type
        self.Start_time = ""  # Start time of data in file
        self.Stop_time = ""  # Start time of data in file
        self.Max_lat = ""  # Maximum latitude of data in file
        self.Min_lat = "" 

class DataL3A:
    def __init__(self):
        self.key = ""  # VMS record key (for 3AL/3LP files only)
        self.Sat_ID = ""  # Satellite Identifier (always set to UARS)
        self.Rec_type = ""  # Record type (always set to '1')
        self.Inst_ID = ""  # Instrument Identifier (HALOE, ISAMS, etc.)
        self.Rec_count = ""  # Physical record count (>= 1)
        self.Spare = ""  # Spare - unused
        self.Total = 0  # Total number of data points in the record
        self.Actual = 0  # Number of actual data points in the record
        self.St_index = 0  # Starting index of first actual point
        self.Time = [0, 0]  # Record time in UDTF format [YYDDD,MMMMMMMM]
        self.Lat = 0.0  # Geodetic latitude (-90 to +90)
        self.Lon = 0.0  # Geodetic longitude (0 to 360)
        self.LST = 0.0  # Local solar time (0 to 24)
        self.SZA = 0.0  # Solar zenith angle (0 to 180)
        self.data = []  # pointer to array of data values
        self.quality = []  # pointer to array of standard deviations


MAXRECSIZE = 8192  # Define the maximum record size
FILL = float("nan")  # Define the fill value for missing data

def ReadL3A(file_name):
    rec_size = MAXRECSIZE
    swap = False

    rec = bytearray(rec_size)  # Array to hold records
    temp = bytearray(rec_size)  # Temporary array for manipulation
    sfdu = bytearray(60 + 1)  # Array to hold SFDU record
    header = bytearray(176 + 1)  # Array to hold header/label record

    # String pointer
    ptr = None

    # Increment counters
    i = 0
    j = 0

    # Size of index key (L3AL only)
    key = 0

    # Total number of data records
    nrecs = 0

    # Number of data array elements
    npts = 0

    # Pointer for 32-bit integers
    i32ptr = None

    # Year, month, day, msec
    year = 0
    month = 0
    day = 0
    msec = 0

    # Pointer for 32-bit floats
    f32ptr = None

    # Pointer for data value array
    data = None
    dptr = None

    # Pointer for data quality array
    qual = None
    qptr = None

    # Structure for SFDU record
    sr = Sfdu()

    # Structure for header record
    hr = HeadL3A()

    # Structure for data record
    dr = None
    dtmp = DataL3A()

    with open(file_name, "rb") as fp:
        # The 0th record contains the SFDU information
        fp.readinto(rec)
        ptr = rec.find(b"CCSD1Z")
        if ptr == -1:
            raise ValueError(f"Error: {file_name} is not a valid UARS data file")
        key = ptr
        sfdu[20 - key: 20 - key + 40 + key] = rec[:40 + key]

        # The 1st record is the header or label
        fp.readinto(rec)
        header[20 - key: 20 - key + 125 + key] = rec[:125 + key]
        if hr.Level.startswith("3AL") or hr.Level.startswith("3LP"):
            header[125 + 20: 125 + 20 + 29] = rec[key + 125:]

        # Get the number of data records in the file
        temp[:8] = hr.N_recs.encode()
        nrecs = int(temp[:8]) - 1  # Minus 1 because the 1st record is the header

        # Get the data array size (number of levels) in the file
        temp[:4] = hr.N_pts.encode()
        npts = int(temp[:4])

        # Allocate memory for the data records
        dr = [DataL3A() for _ in range(nrecs)]

        # Allocate memory for the data


def print_sfd(sfd):
    print("-------------------------------")
    print("    SFDU Descriptor")
    print("-------------------------------")

    # Uncomment if you want to see this value.
    if len(sfd.key) != 0:
        print("Record key : {}".format(sfd.key))

    print("T[z] field : {}".format(sfd.Tz))
    print("L[z] field : {}".format(sfd.Lz))
    print("T[i] field : {}".format(sfd.Ti))
    print("L[i] field : {}\n".format(sfd.Li))

    return

def print_head_l3a(hr):
    print("-------------------------------")
    print("    Header Record")
    print("-------------------------------")

    # Uncomment if you want to see this value.
    if len(hr.key) != 0:
        print("Record key : {}".format(hr.key))

    print("Satellite ID  : {}".format(hr.Sat_ID))
    # Uncomment if you want to see this value.
    print("Record type   : {}".format(hr.Rec_type))

    print("Instrument ID : {}".format(hr.Inst_ID))
    print("Data subtype or species : {}".format(hr.Subtype))
    # Uncomment if you want to see these values.
    print("Format version number   : {}".format(hr.Format_ver))
    print("Physical record count   : {}".format(hr.Rec_count))
    print("Number of continuation records for file label : {}".format(hr.N_cont))
    print("Number of physical records in file : {}".format(hr.N_recs))

    print("File Creation Time : {}".format(hr.Create))
    print("Year of 1st data record  : {}".format(hr.Y_1st))
    print("Day of 1st data record   : {}".format(hr.D_1st))
    print("Milliseconds of 1st data record  : {}".format(hr.Ms_1st))
    print("Year of last data record : {}".format(hr.Y_last))
    print("Day of last data record  : {}".format(hr.D_last))
    print("Milliseconds of last data record : {}".format(hr.Ms_last))
    print("Data level : {}".format(hr.Level))
    print("UARS day number : {}".format(hr.Day))
    # Uncomment if you want to see these values.
    print("Number of data points per record : {}".format(hr.N_pts))
    print("Base index of data point value : {}".format(hr.Bs_index))
    print("Record length in bytes : {}".format(hr.Rec_len))

    if hr.Level.startswith("3AL") or hr.Level.startswith("3LP"):
        print("Minimum latitude for records in file : {}".format(hr.Min_lat))
        print("Maximum latitude for records in file : {}".format(hr.Max_lat))
    print("CCB version number : {}".format(hr.Version))
    print("File cycle number      : {}".format(hr.Cycle))
    # Uncomment if you want to see these values.
    print("Virtual file flag : {}".format(hr.Virt_flag))
    print("Total number of time/version entries in file : {}".format(hr.TV_file))
    print("Number of time/version entries in record : {}".format(hr.TV_rec))

    return

