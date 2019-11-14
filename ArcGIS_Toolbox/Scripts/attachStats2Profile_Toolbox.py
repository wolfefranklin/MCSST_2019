'''
Attaches statistics from a CSV file to the profile lines (statistics calculated within a Jupiter Notebook)
Run as:
    1. Stand alone script with raw input by user of all required inputs.
    2. From an ArcGIS toolbox (created with ArcGIS v10.5.1)
    3. As a test by seting inputs within the script
Notes:
    All the data inputs should include the full path
    Statistics file must be created using the input profile line data
    Stats file must be in right format:
        There are 20 stats values and 2 ages (if provided) exported for each profile
        Stats for each profile is one of the columns in the stats CSV
        The name of the stats column must be composed from "profile_" and the profile ID from the profile lines dataset
        The first column describes statistics stored in each row
        The last column can hold Cumulative value in which case is ignored
Required:
    Profile lines in a feature class or shapefile (with a field holding an unique profile ID)
    Statistics CSV with columns holding stats for profile lines (column name = "profile_" + profileID)
Output:
    Stats data is attached to the original profile line dataset
    Stats data is attched to each profile for which the statistics is calculated
    For profiles in the profile line datsets which don't have statistics calculated, the stats fields will hold 0 or NULL
'''

# Import system modules
import arcpy
import os, sys

#-------------------------------------------------------------------------------------

#Read input parameters from script tool
inProfile = arcpy.GetParameterAsText(0)
profIDField = arcpy.GetParameterAsText(1)
inStatCSV = arcpy.GetParameterAsText(2)

#------------------------------------------------------------------------------------------

# Set enviroments
arcpy.env.overwriteOutput = True

# set dictionary of field types
fieldTypeDict = {"Integer":"LONG", "SmallInteger":"SHORT", "String":"TEXT"}

# get list of the field names in the profile lines field
fldList = [f.name for f in arcpy.ListFields(inProfile)]

# get the type of profile ID field
for f_pair in [(f.name, f.type) for f in arcpy.ListFields(inProfile)]:
    if f_pair[0] == profIDField:
        profIDFieldType = fieldTypeDict[f_pair[1]]

# set dictionary of field names for each row in the stats CSV
fieldName = {1:"SMean", 2:"SMedian", 3:"SMin", 4:"SMax", 5:"S2STD", 6:"SRMean", 7:"SRMedian", 8:"SRMin", 9:"SRMax", 10:"SR2STD",
             11:"SRHMean", 12:"SRHMedian", 13:"SRHMin", 14:"SRHMax", 15:"SRH2STD", 16:"SRVMean", 17:"SRVMedian", 18:"SRVMin",
             19:"SRVMax", 20:"SRV2STD", 21:"AgeMin", 22:"AgeMax"}

# Adding stats fields (delete field first if exists)
arcpy.AddMessage ("Adding statistics fields to the profile lines")
for f_cnt in range(1,23):
    if fieldName[f_cnt] in fldList:
        arcpy.DeleteField_management(inProfile, fieldName[f_cnt])
    arcpy.AddField_management(inProfile, fieldName[f_cnt], "DOUBLE")

# Read statistics from file and attach to the profile line
line_num = 0
lines = open(inStatCSV, "r").readlines() # list of lines in the CSV file
numRows = (len(lines)) # number of outputs to be attached
num_sec = len(lines[0].strip("\n").split(",")) # get number of columns

# read profile by profile  (column by column) - ignore the first column (stats descriptions) and last (if Cumulative)
for i in range(1, num_sec):

    profFileName = lines[0].strip("\n").split(",")[i]   # read profile file name from file
    profName = profFileName.replace("profile_", "") # strip profile from the name (added when the createProfileCSV.py was run)

    if profName != "Cumulative":

        # Attaching statistics to section lines
        arcpy.AddMessage ("Profile: " + profName)

        # create filtering clause
        if profIDFieldType in ["LONG", "SHORT"]:
            where_clause = '"' + profIDField + '" = ' + profName  # set query for serching prof lines if ID is integer
        elif profIDFieldType in ["TEXT"]:
            where_clause = '"' + profIDField + '" = ' + "'" + profName + "'"  # set query for serching prof lines if ID is string
        else:
            sys.exit("!!!Wrong field type for Profile ID!!!!")

        # selecting the profile
        arcpy.MakeFeatureLayer_management(inProfile, "sec", where_clause) # select section
        # check if profile exists
        if int(arcpy.GetCount_management("sec").getOutput(0)) > 0:
            arcpy.AddMessage (" Attaching statistics from file")
            for f_cnt in range(1,numRows):
                arcpy.CalculateField_management("sec", fieldName[f_cnt], lines[f_cnt].strip("\n").split(",")[i])

        else:
            arcpy.AddMessage (" Profile " + profName + " doesn't exist")

arcpy.Delete_management("sec")


