'''
Creates points along defined profile lines at equal intervals, extracts elevation and geology info and exports data to a separate CSV file for each profile
Run as:
    From an ArcGIS toolbox (created with ArcGIS v10.5.1)
Notes:
    All the data sources should include the full path.
    Input field names have to exist within the relevant dataset
    Tool requires ArcGIS Spatial Analyst licence
Required:
    The working folder (full path) where processing is done and the output folder is created
    Profile lines in a feature class or shapefile (with a field holding an unique profile ID)
    DEM
    geology layer if exists (holding info on geology unit name and min and max unit age)
    point distance (distance of points along a profile line ued to extract info)
Outputs:
    CSV file for each profile line
    CSV file names are composed from "profile_" and the profile line ID (e.g. profile_1.csv)
    The output folder name is composed from "Outputs_" and the profile file name (e.g. Outputs_ProfileLines)
'''

# Import system modules
import arcpy
import os, sys
from arcpy.sa import *
import uuid

# licences
arcpy.CheckOutExtension('Spatial')

#-------------------------------------------------------------------------------------

#Read input parameters from script tool
workDir = arcpy.GetParameterAsText(0)
inDEM = arcpy.GetParameterAsText(1)
geolExist = arcpy.GetParameterAsText(2)
inGeol = arcpy.GetParameterAsText(3)
ageMinField = arcpy.GetParameterAsText(4)
ageMaxField = arcpy.GetParameterAsText(5)
geolField = arcpy.GetParameterAsText(6)
inProfile = arcpy.GetParameterAsText(7)
profIDField = arcpy.GetParameterAsText(8)
pointDist = int(arcpy.GetParameterAsText(9))

#-------------------------------------------------------------------------------------

# Outputs
profName = os.path.basename(inProfile.replace(".shp", ""))
outCSVFolder = os.path.join(workDir, "Outputs_" + profName)
outProfPt = profName + "_pt"
distField = "Distance"

# set dictionary of field types
fieldTypeDict = {"Integer":"LONG", "SmallInteger":"SHORT", "String":"TEXT"}

# Define function for createing points along line

def points_along_line(line_lyr, pnt_layer, prof_id, pnt_dist, dist_f):
    """
    line_lyr (feature layer) - Single part line
    pnt_layer (feature layer) - Path to point feature class
    prof_id (attribute field) - Name of the id field (to be  added to points)
    pnt_dist (integer) - Interval distance in map units between points
    dist_f (attribute field) - Field to store distance along line
    """

    search_cursor = arcpy.da.SearchCursor(line_lyr, ['SHAPE@', prof_id])
    insert_cursor = arcpy.da.InsertCursor(pnt_layer, ['SHAPE@', prof_id, dist_f])

    for row in search_cursor:
        distance = 0
        for dist in range(0, int(row[0].length) + 1, pnt_dist):
            point = row[0].positionAlongLine(dist).firstPoint
            point_id = row[1]
            insert_cursor.insertRow([point, point_id, distance])
            distance = distance + pointDist
    del search_cursor
    del insert_cursor

# Create output folders if they don't exist
if not os.path.exists(outCSVFolder):
    arcpy.AddMessage ("Creating output folder for CSV files")
    os.makedirs(outCSVFolder)

# Create working GDB
arcpy.AddMessage ("Creating temporary database")
gdbName = "temp_" + str(uuid.uuid4()) + ".gdb"
arcpy.CreateFileGDB_management(workDir, gdbName)
workGDB = os.path.join(workDir, gdbName)

# Set enviroments
arcpy.env.workspace = workGDB
arcpy.env.overwriteOutput = True

# Get the field type for ProfileID field
for f_pair in [(f.name, f.type) for f in arcpy.ListFields(inProfile)]:
    if f_pair[0] == profIDField:
        profIDFieldType = fieldTypeDict[f_pair[1]]

# Create point layer
arcpy.AddMessage ("Creating point layer")
if arcpy.Exists(outProfPt):
    arcpy.Delete_management(outProfPt)
arcpy.CreateFeatureclass_management(workGDB, outProfPt, 'POINT', '#', '#', '#', inProfile)
arcpy.AddField_management(outProfPt, profIDField, profIDFieldType)
arcpy.AddField_management(outProfPt, distField, "DOUBLE")

arcpy.AddMessage ("Creating points along profile lines")
points_along_line(inProfile, outProfPt, profIDField, pointDist, distField)

# Extracting elevation
arcpy.AddMessage ("Extracting elevation")
ExtractValuesToPoints(outProfPt, inDEM, "temp_elev")
arcpy.AddMessage (" Attaching elevation to profile points")
arcpy.JoinField_management(outProfPt, "OBJECTID", "temp_elev", "OBJECTID", "RASTERVALU")
arcpy.AddMessage (" Renaming field")
arcpy.AlterField_management(outProfPt, "RASTERVALU", "Elevation", "#", "#", "#", "#", "True")

# If geology layer is provided
if geolExist == "true":
    # Attaching geology age
    arcpy.AddMessage ("Attaching geology info")
    # delete fields if exist
    for fld in [geolField, ageMinField, ageMaxField]:
        if fld in [f.name for f in arcpy.ListFields(outProfPt)]:
            arcpy.AddMessage (" Deleting " + fld + " field")
            arcpy.DeleteField_management(outProfPt, fld)
    arcpy.AddMessage (" Performing Identity")
    arcpy.Identity_analysis(outProfPt, inGeol, "temp_" + fld)
    arcpy.AddMessage (" Joining geology info to summary points")
    arcpy.JoinField_management(outProfPt, "OBJECTID", "temp_" + fld, "FID_" + outProfPt, [geolField, ageMinField, ageMaxField])

# Writing profile files
arcpy.AddMessage ("Writing files")

# create a header list for the output
headerList = ["PointID", "ProfID", "Distance", "Geology", "Age_Min", "Age_Max", "Elevation"]
fieldList = [profIDField, distField, geolField, ageMinField, ageMaxField, "Elevation"]
if geolExist == "false":  # no geology exists
    headerList = ["PointID", "ProfID", "Distance", "Elevation"]
    fieldList = [profIDField, distField, "Elevation"]

with arcpy.da.SearchCursor(inProfile, [profIDField]) as line_cursor:
    for line_row in line_cursor:
        profID = line_row[0] # get profile id
        arcpy.AddMessage (" Profile: " + str(profID))
        outFile = os.path.join(outCSVFolder, "profile_" + str(profID) + ".csv") # set output file
        out_f = open(outFile, "w") # open file for writing

        # writing header
        header = ",".join(headerList) + "\n"
        out_f.write(header)

        if profIDFieldType in ["LONG", "SHORT"]:
            where_clause = '"' + profIDField + '" = ' + str(profID)  # set query for serching points
        elif profIDFieldType in ["TEXT"]:
            where_clause = '"' + profIDField + '" = ' + "'" + str(profID) + "'"  # set query for serching points
        else:
            sys.exit("!!!Wrong field type for Profile ID!!!")

        pt_cnt = 0 # point counter
        arcpy.MakeFeatureLayer_management(outProfPt, "pt_lyr", where_clause)
        with arcpy.da.SearchCursor("pt_lyr", fieldList) as pt_cursor:
            for pt_row in pt_cursor:
                out_f.write(str(pt_cnt) + "," + ",".join(map(str,[val for val in pt_row])) + "\n")
                pt_cnt += 1
        out_f.close()
        del pt_cursor
del line_cursor

# delete layer
arcpy.Delete_management("pt_lyr")

arcpy.AddMessage ("Cleaning workspace")
if arcpy.Exists(workGDB):
    arcpy.Delete_management(workGDB)
