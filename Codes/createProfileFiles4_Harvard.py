# Creates points alon section lines at equal intervals, extracts elevation and geology age and exports data to separate files for each section
# Required: Section line feature class (with section ID field), DEM, geology layer (with age min and age max fields), point distance
# Outputs: Point feature class and CSV files for each section

# Import system modules
import arcpy
import os, sys
from arcpy.sa import *

# licences
arcpy.CheckOutExtension('Spatial')

#execfile('ArcGIScodes/createProfileFiles4_Harvard.py')

#-------------------------------------------------------------------------------------

# Workspaces
#workDir = r'.'
workDir = r'C:\Users\frw313\Desktop\GNS'
workGDB = os.path.join(workDir, "DrumMtn.gdb")

# Set inputs
inDEM = os.path.join(workGDB, "sevier_comb")
#inDEM = r'DrumMtn.gdb\sevier_comb'
#inGeol = r'DrumMtn.gdb\geologicunits_reproj_gd'
inGeol = os.path.join(workGDB, "geologicunits_reproj_gd")
#C:\Users\frw313\Desktop\GNS\Reproj\Reproj-20181114T170027Z-001\Reproj\geologicunits_reproj.shp
ageMinField = "AGE_Min"
ageMaxField = "AGE_Max"

inSection = raw_input("Enter name of section shapefile: ")
secIDField = "SectionID"

# point distance
pointDist = 1

#------------------------------------------------------------------------------------------
# Outputs
outFolder = os.path.join(workDir, "Outputs_" + inSection)
outSecPoints = inSection + "_points"
distField = "Distance"

# set dictionary of field types
fieldTypeDict = {"Integer":"LONG", "SmallInteger":"SHORT", "String":"TEXT"}

# Define function for createing points along line
def points_along_line(line_lyr, pnt_layer, sec_id, pnt_dist, dist_f):
    """
    line_lyr (feature layer) - Single part line
    pnt_layer (feature layer) - Path to point feature class
    sec_id (attribute field) - Name of the id field (to be  added to points)
    pnt_dist (integer) - Interval distance in map units between points
    dist_f (attribute field) - Field to store distance along line
    """

    search_cursor = arcpy.da.SearchCursor(line_lyr, ['SHAPE@', sec_id])
    insert_cursor = arcpy.da.InsertCursor(pnt_layer, ['SHAPE@', sec_id, dist_f])

    for row in search_cursor:
        distance = 0
        for dist in range(0, int(row[0].length), pnt_dist):
            point = row[0].positionAlongLine(dist).firstPoint
            point_id = row[1]
            insert_cursor.insertRow([point, point_id, distance])
            distance = distance + pointDist

# Set enviroments
arcpy.env.workspace = workGDB
arcpy.env.overwriteOutput = True

# Create output folder if it doesnt exist
if not os.path.exists(outFolder):
    print "Creating output directory"
    os.makedirs(outFolder)

# Get the field type for SectionID field
for f_pair in [(f.name, f.type) for f in arcpy.ListFields(inSection)]:
    if f_pair[0] == secIDField:
        secIDFieldType = fieldTypeDict[f_pair[1]]

# Create point layer
print "Creating point layer"
if arcpy.Exists(outSecPoints):
    arcpy.Delete_management(outSecPoints)
arcpy.CreateFeatureclass_management(workGDB, outSecPoints, 'POINT', '#', '#', '#', inSection)
arcpy.AddField_management(outSecPoints, secIDField, secIDFieldType)
arcpy.AddField_management(outSecPoints, distField, "DOUBLE")

print "Creating points along section lines"
points_along_line(inSection, outSecPoints, secIDField, pointDist, distField)

# Extracting elevation
print "Extracting elevation"
ExtractValuesToPoints(outSecPoints, inDEM, "temp")
print " Attaching elevation to section points"
arcpy.JoinField_management(outSecPoints, "OBJECTID", "temp", "OBJECTID", "RASTERVALU")
print " Renaming field"
arcpy.AlterField_management(outSecPoints, "RASTERVALU", "Elevation", "#", "#", "#", "#", "True")
print " Deleting temp layer"
arcpy.Delete_management("temp")

# Attaching geology age
print "Attaching geology age"
if ageMinField in [f.name for f in arcpy.ListFields(outSecPoints)]:
    print " Deleting " + ageMinField +" field"
    arcpy.DeleteField_management(outSecPoints, ageMinField)
if ageMaxField in [f.name for f in arcpy.ListFields(outSecPoints)]:
    print " Deleting " + ageMaxField +" field"
    arcpy.DeleteField_management(outSecPoints, ageMaxField)
print " Performing Identity"
arcpy.Identity_analysis(outSecPoints, inGeol, "temp")
print " Joining geology age fields to summary points"
arcpy.JoinField_management(outSecPoints, "OBJECTID", "temp", "FID_" + outSecPoints, [ageMinField, ageMaxField])
print " Deleting temp layer"
arcpy.Delete_management("temp")

# Writing section files
print "Writing files"
with arcpy.da.SearchCursor(inSection, [secIDField]) as line_cursor:
    for line_row in line_cursor:
        secID = line_row[0] # get section id
        outFile = os.path.join(outFolder, "Section_" + str(secID) + ".csv") # set output file
        out_f = open(outFile, "w") # open file for writing
        #out_f.write("{0},{1},{2},{3},{4}\n".format(secIDField, distField, "Elevation", "Age_Min_KY", "Age_Max_KY")) # write header

        if secIDFieldType in ["LONG", "SHORT"]:
            where_clause = '"' + secIDField + '" = ' + str(secID)  # set query for serching points
        elif secIDFieldType in ["TEXT"]:
            where_clause = '"' + secIDField + '" = ' + "'" + str(secID) + "'"  # set query for serching points
        else:
            sys.exit("Wrong field type for Section ID")

        arcpy.MakeFeatureLayer_management(outSecPoints, "pt_lyr", where_clause)
        with arcpy.da.SearchCursor("pt_lyr", [secIDField, distField, "Elevation", ageMinField, ageMaxField]) as pt_cursor:
            for pt_row in pt_cursor:
                age_min = float("{0:.3f}".format(pt_row[3]))
                age_max = float("{0:.3f}".format(pt_row[4]))
                out_f.write("{0},{1},{2},{3:.3f},{4:.3f}\n".format(pt_row[0], pt_row[1], pt_row[2], age_min, age_max))
        out_f.close()
        del pt_cursor
del line_cursor

