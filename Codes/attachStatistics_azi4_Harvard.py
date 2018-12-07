# Attaches statistics from CSV files
# Section name field added and populated with name: Section"SECTION ID"
# Extracts a point with the max slope along each section line
# CSV files in InputStatistics folder
# CSV file name format: Slip_Statistics_"SECTION FILE NAME"

# Import system modules
import arcpy
import os, sys
import operator
from arcpy.sa import *
import math

# licences
arcpy.CheckOutExtension('Spatial')

#execfile('ArcGIScodes/attachStatistics_azi4_Harvard.py')
#-------------------------------------------------------------------------------------

# Workspaces
#workDir = r'.'
workDir = r'C:\Users\frw313\Desktop\GNS'
workGDB = os.path.join(workDir, "DrumMtn.gdb")
scratchGDB =  os.path.join(workDir, "scratch.gdb")

# Set inputs
inDEM = os.path.join(workGDB, "sevier_comb")

inSection = raw_input("Enter name of section GIS shapefile: ")
secIDField = "SectionID"
inCSV = raw_input("Enter name of section csv - between Statistics_ and .csv: ")
##inCSV = "profiles1"

inStatCSV = os.path.join(workDir, "OutputTables", "Slip_Statistics_" + inCSV + ".csv")
inSlipField = "SlipRate_avg_horiz"

#------------------------------------------------------------------------------------------
# Outputs
outSection = inCSV + "_Statistics"
outLabelPt = inCSV + "_StatsLabel"
secNameField = "SectionName"

#Outputs
slipWEField = inSlipField + "_WE"
slipSNField = inSlipField + "_SN"
aziField = 'Azimuth'
arrowField = "ArrowDir"

# set dictionary of field types
fieldName = {1:"SlipMean", 2:"SlipMedian", 3:"SlipMode",
             4: "Slip_CIlow", 5: "Slip_CIhigh",
             6: "SlipRate_mean", 7: "SlipRate_CIlow", 8: "SlipRate_CIhigh", 9: "SlipRate_avg_horiz",
             10: "SlipRate_CILow_horiz", 11: "SlipRate_CIhigh_horiz", 12: "SlipRate_avg_vert",
             13: "SlipRate_CIlow_vert", 14: "SlipRate_CIhigh_vert", 15:"DownToEast", 16: "AgeMin", 17:"AgeMax",
             18: "Azimuth", 19: "SlipRate_avg_horiz_WE",
             20:"SlipRate_avg_horiz_SN", 21:"ArrowDir"}

# create scratch GDB if it doesn't exist
if not arcpy.Exists(scratchGDB):
    print "Creating scratch GDB"
    arcpy.CreateFileGDB_management(workDir, os.path.basename(scratchGDB))

# Set enviroments
arcpy.env.workspace = workGDB
arcpy.env.scratchWorkspace = scratchGDB
arcpy.env.overwriteOutput = True

# Copy section lines and add necessary fields
print "Creating output section feature class"
print " Copying features"
arcpy.CopyFeatures_management(inSection, outSection)
print " Adding ID"
arcpy.AddField_management(outSection, secNameField, "TEXT")
print " Adding other fields"
for f_cnt in range(1,18):
    arcpy.AddField_management(outSection, fieldName[f_cnt], "DOUBLE")
print "Calculating section name"
#arcpy.CalculateField_management(outSection, secNameField, '"proeedSetion_" + str(!' + secIDField + '!)', "PYTHON")
arcpy.CalculateField_management(outSection, secNameField, '"proeed_Setion_" + str(!' + secIDField + '!)', "PYTHON")

# Creating an empty point feature class to hold max slope points"
print "Creating an empty point feature class to hold max slope points"
if arcpy.Exists(outLabelPt):
    arcpy.Delete_management(outLabelPt)
arcpy.CreateFeatureclass_management(workGDB, outLabelPt, 'POINT', '#', '#', '#', inSection)
arcpy.AddField_management(outLabelPt, secNameField, "TEXT")

# Read statistics from file and extract max slope
line_num = 0
lines = open(inStatCSV, "r").readlines() # list of lines in the CSV file
num_sec = len(lines[0].strip("\n").split(",")) # get number of section
for i in range(0, num_sec): # read section by section  (column by column)

    SectionName = lines[0].strip("\n").split(",")[i]   # read section name from file

    # Attaching statistics to section lines
    print SectionName
    where_clause = secNameField + " = '" + SectionName + "'"
    arcpy.MakeFeatureLayer_management(outSection, "sec", where_clause) # select section
    # check if section exists
    if int(arcpy.GetCount_management("sec").getOutput(0)) > 0:
        print " Attcahing statistics from file"
        for f_cnt in range(1,18):
            arcpy.CalculateField_management("sec", fieldName[f_cnt], lines[f_cnt].strip("\n").split(",")[i])

##        # Extract max slope
##        print " Converting line to raster"
##        arcpy.ClearEnvironment("cellSize")
##        arcpy.FeatureToRaster_conversion("sec", secIDField, "temp_rast", "100")
##        print " Extracting elevation along section line"
##        arcpy.env.cellSize = 20
##        outCon = Con("temp_rast", inDEM)
##        print " Calculating slope"
##        outSlope = Slope(outCon)
##        print " Converting slope to points"
##        arcpy.RasterToPoint_conversion(outSlope, "temp_pt")
##        print " Adding section name"
##        arcpy.AddField_management("temp_pt", secNameField, "TEXT")
##        arcpy.CalculateField_management("temp_pt", secNameField, '"' + SectionName + '"')
##        print " Performing near function"
##        arcpy.Near_analysis("temp_pt", outSection, 5)
##        print " Selecting points close to the line and extracting the max slope"
##        arcpy.MakeFeatureLayer_management("temp_pt", "temp_pt2", '"NEAR_FID" <> -1')
##        #Build a dictionary of objectids and values:
##        all_entries = {key:value for (key, value) in arcpy.da.SearchCursor("temp_pt2",['OID@','grid_code'])} #Change VALUECOLUMN to the name of your column
##        #Find max value and extract objectid
##        id_max = max(all_entries.iteritems(), key=operator.itemgetter(1))[0]
##        print " Adding point with max slope to the output point feature class"
##        arcpy.MakeFeatureLayer_management("temp_pt2", "max_pt", '"OBJECTID" = ' + str(id_max))
##        arcpy.Append_management("max_pt", outLabelPt, "NO_TEST")
##    else:
##        print " " + SectionName + " doesn't exist"

print "Removing empty records from section lines"
arcpy.MakeFeatureLayer_management(outSection, "sec_del", fieldName[1] + ' IS NULL')
if int(arcpy.GetCount_management("sec_del").getOutput(0)) > 0:
    arcpy.DeleteRows_management("sec_del")

# Adding fields for azimuth and vector components of slip
print "Adding fields"
arcpy.AddField_management(outSection, aziField, "DOUBLE")
arcpy.AddField_management(outSection, slipWEField, "DOUBLE")
arcpy.AddField_management(outSection, slipSNField, "DOUBLE")
arcpy.AddField_management(outSection, arrowField, "DOUBLE")

# Calculating azimuth and rates in vector components of slip"
print "Calculating.."
azi_calc = "180  + math.atan2((!Shape.firstpoint.x! - !Shape.lastpoint.x! ),(!Shape.firstpoint.y! - !Shape.lastpoint.y! ) ) * (180 / math.pi)"
arcpy.CalculateField_management(outSection, aziField, azi_calc, "PYTHON_9.3")
we_calc = "abs(" "!" + inSlipField + "!  * math.sin( (!" + aziField + "!*math.pi)/180))"
arcpy.CalculateField_management(outSection, slipWEField, we_calc, "PYTHON")
sn_calc = "abs(" "!" + inSlipField + "!  * math.cos( (!" + aziField + "!*math.pi)/180))"
arcpy.CalculateField_management(outSection, slipSNField, sn_calc,  "PYTHON")

# Calculating arrow direction
print "Arrow Direction"
arcpy.MakeFeatureLayer_management(outSection, "arrow")
where_clause = fieldName[15] + " = 1"
arcpy.SelectLayerByAttribute_management("arrow", "NEW_SELECTION", where_clause)
arcpy.CalculateField_management("arrow", arrowField, "!" + aziField + "!", "PYTHON")
arcpy.SelectLayerByAttribute_management("arrow", "SWITCH_SELECTION")
arcpy.CalculateField_management("arrow", arrowField, "!" + aziField + "! + 180", "PYTHON")

print "Attaching statistics to points"
fieldList = [fieldName[i] for i in range(1,22)]
arcpy.JoinField_management(outLabelPt, secNameField, outSection, secNameField, fieldList)

print "Cleaning workspace"
arcpy.Delete_management("temp_rast")
arcpy.Delete_management("temp_pt")

sys.exit()

