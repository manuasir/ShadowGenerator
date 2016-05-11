import arcpy
import os
import math
import multiprocessing
import time

############################    Configuration:    ##############################
# Specify scratch workspace
scratchws = r"c:\.." # MUST be a folder, not a geodatabase!

# Specify output field name for the original FID
origfidfield = "ORIG_FID"

# Specify the number of processors (CPU cores) to use (0 to use all available)
cores = 1

# Specify per-process feature count limit, tune for optimal
# performance/memory utilization (0 for input row count divided by cores)
procfeaturelimit = 0

# TIP: Set 'cores' to 1 and 'procfeaturelimit' to 0 to avoid partitioning and
# multiprocessing completely
################################################################################

def message(msg, severity=0):
    print msg

    try:
        for string in msg.split('\n'):
            if severity == 0:
                arcpy.AddMessage(string)
            elif severity == 1:
                arcpy.AddWarning(string)
            elif severity == 2:
                arcpy.AddError(string)
    except:
        pass

def getOidRanges(inputFC, oidfield, count):
    oidranges = []
    if procfeaturelimit > 0:
        rows = arcpy.SearchCursor(inputFC, "", "", oidfield, "%s A" % oidfield)
        minoid = -1
        maxoid = -1
        for r, row in enumerate(rows):
            interval = r % procfeaturelimit
            if minoid < 0 and (interval == 0 or r == count - 1):
                minoid = row.getValue(oidfield)
            if maxoid < 0 and (interval == procfeaturelimit - 1 or r == count - 1):
                maxoid = row.getValue(oidfield)
            if minoid >= 0 and maxoid >= 0:
                oidranges.append([minoid, maxoid])
                minoid = -1
                maxoid = -1
        del row, rows
    return oidranges

def computeShadows(inputFC, outputFC, oidfield, shapefield, heightfield, azimuth, altitude, outputSR="", whereclause=""):
    # Set outputs to be overwritten just in case; each subprocess gets its own environment settings
    arcpy.env.overwriteOutput=True

    # Create in-memory feature class for holding the shadow polygons
    tempshadows = r"c:/users/natalia/desktop/temp/tempshadows.shp"
    arcpy.CreateFeatureclass_management("c:/../temp","tempshadows","POLYGON", "", "", "",outputSR)
    arcpy.AddField_management(tempshadows, origfidfield, "LONG")

    # Open a cursor on the input feature class
    searchfields = ",".join([heightfield, oidfield, shapefield])
    rows = arcpy.SearchCursor(inputFC, whereclause, "", searchfields)

    # Open an insert cursor on the in-memory feature class
    tempinscur = arcpy.InsertCursor(tempshadows)

    # Create array for holding shadow polygon vertices
    arr = arcpy.Array()

    # Compute the shadow offsets.
    spread = 1/math.tan(altitude)

    # Read the input features
    for row in rows:
        oid = int(row.getValue(oidfield))
        shape = row.getValue(shapefield)
        height = float(row.getValue(heightfield))

        # Compute the shadow offsets.
        x = -height * spread * math.sin(azimuth)
        y = -height * spread * math.cos(azimuth)

        # Copy the original shape as a new feature
        tempnewrow = tempinscur.newRow()
        tempnewrow.setValue(origfidfield, oid) # Copy the original FID value to the new feature
        tempnewrow.shape = shape
        tempinscur.insertRow(tempnewrow)

        # Compute the wall shadow polygons and insert them into the in-memory feature class
        for part in shape:
            for i,j in enumerate(range(1,part.count)):
                pnt0 = part[i]
                pnt1 = part[j]
                if pnt0 is None or pnt1 is None:
                    continue # skip null points so that inner wall shadows can also be computed

                # Compute the shadow offset points
                pnt0offset = arcpy.Point(pnt0.X+x,pnt0.Y+y)
                pnt1offset = arcpy.Point(pnt1.X+x,pnt1.Y+y)

                # Construct the shadow polygon and insert it to the in-memory feature class
                [arr.add(pnt) for pnt in [pnt0,pnt1,pnt1offset,pnt0offset,pnt0]]
                tempnewrow.shape = arr
                tempnewrow.setValue(origfidfield, oid) # Copy the original FID value to the new feature
                tempinscur.insertRow(tempnewrow)
                arr.removeAll() # Clear the array so it can be reused

    # Clean up the insert cursor
    del tempnewrow, tempinscur

    # Dissolve the shadow polygons to the output feature class
    dissolved = arcpy.Dissolve_management(tempshadows, outputFC, origfidfield).getOutput(0)

    # Clean up the in-memory workspace
    arcpy.Delete_management("in_memory")

    return dissolved

if __name__ == "__main__":
    arcpy.env.overwriteOutput=True

    # Read in parameters
    inputFC = r"c:\..\buildings.shp"
    outputFC = r"c:\..\sombras.shp"
    heightfield = "alt" #Must be in the same units as the coordinate system!
    azimuth = math.radians(float(200)) #Must be in degrees
    altitude = math.radians(float(35)) #Must be in degrees

    # Get field names
    desc = arcpy.Describe(inputFC)
    shapefield = desc.shapeFieldName
    oidfield = desc.oidFieldName

    count = int(arcpy.GetCount_management(inputFC).getOutput(0))
    message("Figuras a procesar: %d" % count)

    #Export output spatial reference to string so it can be pickled by multiprocessing
    if arcpy.env.outputCoordinateSystem:
        outputSR = arcpy.env.outputCoordinateSystem.exportToString()
    elif desc.spatialReference:
        outputSR = desc.spatialReference.exportToString()
    else:
        outputSR = ""

    # Configure partitioning
    if cores == 0:
        cores = multiprocessing.cpu_count()
    if cores > 1 and procfeaturelimit == 0:
        procfeaturelimit = int(math.ceil(float(count)/float(cores)))

     # Start timing
    start = time.clock()

    # Partition row ID ranges by the per-process feature limit
    oidranges = getOidRanges(inputFC, oidfield, count)

    if len(oidranges) > 0: # Use multiprocessing
        message("Calculando sombras. usando multiprocesamiento (%d procesos, %d trabajos de %d figuras maximas cada uno) ..." % (cores, len(oidranges), procfeaturelimit))

        # Create a Pool of subprocesses
        pool = multiprocessing.Pool(cores)
        jobs = []

        # Get the appropriately delmited field name for the OID field
        oidfielddelimited = arcpy.AddFieldDelimiters(inputFC, oidfield)

        # Ensure the scratch workspace folder exists
        if not os.path.exists(scratchws):
            os.mkdir(scratchws)

        for o, oidrange in enumerate(oidranges):
            # Build path to temporary output feature class (dissolved shadow polygons)
            # Named e.g. <scratchws>\dissolvedshadows0000.shp
            tmpoutput = os.path.join(scratchws, "%s%04d.shp" % ("dissolvedshadows", o))

            # Build a where clause for the given OID range
            whereclause = "%s >= %d AND %s <= %d" % (oidfielddelimited, oidrange[0], oidfielddelimited, oidrange[1])

            # Add the job to the multiprocessing pool asynchronously
            jobs.append(pool.apply_async(computeShadows, (inputFC, tmpoutput, oidfield, shapefield, heightfield, azimuth, altitude, outputSR, whereclause)))

        # Clean up worker pool; waits for all jobs to finish
        pool.close()
        pool.join()

         # Get the resulting outputs (paths to successfully computed dissolved shadow polygons)
        results = [job.get() for job in jobs]

        try:
            # Merge the temporary outputs
            message("Uniendo capas temporales %s ..." % outputFC)
            arcpy.Merge_management(results, outputFC)
        finally:
            # Clean up temporary data
            message("Eliminando datos temporales ...")
            for result in results:
                message("Borrando %s" % result)
                try:
                    arcpy.Delete_management(result)
                except:
                    pass
    else: # Use a single process
        message("Calculando sombras ...")
        computeShadows(inputFC, outputFC, oidfield, shapefield, heightfield, azimuth, altitude, outputSR)

    # Stop timing and report duration
    end = time.clock()
    duration = end - start
    hours, remainder = divmod(duration, 3600)
    minutes, seconds = divmod(remainder, 60)
    message("Completado en %d:%d:%f" % (hours, minutes, seconds))