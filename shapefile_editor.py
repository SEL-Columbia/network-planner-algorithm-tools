def addNewColumnsToShapeFile(shapeFile,nearestFacDist, nearestFacID):

    source = ogr.Open(shapeFile, 1)
    layer = source.GetLayer()

    newField1 = ogr.FieldDefn("NearFacDist", ogr.OFTReal)
    newField1.SetWidth(14) 
    newField1.SetPrecision(6)
    layer.CreateField(newField1)

    newField2=ogr.FieldDefn("NearFacID", ogr.OFTInteger)
    layer.CreateField(newField2)

    source = None

    ds = ogr.Open(shapeFile)
    layer = ds.GetLayer(0)


    feat = layer.GetNextFeature()
    while feat is not None:
        feat.GetField("NearFacID")
        feat.SetField("NearFacID", 100)
        layer.SetFeature(feat)
        feat = layer.GetNextFeature()
        
    ds.Destroy()
    # Close the Shapefile
    source = None
 




def writeFieldToShp(shapefile, nearestFacDistDict):
    'Writes a field (provided as a dictionary by FID) to a shapefile.'
    #TODO: fix to allow creation of a new field in a shapefile that already has
    #features... if OGR will allow it. Can copy into a new shapefle as in the
    #soundg examples.
    field="Near"
    ds = ogr.Open(shapefile)
    layer = ds.GetLayer(0)
    feat = layer.GetNextFeature()
    fieldIndex = feat.GetFieldIndex(field)
    if fieldIndex == -1:
        fieldType = getFieldType(nearestFacDistDict.values()[0])
        fieldDefn = ogr.FieldDefn(field, fieldType)
        layer.CreateField(fieldDefn)
    while feat is not None:
        FID = feat.GetFID()
        try:
            fieldValue = nearestFacDistDict[FID]
        except KeyError:
            pass
        else:
            feat.SetField(fieldIndex, fieldValue)
            feat = layer.GetNextFeature()
    #ds=None
    ds.Destroy()
    
    return 0




def getFieldType(fieldValue):
    'Returns OGR field type appropriate for the given value'
    if type(fieldValue) == float:
        return ogr.OFTReal
    elif type(fieldValue) == int:
        return ogr.OFTInteger
    elif type(fieldValue) == str:
        return ogr.OFTString