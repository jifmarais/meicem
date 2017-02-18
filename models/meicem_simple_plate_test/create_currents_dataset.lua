-- Create a currents dataset that matches the triangles.

local function myToNumber(val)
    if (val == "NaN" or val== "-NaN" or val == "nan" or val == "-nan") then
        val = 0
    end
    return tonumber(val)
end

local app = pf.GetApplication()
printlist(pf.SurfaceCurrentsAndCharges.GetNames())
local surfaceCurrents = app.Models[1].Configurations[1].SurfaceCurrents
-- local ds = surfaceCurrents[1]:GetDataSet()

-- print(ds)

frequency = 1
triangleIndex = 1
nodeIndex = 1

-- This is not general - assumes a single face (just add a FOR loop)
numTriangles = app.Models[1].Configurations[1].Mesh.TriangleFaces[1].Triangles.Count
print(numTriangles)
local indeces = {}
for ii = 1, numTriangles do
    table.insert(indeces, ii)
    table.insert(indeces, ii)
    table.insert(indeces, ii)
end

ds = pf.DataSet.New()
ds.Axes:Add( pf.Enums.DataSetAxisEnum.Frequency, "GHz", 0.1, 0.1, 1)
ds.Axes:Add( pf.Enums.DataSetAxisEnum.MeshIndex,"m", indeces)

ds.Quantities:Add( "ElectricX", pf.Enums.DataSetQuantityTypeEnum.Complex, "A/m")
ds.Quantities:Add( "ElectricY", pf.Enums.DataSetQuantityTypeEnum.Complex, "A/m")
ds.Quantities:Add( "ElectricZ", pf.Enums.DataSetQuantityTypeEnum.Complex, "A/m")

print(ds)

--- Read in the OS file
plfile = require("pl.file")
-- fileAsString = plfile.read([[/home/jif/Dropbox/Altair/tmp/meicem_simple_plate_test/simple_plate_test.os]])
fileAsString = plfile.read([[/home/jif/gitRepositories/meicem/bin/test1.os]])
-- print(fileAsString)

stringx = require("pl.stringx")

stringio = require("pl.stringio")
f = stringio.open(fileAsString)
-- l1 = f:read()  -- read first line
-- n,m = f:read ('*n','*n') -- read two numbers
ii = 1
processLine = true;
for line in f:lines() do
    local processLine = true
    line = stringx.strip(line)
    
--     if stringx.lfind(line, "##") == 1 
--     then
--         break
--     end
    
    if stringx.startswith(line, "##") 
    then
        processLine = false
    end

    if stringx.startswith(line, "**") 
    then
        processLine = false
    end

    if stringx.startswith(line, "#") 
    then
        processLine = false
    end

    if processLine
    then
        local fields = stringx.split(line)
--         print(#fields)
        if #fields == 31 
        then
            ds[frequency][(ii-1)*3 + 1].ElectricX = myToNumber(fields[14]) + i*myToNumber(fields[15])
            ds[frequency][(ii-1)*3 + 1].ElectricY = myToNumber(fields[16]) + i*myToNumber(fields[17])
            ds[frequency][(ii-1)*3 + 1].ElectricZ = myToNumber(fields[18]) + i*myToNumber(fields[19])
            
            ds[frequency][(ii-1)*3 + 2].ElectricX = myToNumber(fields[20]) + i*myToNumber(fields[21])
            ds[frequency][(ii-1)*3 + 2].ElectricY = myToNumber(fields[22]) + i*myToNumber(fields[23])
            ds[frequency][(ii-1)*3 + 2].ElectricZ = myToNumber(fields[24]) + i*myToNumber(fields[25])
            
            ds[frequency][(ii-1)*3 + 3].ElectricX = myToNumber(fields[26]) + i*myToNumber(fields[27])
            ds[frequency][(ii-1)*3 + 3].ElectricY = myToNumber(fields[28]) + i*myToNumber(fields[29])
            ds[frequency][(ii-1)*3 + 3].ElectricZ = myToNumber(fields[30]) + i*myToNumber(fields[31])
            
            ii = ii+1
        end
    end
end -- iterate over all lines

--- Store new DataSet
storedCurrents = ds:StoreData(pf.Enums.StoredDataTypeEnum.SurfaceCurrentsAndCharges)
