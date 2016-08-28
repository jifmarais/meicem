-- Create a currents dataset that matches the triangles.

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
--             print(fields[1], 
--     --               fields[2], 
--     --               fields[3], 
--     --               fields[4], 
--     --               fields[5], 
--     --               fields[6], 
--     --               fields[7], 
--     --               fields[8], 
--     --               fields[9], 
--     --               fields[10], 
--     --               fields[11], 
--     --               fields[11],
--     --               fields[12], 
--     --               fields[13], 
--                   fields[14], -- Re (Jx C1)
--                   fields[15], -- Im (Jx C1)
--                   fields[16], -- Re (Jy C1)
--                   fields[17], -- Im (Jy C1) 
--                   fields[18], -- Re (Jz C1)
--                   fields[19], -- Im (Jz C1)
--                   fields[20], -- Re (Jx C2)
--                   fields[21], -- Im (Jx C2)
--                   fields[22], -- Re (Jy C2)
--                   fields[23], -- Im (Jy C2)
--                   fields[24], -- Re (Jz C2)
--                   fields[25], -- Im (Jz C2)
--                   fields[26], -- Re (Jx C3)
--                   fields[27], -- Im (Jx C3)
--                   fields[28], -- Re (Jy C3)
--                   fields[29], -- Im (Jy C3)
--                   fields[30], -- Re (Jz C3)
--                   fields[31]) -- Im (Jz C3)
            ds[frequency][(ii-1)*3 + 1].ElectricX = tonumber(fields[14]) + i*tonumber(fields[15])
            ds[frequency][(ii-1)*3 + 1].ElectricY = tonumber(fields[16]) + i*tonumber(fields[17])
            ds[frequency][(ii-1)*3 + 1].ElectricZ = tonumber(fields[18]) + i*tonumber(fields[19])
            
            ds[frequency][(ii-1)*3 + 2].ElectricX = tonumber(fields[20]) + i*tonumber(fields[21])
            ds[frequency][(ii-1)*3 + 2].ElectricY = tonumber(fields[22]) + i*tonumber(fields[23])
            ds[frequency][(ii-1)*3 + 2].ElectricZ = tonumber(fields[24]) + i*tonumber(fields[25])
            
            ds[frequency][(ii-1)*3 + 3].ElectricX = tonumber(fields[26]) + i*tonumber(fields[27])
            ds[frequency][(ii-1)*3 + 3].ElectricY = tonumber(fields[28]) + i*tonumber(fields[29])
            ds[frequency][(ii-1)*3 + 3].ElectricZ = tonumber(fields[30]) + i*tonumber(fields[31])
            
            ii = ii+1
        end
    end
end -- iterate over all lines

--- Store new DataSet
storedCurrents = ds:StoreData(pf.Enums.StoredDataTypeEnum.SurfaceCurrentsAndCharges)
