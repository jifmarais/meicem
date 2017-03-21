-- CADFEKO v2017-293043 (x64)
app = cf.GetApplication()
project = app.Project

-- New project
project = app:NewProject()

-- Created geometry: rectangle "Rectangle1"
properties = cf.Rectangle.GetDefaultProperties()
properties.DefinitionMethod = cf.Enums.RectangleDefinitionMethodEnum.BaseAtCentre
properties.Depth = "4"
properties.Label = "Rectangle1"
properties.Width = "4"
Rectangle1 = project.Geometry:AddRectangle(properties)

-- Add scale transform
properties = cf.Scale.GetDefaultProperties()
properties.ScaleFactor = "0.2"
Rectangle1 = project.Geometry["Rectangle1"]
Scale1 = Rectangle1.Transforms:AddScale(properties)

-- Set the frequency to single frequency.
StandardConfiguration1 = project.SolutionConfigurations["StandardConfiguration1"]
FrequencyRange1 = StandardConfiguration1.Frequency
properties = FrequencyRange1:GetProperties()
properties.Start = "1e8"
FrequencyRange1:SetProperties(properties)

-- Updating mesh parameters
MeshSettings = project.Mesher.Settings
properties = MeshSettings:GetProperties()
properties.Advanced.MinElementSize = 37.0511713132585
properties.Advanced.RefinementFactor = 62.4196350581785
MeshSettings:SetProperties(properties)

-- Mesh the model
project.Mesher:Mesh()

-- Created solution entity: PlaneWaveSource1
properties = cf.PlaneWave.GetDefaultProperties()
properties.Label = "PlaneWaveSource1"
PlaneWaveSource1 = project.SolutionConfigurations["StandardConfiguration1"].Sources:AddPlaneWave(properties)

-- Created solution entity: Currents1
properties = cf.Currents.GetDefaultProperties()
properties.ExportSettings.ASCIIEnabled = true
properties.Label = "Currents1"
Currents1 = project.SolutionConfigurations["StandardConfiguration1"].Currents:Add(properties)

-- Created solution configuration: StandardConfiguration2
StandardConfiguration2 = project.SolutionConfigurations:AddStandardConfiguration()

-- Per configuration settings for excitations changed.
project.SolutionConfigurations:SetSourcesPerConfiguration()

-- Modified solution entity: PlaneWaveSource1
PlaneWaveSource1_1 = StandardConfiguration2.Sources["PlaneWaveSource1"]
properties = PlaneWaveSource1_1:GetProperties()
properties.StartPhi = "20"
properties.StartTheta = "30"
PlaneWaveSource1_1:SetProperties(properties)

-- Created solution entity: Currents1
properties = cf.Currents.GetDefaultProperties()
properties.ExportSettings.ASCIIEnabled = true
properties.Label = "Currents1"
Currents1 = project.SolutionConfigurations["StandardConfiguration2"].Currents:Add(properties)

-- Save project
app:SaveAs([[test3_plate.cfx]])

--- Mesh export ---
project.Exporter.Mesh.ExportFileFormat = cf.Enums.ExportFileFormatEnum.NASTRAN

project.Exporter.Mesh.ExportMeshType = cf.Enums.ExportMeshTypeEnum.SimulationMesh

project.Exporter.Mesh.ExportOnlyBoundingFacesEnabled = false

geometryTargets = { project.Geometry["Rectangle1"] }
project.Exporter.Mesh:ExportParts("test3_plate.nas",geometryTargets,{})

-- Add rotate transform
properties = cf.Rotate.GetDefaultProperties()
properties.Angle = "30"
properties.Axis.N = "0"
properties.Axis.V = "1"
Union1 = project.Geometry["Rectangle1"]
Rotate1 = Union1.Transforms:AddRotate(properties)

-- Add rotate transform
properties = cf.Rotate.GetDefaultProperties()
properties.Angle = "20"
Rotate2 = Union1.Transforms:AddRotate(properties)

-- Mesh the model
project.Mesher:Mesh()

Union1 = project.Geometry["Rectangle1"]
geometryTargets = { Union1 }
project.Exporter.Mesh:ExportParts("test3_plate_rotate.nas",geometryTargets,{})

-- Save project
app:SaveAs([[test3_plate_rotate.cfx]])

Rotate2:Delete()
Rotate1:Delete()
