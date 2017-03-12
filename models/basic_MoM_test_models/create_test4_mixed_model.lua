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

-- Created geometry: rectangle "Rectangle2"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "3"
properties.Label = "Rectangle2"
properties.Origin.N = "0"
properties.Origin.U = "-1.6"
properties.Origin.V = "-1.6"
properties.Width = "0.2"
Rectangle2 = project.Geometry:AddRectangle(properties)

-- Created geometry: rectangle "Rectangle3"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "0.4"
properties.Label = "Rectangle3"
properties.Origin.N = "0"
properties.Origin.U = "-1.4"
properties.Origin.V = "-1.6"
properties.Width = "1.2"
Rectangle3 = project.Geometry:AddRectangle(properties)

-- Created geometry: ellipse "Ellipse1"
properties = cf.Ellipse.GetDefaultProperties()
properties.Centre.N = "0"
properties.Centre.U = "0.4"
properties.Centre.V = "1"
properties.Label = "Ellipse1"
properties.RadiusU = "0.4"
properties.RadiusV = "0.5656854249"
Ellipse1 = project.Geometry:AddEllipse(properties)

-- Created geometry: rectangle "Rectangle4"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "0.8"
properties.Label = "Rectangle4"
properties.Origin.N = "0"
properties.Origin.U = "0.6"
properties.Origin.V = "-1.6"
properties.Width = "1"
Rectangle4 = project.Geometry:AddRectangle(properties)

-- Created geometry: union "Union1"
targets = { Rectangle3, Ellipse1, Rectangle2, Rectangle4 }
project.Geometry:Union(targets)

-- Created geometry: subtract "Subtract1"
Union1 = project.Geometry["Union1"]
targets = { Union1 }
Rectangle1 = project.Geometry["Rectangle1"]
project.Geometry:Subtract(Rectangle1, targets)

-- Created geometry: sphere "Sphere1"
properties = cf.Spheroid.GetDefaultProperties()
properties.Centre.N = "0.3"
properties.Centre.U = "-3.4"
properties.Centre.V = "0"
properties.Label = "Sphere1"
properties.Radius = "1"
Sphere1 = project.Geometry:AddSpheroid(properties)

-- Created geometry: split "Split_back1"
properties = cf.Split.GetDefaultProperties()
properties.Origin.N = "0"
properties.Origin.U = "-4.2"
properties.Origin.V = "-0.3"
properties.Plane = cf.Enums.SplitPlanesEnum.VN
Sphere1 = project.Geometry["Sphere1"]
project.Geometry:Split(Sphere1, properties)

-- Deleted geometry: split "Split_back1"
Split_back1 = project.Geometry["Split_back1"]
Split_back1:Delete()

-- Deleting geometry entities
Split_front1 = project.Geometry["Split_front1"]
Face8 = Split_front1.Faces["Face8"]
Face8:Delete()

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
app:SaveAs([[test4_mixed.cfx]])

--- Mesh export ---
project.Exporter.Mesh.ExportFileFormat = cf.Enums.ExportFileFormatEnum.NASTRAN

project.Exporter.Mesh.ExportMeshType = cf.Enums.ExportMeshTypeEnum.SimulationMesh

project.Exporter.Mesh.ExportOnlyBoundingFacesEnabled = false

project.Exporter.Mesh:ExportParts("test4_mixed.nas",project.Geometry:Items(),{})
