import paraview.simple as pv

# 加载数据文件
filename = "./five_layer_fem.vtu"
reader = pv.XMLUnstructuredGridReader(FileName=[filename])

# 创建视图并显示数据
view = pv.CreateRenderView()
view.Background = [1, 1, 1]  # 设置背景颜色为白色
view.ViewSize = [800, 800]   # 设置视图大小
view.DisplayAxes = 1         # 显示坐标轴

# 使用 Jet 颜色预设进行可视化
coloring = pv.ColorBy(reader.PointData.scalars)
lut = pv.GetColorTransferFunction("yourcolorpreset")  # yourcolorpreset替换为你想要使用的颜色预设名称
lut.RGBPoints = [0, 0, 0, 1, 1, 1, 0, 0]
coloring.LookupTable = lut

# 创建渲染器并将数据添加到渲染器中
render = pv.CreateRepresentation()
render.Representation = "Surface With Edges"
render.ColorArrayName = ["POINTS", "yourcolorarray"]  # yourcolorarray替换为你想要使用的颜色数据的名称
render.Visibility = 1

# 渲染并保存结果
pv.Render()
pv.SaveScreenshot("output.png", view, ImageResolution=[800, 800])

