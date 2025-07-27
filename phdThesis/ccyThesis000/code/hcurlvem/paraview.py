import paraview.simple as pv

def save_screenshot(file, datatype='CELLS'):
    fname = file.split('.')[0]
    reader = pv.OpenDataFile(file)

    # 获取渲染视图
    render_view = pv.GetRenderView()

    # 创建一个表示单元数据的渲染器
    dis = pv.Show(reader, render_view, "UnstructuredGridRepresentation")

    readerRep = pv.GetRepresentation()

    # 设置相机位置
    camera = render_view.GetActiveCamera()
    camera.LightFollowCamera = 0
    camera.AmbientColor = [1, 1, 1]
    camera.LightHeadlight = 0
    camera.LightSwitch = 0

    pv.ColorBy(readerRep, (datatype, "cval"))
    pv.Render()
    screenshot_path = fname+"_cuh.png"
    pv.SaveScreenshot(screenshot_path, render_view)

    pv.ColorBy(readerRep, (datatype, "val", "Magnitude"))
    pv.Render()
    screenshot_path = fname+"_uh.png"
    pv.SaveScreenshot(screenshot_path, render_view)

    pv.ColorBy(readerRep, (datatype, "val", "X"))
    pv.Render()
    screenshot_path = fname+"_uhx.png"
    pv.SaveScreenshot(screenshot_path, render_view)

    pv.ColorBy(readerRep, (datatype, "val", "Y"))
    pv.Render()
    screenshot_path = fname+"_uhy.png"
    pv.SaveScreenshot(screenshot_path, render_view)

ff = ["double_band_vem", "five_band_vem", "double_circle_vem", "poly_vem"]
ff = ["two_band_fem", "five_band_fem", "two_circle_fem", "poly_fem"]
for fff in ff:
    save_screenshot(fff + '.vtu', 'POINTS')
