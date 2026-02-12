# trace generated using paraview version 6.0.1
#import paraview
#paraview.compatibility.major = 6
#paraview.compatibility.minor = 0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
ResetSession()

import os

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Set paths relative to the script directory
validation_file_path = os.path.join(script_dir, 'example', 'biv_with_outflow_tracts', 'validation_doste_combined.vtu')
png_output_path = os.path.join(script_dir, 'example', 'biv_with_outflow_tracts')

fiber_families = ['f', 's', 'n']
fiber_family_names = {'f': 'fiber', 's': 'sheet', 'n': 'sheet-normal'}

# create a new 'XML Unstructured Grid Reader'
validation_doste_combinedvtu = XMLUnstructuredGridReader(registrationName='validation_doste_combined.vtu', FileName=[validation_file_path])

# Properties modified on validation_doste_combinedvtu
validation_doste_combinedvtu.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# Use a solid white background for all screenshots
renderView1.UseColorPaletteForBackground = 0
renderView1.BackgroundColorMode = "Single Color"
renderView1.Background = [1.0, 1.0, 1.0] 

# Make axis titles/labels visible on white background
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# Make orientation axes (X/Y/Z) labels visible on white background
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]

# show data in view
validation_doste_combinedvtuDisplay = Show(validation_doste_combinedvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
validation_doste_combinedvtuDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
validation_doste_combinedvtuDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')

# Ensure scalar bar text is visible on white background
fLUTColorBar = GetScalarBar(fLUT, renderView1)
fLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
fLUTColorBar.LabelColor = [0.0, 0.0, 0.0]

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')

# get 2D transfer function for 'f'
fTF2D = GetTransferFunction2D('f')


# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=validation_doste_combinedvtu)

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Surface'

# hide data in view
Hide(validation_doste_combinedvtu, renderView1)

# show color bar/color legend
cellDatatoPointData1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=cellDatatoPointData1,
    SeedType='Line')

# toggle interactive widget visibility (only when running from the GUI)
ShowInteractiveWidgets(proxy=streamTracer1.SeedType)

# Properties modified on streamTracer1
streamTracer1.SeedType = 'Point Cloud'

# Properties modified on streamTracer1.SeedType
streamTracer1.SeedType.Radius = 60.0
streamTracer1.SeedType.NumberOfPoints = 200000

# show data in view
streamTracer1Display = Show(streamTracer1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
streamTracer1Display.Representation = 'Surface'

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# show color bar/color legend
streamTracer1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on streamTracer1Display
streamTracer1Display.RenderLinesAsTubes = 1

# Properties modified on streamTracer1Display
streamTracer1Display.LineWidth = 2.0

# set scalar coloring
ColorBy(streamTracer1Display, ('POINTS', 'alpha_combined'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(fLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
streamTracer1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
streamTracer1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'alpha_combined'
alpha_combinedLUT = GetColorTransferFunction('alpha_combined')

# get opacity transfer function/opacity map for 'alpha_combined'
alpha_combinedPWF = GetOpacityTransferFunction('alpha_combined')

# get 2D transfer function for 'alpha_combined'
alpha_combinedTF2D = GetTransferFunction2D('alpha_combined')

# set active source
SetActiveSource(validation_doste_combinedvtu)

# toggle interactive widget visibility (only when running from the GUI)
HideInteractiveWidgets(proxy=streamTracer1.SeedType)

# get color legend/bar for alpha_combinedLUT in view renderView1
alpha_combinedLUTColorBar = GetScalarBar(alpha_combinedLUT, renderView1)

# Ensure scalar bar text is visible on white background
alpha_combinedLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
alpha_combinedLUTColorBar.LabelColor = [0.0, 0.0, 0.0]

# change scalar bar placement
alpha_combinedLUTColorBar.Set(
    WindowLocation='Any Location',
    Position=[0.8227593152064452, 0.22451317296678122],
    ScalarBarLength=0.32999999999999996,
)

# change scalar bar placement
alpha_combinedLUTColorBar.Set(
    Position=[0.8277945619335347, 0.3115693012600229],
    ScalarBarLength=0.3299999999999996,
)

# change scalar bar placement
alpha_combinedLUTColorBar.Set(
    Position=[0.8529707955689829, 0.320067884829428],
    ScalarBarLength=0.32999999999999957,
)

# set active source
SetActiveSource(streamTracer1)

# toggle interactive widget visibility (only when running from the GUI)
ShowInteractiveWidgets(proxy=streamTracer1.SeedType)

# set active source
SetActiveSource(validation_doste_combinedvtu)

# toggle interactive widget visibility (only when running from the GUI)
HideInteractiveWidgets(proxy=streamTracer1.SeedType)

# show data in view
validation_doste_combinedvtuDisplay = Show(validation_doste_combinedvtu, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
validation_doste_combinedvtuDisplay.SetScalarBarVisibility(renderView1, True)

# Properties modified on validation_doste_combinedvtuDisplay
validation_doste_combinedvtuDisplay.Opacity = 0.4

# Properties modified on validation_doste_combinedvtuDisplay
validation_doste_combinedvtuDisplay.Opacity = 0.1

# turn off scalar coloring
ColorBy(validation_doste_combinedvtuDisplay, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(fLUT, renderView1)

# set active source
SetActiveSource(streamTracer1)

# toggle interactive widget visibility (only when running from the GUI)
ShowInteractiveWidgets(proxy=streamTracer1.SeedType)

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(validation_doste_combinedvtu)

# toggle interactive widget visibility (only when running from the GUI)
HideInteractiveWidgets(proxy=streamTracer1.SeedType)

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(993, 706)

# current camera placement for renderView1
renderView1.Set(
    CameraPosition=[52.04115429873035, -45.38613870538731, 216.54244474035633],
    CameraFocalPoint=[-21.838666915893594, -60.56520140171051, 48.49395084381101],
    CameraViewUp=[-0.11136472699136057, 0.9929446544942593, -0.040728499768785655],
    CameraParallelScale=102.19338444140465,
)


# save screenshot
SaveScreenshot(filename=os.path.join(png_output_path, 'doste_fiber.png'), viewOrLayout=renderView1, location=16, ImageResolution=[993, 706], TransparentBackground=0)

# set active source
SetActiveSource(streamTracer1)

# toggle interactive widget visibility (only when running from the GUI)
ShowInteractiveWidgets(proxy=streamTracer1.SeedType)

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=validation_doste_combinedvtu)

# Properties modified on slice1.SliceType
slice1.SliceType.Set(
    Origin=[-25.530508094611232, -64.79965909808317, 42.457392775039196],
    Normal=[0.3179741978250393, 0.5256304018456659, 0.7890532872836198],
)

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Hide streamTracer1 before the loop so it doesn't appear in screenshots
Hide(streamTracer1, renderView1)

# Loop through fiber families to generate screenshots
for fiber_family in fiber_families:
    family_name = fiber_family_names[fiber_family]
    
    # create a new 'Stream Tracer' for this fiber family
    streamTracer_current = StreamTracer(registrationName=f'StreamTracer_{family_name}', Input=cellDatatoPointData1,
        SeedType='Line')
    
    # Properties modified on streamTracer_current
    streamTracer_current.SeedType = 'Point Cloud'
    streamTracer_current.Vectors = ['POINTS', fiber_family]
    
    # Properties modified on streamTracer_current.SeedType
    streamTracer_current.SeedType.Radius = 60.0
    streamTracer_current.SeedType.NumberOfPoints = 200000
    
    # show data in view
    streamTracer_currentDisplay = Show(streamTracer_current, renderView1, 'GeometryRepresentation')
    
    # trace defaults for the display properties.
    streamTracer_currentDisplay.Representation = 'Surface'
    streamTracer_currentDisplay.RenderLinesAsTubes = 1
    streamTracer_currentDisplay.LineWidth = 2.0
    
    # set scalar coloring
    ColorBy(streamTracer_currentDisplay, ('POINTS', 'alpha_combined'))
    
    # rescale color and/or opacity maps used to include current data range
    streamTracer_currentDisplay.RescaleTransferFunctionToDataRange(True, False)
    
    # show color bar/color legend
    streamTracer_currentDisplay.SetScalarBarVisibility(renderView1, True)
    
    # create a new 'Glyph'
    glyph1 = Glyph(registrationName='Glyph1', Input=slice1,
        GlyphType='Arrow')

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=glyph1.GlyphType)

    # Properties modified on glyph1
    glyph1.Set(
        GlyphType='Line',
        OrientationArray=['CELLS', fiber_family],
        ScaleArray=['POINTS', 'No scale array'],
        GlyphMode='Uniform Spatial Distribution (Surface Sampling)',
    )

    # show data in view
    glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    glyph1Display.Representation = 'Surface'

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(glyph1Display, ('POINTS', 'alpha_combined'))

    # rescale color and/or opacity maps used to include current data range
    glyph1Display.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, True)

    # Properties modified on glyph1Display
    glyph1Display.LineWidth = 2.0

    # Properties modified on glyph1Display
    glyph1Display.RenderLinesAsTubes = 1

    # set active source
    SetActiveSource(validation_doste_combinedvtu)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=glyph1.GlyphType)

    # set active source
    SetActiveSource(cellDatatoPointData1)

    # set active source
    SetActiveSource(slice1)

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=slice1.SliceType)

    # turn off scalar coloring
    ColorBy(slice1Display, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(fLUT, renderView1)

    # hide data in view
    Hide(validation_doste_combinedvtu, renderView1)

    # set active source
    SetActiveSource(validation_doste_combinedvtu)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice1.SliceType)

    # set active source
    SetActiveSource(glyph1)

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=glyph1.GlyphType)

    # Properties modified on glyph1
    glyph1.MaximumNumberOfSamplePoints = 10000

    # update the view to ensure updated data information
    renderView1.Update()

    # change scalar bar placement
    alpha_combinedLUTColorBar.Set(
        Position=[0.879154078549849, 0.3583115108917509],
        ScalarBarLength=0.32999999999999935,
    )


    # hide data in view
    Hide(validation_doste_combinedvtu, renderView1)

    # hide data in view
    Hide(streamTracer_current, renderView1)

    # layout/tab size in pixels
    layout1.SetSize(993, 706)

    # current camera placement for renderView1
    renderView1.Set(
        CameraPosition=[58.74149269073567, -105.91797844910512, 207.98796445502148],
        CameraFocalPoint=[-22.22441909870438, -68.10289751255294, 46.918200036593994],
        CameraViewUp=[-0.015152410940450833, 0.9716964855848043, 0.2357463559523924],
        CameraParallelScale=102.19338444140465,
    )

    # save screenshot
    SaveScreenshot(filename=os.path.join(png_output_path, f'doste_{family_name}_slice.png'), viewOrLayout=renderView1, location=16, ImageResolution=[993, 706], TransparentBackground=0)

    # set active source
    SetActiveSource(validation_doste_combinedvtu)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=glyph1.GlyphType)

    # show data in view
    validation_doste_combinedvtuDisplay = Show(validation_doste_combinedvtu, renderView1, 'UnstructuredGridRepresentation')

    # hide data in view
    Hide(validation_doste_combinedvtu, renderView1)

    # Delete the glyph for this iteration before creating the next one
    Delete(glyph1)
    del glyph1
    
    # Show streamlines and save second screenshot for this fiber family
    # set active source
    SetActiveSource(streamTracer_current)

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=streamTracer_current.SeedType)

    # show data in view
    streamTracer_currentDisplay = Show(streamTracer_current, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    streamTracer_currentDisplay.SetScalarBarVisibility(renderView1, True)

    # set active source
    SetActiveSource(validation_doste_combinedvtu)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=streamTracer_current.SeedType)

    # show data in view
    validation_doste_combinedvtuDisplay = Show(validation_doste_combinedvtu, renderView1, 'UnstructuredGridRepresentation')

    # hide data in view
    Hide(slice1, renderView1)

    # layout/tab size in pixels
    layout1.SetSize(993, 706)

    # current camera placement for renderView1
    renderView1.Set(
        CameraPosition=[52.04115429873035, -45.38613870538731, 216.54244474035633],
        CameraFocalPoint=[-21.838666915893594, -60.56520140171051, 48.49395084381101],
        CameraViewUp=[-0.11136472699136057, 0.9929446544942593, -0.040728499768785655],
        CameraParallelScale=102.19338444140465,
    )

    # save screenshot
    SaveScreenshot(filename=os.path.join(png_output_path, f'doste_{family_name}.png'), viewOrLayout=renderView1, location=16, ImageResolution=[993, 706], TransparentBackground=0)
    
    # Delete the stream tracer for this iteration before creating the next one
    Delete(streamTracer_current)
    del streamTracer_current
    

# Final screenshots outside the loop
# layout/tab size in pixels
layout1.SetSize(993, 706)

# current camera placement for renderView1
renderView1.Set(
    CameraPosition=[52.04115429873035, -45.38613870538731, 216.54244474035633],
    CameraFocalPoint=[-21.838666915893594, -60.56520140171051, 48.49395084381101],
    CameraViewUp=[-0.11136472699136057, 0.9929446544942593, -0.040728499768785655],
    CameraParallelScale=102.19338444140465,
)


# After loop, show streamlines and save final screenshots
# set active source
SetActiveSource(streamTracer1)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(993, 706)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.Set(
    CameraPosition=[52.04115429873035, -45.38613870538731, 216.54244474035633],
    CameraFocalPoint=[-21.838666915893594, -60.56520140171051, 48.49395084381101],
    CameraViewUp=[-0.11136472699136057, 0.9929446544942593, -0.040728499768785655],
    CameraParallelScale=102.19338444140465,
)



##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/nightly/python/
##--------------------------------------------