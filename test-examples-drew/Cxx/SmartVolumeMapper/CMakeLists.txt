
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(SmartVolumeMapper)

find_package(VTK COMPONENTS 
  vtkCommonCore
  vtkCommonDataModel
  vtkIOXML
  vtkImagingCore
  vtkImagingHybrid
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  vtkRenderingVolumeOpenGL2
  QUIET
)

if (NOT VTK_FOUND)
  message("Skipping SmartVolumeMapper: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(SmartVolumeMapper MACOSX_BUNDLE SmartVolumeMapper.cxx )
  target_link_libraries(SmartVolumeMapper PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(SmartVolumeMapper MACOSX_BUNDLE SmartVolumeMapper.cxx )
  target_link_libraries(SmartVolumeMapper PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS SmartVolumeMapper
    MODULES ${VTK_LIBRARIES}
    )
endif ()
