#include <vtkActor.h>
#include <vtkAreaPicker.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPlanes.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>
#include <vtkVertexGlyphFilter.h>

#if VTK_VERSION_NUMBER >= 89000000000ULL
#define VTK890 1
#endif

#include <vtkBYUReader.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtksys/SystemTools.hxx>

#define VTKISRBP_ORIENT 0
#define VTKISRBP_SELECT 1

namespace {
// Define interaction style
class HighlightInteractorStyle : public vtkInteractorStyleRubberBandPick
{
public:
  static HighlightInteractorStyle* New();
  vtkTypeMacro(HighlightInteractorStyle, vtkInteractorStyleRubberBandPick);

  HighlightInteractorStyle() : vtkInteractorStyleRubberBandPick()
  {
    this->SelectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    this->SelectedActor = vtkSmartPointer<vtkActor>::New();
    this->SelectedActor->SetMapper(SelectedMapper);
  }

  virtual void OnLeftButtonUp() override
  {
    // Forward events
    vtkInteractorStyleRubberBandPick::OnLeftButtonUp();

    if (this->CurrentMode == VTKISRBP_SELECT)
    {
      vtkNew<vtkNamedColors> colors;

      vtkPlanes* frustum =
          static_cast<vtkAreaPicker*>(this->GetInteractor()->GetPicker())
              ->GetFrustum();

      vtkNew<vtkExtractPolyDataGeometry> extractPolyDataGeometry;
      extractPolyDataGeometry->SetInputData(this->PolyData);
      extractPolyDataGeometry->SetImplicitFunction(frustum);
      extractPolyDataGeometry->Update();

      std::cout << "Extracted "
                << extractPolyDataGeometry->GetOutput()->GetNumberOfCells()
                << " cells." << std::endl;
      this->SelectedMapper->SetInputData(extractPolyDataGeometry->GetOutput());

      //        vtkIdTypeArray* ids =
      //        dynamic_cast<vtkIdTypeArray*>(selected->GetPointData()->GetArray("OriginalIds"));

      this->SelectedActor->GetProperty()->SetColor(
          colors->GetColor3d("Tomato").GetData());
      this->SelectedActor->GetProperty()->SetPointSize(5);

      this->GetInteractor()
          ->GetRenderWindow()
          ->GetRenderers()
          ->GetFirstRenderer()
          ->AddActor(this->SelectedActor);

      this->GetInteractor()->GetRenderWindow()->Render();
      this->HighlightProp(NULL);
    }
  }

  void SetPolyData(vtkSmartPointer<vtkPolyData> polyData)
  {
    this->PolyData = polyData;
  }

private:
  vtkSmartPointer<vtkPolyData> PolyData;
  vtkSmartPointer<vtkActor> SelectedActor;
  vtkSmartPointer<vtkDataSetMapper> SelectedMapper;
};
vtkStandardNewMacro(HighlightInteractorStyle);

vtkSmartPointer<vtkPolyData> ReadPolyData(const char* fileName);
} // namespace

int main(int argc, char* argv[])
{
  // auto polyData = ReadPolyData(argc > 1 ? argv[1] : "");

  std::string strfileName ("/Users/dsvandet/Software/Private/vtk/Testing/Data/cow.vtp");
  const char * fileName = strfileName.c_str();
  vtkSmartPointer<vtkPolyData> polyData;
  vtkNew<vtkXMLPolyDataReader> reader;
  reader->SetFileName(fileName);
  reader->Update();
  polyData = reader->GetOutput();

  vtkNew<vtkNamedColors> colors;

  vtkNew<vtkIdFilter> idFilter;
  idFilter->SetInputData(polyData);
#if VTK890
  idFilter->SetCellIdsArrayName("OriginalIds");
  idFilter->SetPointIdsArrayName("OriginalIds");
#else
  idFilter->SetIdsArrayName("OriginalIds");
#endif
  idFilter->Update();

  // This is needed to convert the ouput of vtkIdFilter (vtkDataSet) back to
  // vtkPolyData
  vtkNew<vtkDataSetSurfaceFilter> surfaceFilter;
  surfaceFilter->SetInputConnection(idFilter->GetOutputPort());
  surfaceFilter->Update();

  vtkPolyData* input = surfaceFilter->GetOutput();

  // Create a mapper and actor
  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputData(polyData);
  mapper->ScalarVisibilityOff();

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  actor->GetProperty()->SetPointSize(5);
  actor->GetProperty()->SetDiffuseColor(
      colors->GetColor3d("Peacock").GetData());
  // Visualize
  vtkNew<vtkRenderer> renderer;
  renderer->UseHiddenLineRemovalOn();

  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(640, 480);
  renderWindow->SetWindowName("HighlightSelection");

  vtkNew<vtkAreaPicker> areaPicker;
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetPicker(areaPicker);
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  renderer->SetBackground(colors->GetColor3d("Tan").GetData());

  renderWindow->Render();

  vtkNew<HighlightInteractorStyle> style;
  style->SetPolyData(input);
  renderWindowInteractor->SetInteractorStyle(style);

  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
namespace {
vtkSmartPointer<vtkPolyData> ReadPolyData(const char* fileName)
{
  vtkSmartPointer<vtkPolyData> polyData;
  std::string extension =
      vtksys::SystemTools::GetFilenameLastExtension(std::string(fileName));
  if (extension == ".ply")
  {
    vtkNew<vtkPLYReader> reader;
    reader->SetFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else if (extension == ".vtp")
  {
    vtkNew<vtkXMLPolyDataReader> reader;
    reader->SetFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else if (extension == ".obj")
  {
    vtkNew<vtkOBJReader> reader;
    reader->SetFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else if (extension == ".stl")
  {
    vtkNew<vtkSTLReader> reader;
    reader->SetFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else if (extension == ".vtk")
  {
    vtkNew<vtkPolyDataReader> reader;
    reader->SetFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else if (extension == ".g")
  {
    vtkNew<vtkBYUReader> reader;
    reader->SetGeometryFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else
  {
    vtkNew<vtkSphereSource> source;
    source->SetPhiResolution(21);
    source->SetThetaResolution(40);
    source->Update();
    polyData = source->GetOutput();
  }
  return polyData;
}
} // namespace
