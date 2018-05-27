#define _USE_MATH_DEFINES
#include <math.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/unproject_onto_mesh.h>
#include <Eigen/Sparse>
#include <igl/edges.h>
#include <igl/cat.h>
#include <igl/invert_diag.h>

#include "main.h"

// Setup constants
typedef Eigen::Triplet<double> T;
int HANDLE_MODE = 0, NO_MODE = 1;
Eigen::RowVector3d RED(255.0/255.0, 0.0/255.0, 0.0/255.0);
Eigen::RowVector3d WHITE(255.0/255.0, 255.0/255.0, 255.0/255.0);
Eigen::RowVector3d GREY(125.0/255.0, 125.0/255.0, 125.0/255.0);
Eigen::RowVector3d GOLD(255.0/255.0, 228.0/255.0, 58.0/255.0);
Eigen::RowVector3d GREEN(70./255.,252./255.,167./255.);

// Setup global variables
Eigen::MatrixXd C, V, origV, boxV(8, 3), AA, AN, ANK, H, K;
Eigen::MatrixXi E, F, N, boxE(12, 2);
Eigen::RowVector3d boundMin, boundMax, baryC, prevHandle, handle, posROI;
Eigen::VectorXd numberNeighbours;
std::string offModel = "cow";
std::vector<int> v_roi;
std::vector<T> m_coeffs, c_coeffs, l_coeffs;

float boundDiag;
bool firstRun = true, cotanMode = true, simpleCostMode;
int fid, mode, discretize, prevFid, handleInd;

igl::opengl::glfw::Viewer viewer;
igl::opengl::glfw::imgui::ImGuiMenu menu;

void ComputeAreasAndAngles() {
  double x, y, z;

  // Initialize storage matrices:
  //  AA -> Areas, AN -> Alpha and Beta Angles, ANK -> Theta Angles
  AA = Eigen::MatrixXd::Zero(V.rows(), 1);
  AN = Eigen::MatrixXd::Zero(F.rows(), V.cols());
  ANK = Eigen::MatrixXd::Zero(V.rows(), 1);

  for (int f = 0; f < F.rows(); f++) {
    // Compute the magnitudes of each face edge
    x = ((V.row(F.row(f)[1])) - (V.row(F.row(f)[2]))).norm();
    y = ((V.row(F.row(f)[0])) - (V.row(F.row(f)[2]))).norm();
    z = ((V.row(F.row(f)[0])) - (V.row(F.row(f)[1]))).norm();

    // Compute and store the angles of each constituent face angle
    AN.row(f)[0] = acos((pow(y, 2) + pow(z, 2) - pow(x, 2)) / (2 * y * z));
    AN.row(f)[1] = acos((pow(x, 2) + pow(z, 2) - pow(y, 2)) / (2 * x * z));
    AN.row(f)[2] = acos((pow(x, 2) + pow(y, 2) - pow(z, 2)) / (2 * x * y));

    //  Separately store the vertex co-joint angles
    ANK(F(f, 0), 0) = ANK(F(f, 0), 0) + AN.row(f)[0];
    ANK(F(f, 1), 0) = ANK(F(f, 1), 0) + AN.row(f)[1];
    ANK(F(f, 2), 0) = ANK(F(f, 2), 0) + AN.row(f)[2];

    // Compute and store the areas of the mesh faces
    AA(F(f, 0), 0) = AA(F(f, 0), 0) + (0.5 * y * z * sin(AN.row(f)[0])) / 3;
    AA(F(f, 1), 0) = AA(F(f, 1), 0) + (0.5 * x * z * sin(AN.row(f)[1])) / 3;
    AA(F(f, 2), 0) = AA(F(f, 2), 0) + (0.5 * x * y * sin(AN.row(f)[2])) / 3;
  }
}

void ComputeCotanAngles() {
  double x, y, z;

  for (int f = 0; f < F.rows() ; f++) {
    // Compute the magnitudes of each face edge
    x = ( (V.row(F.row(f)[1])) - (V.row(F.row(f)[2])) ).norm();
    y = ( (V.row(F.row(f)[0])) - (V.row(F.row(f)[2])) ).norm();
    z = ( (V.row(F.row(f)[0])) - (V.row(F.row(f)[1])) ).norm();

    // Compute and store all M diagonal sparse entries
    m_coeffs.push_back(T(F(f, 0), F(f, 0), ( 0.5 * y * z * sin(AN.row(f)[0]) ) / 3 ));
    m_coeffs.push_back(T(F(f, 1), F(f, 1), ( 0.5 * x * z * sin(AN.row(f)[1]) ) / 3 ));
    m_coeffs.push_back(T(F(f, 2), F(f, 2), ( 0.5 * x * y * sin(AN.row(f)[2]) ) / 3 ));

    // Compute and store all C sparse entries
    c_coeffs.push_back(T(F(f, 0), F(f, 0), -0.5 * ((cos(AN.row(f)[1]) / sin(AN.row(f)[1])) + (cos(AN.row(f)[2]) / sin(AN.row(f)[2])) )));
    c_coeffs.push_back(T(F(f, 1), F(f, 1), -0.5 * ((cos(AN.row(f)[0]) / sin(AN.row(f)[0])) + (cos(AN.row(f)[2]) / sin(AN.row(f)[2])) )));
    c_coeffs.push_back(T(F(f, 2), F(f, 2), -0.5 * ((cos(AN.row(f)[0]) / sin(AN.row(f)[0])) + (cos(AN.row(f)[1]) / sin(AN.row(f)[1])) )));
    c_coeffs.push_back(T(F(f, 0), F(f, 1), 0.5 * (cos(AN.row(f)[2]) / sin(AN.row(f)[2])) ));
    c_coeffs.push_back(T(F(f, 1), F(f, 0), 0.5 * (cos(AN.row(f)[2]) / sin(AN.row(f)[2])) ));
    c_coeffs.push_back(T(F(f, 1), F(f, 2), 0.5 * (cos(AN.row(f)[0]) / sin(AN.row(f)[0])) ));
    c_coeffs.push_back(T(F(f, 2), F(f, 1), 0.5 * (cos(AN.row(f)[0]) / sin(AN.row(f)[0])) ));
    c_coeffs.push_back(T(F(f, 0), F(f, 2), 0.5 * (cos(AN.row(f)[1]) / sin(AN.row(f)[1])) ));
    c_coeffs.push_back(T(F(f, 2), F(f, 0), 0.5 * (cos(AN.row(f)[1]) / sin(AN.row(f)[1])) ));
  }
}

Eigen::MatrixXi getNeighbours(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXi &E, Eigen::VectorXd &number){

  // initialising of matrices
  Eigen::MatrixXi N = Eigen::MatrixXi::Zero(V.rows(), V.rows());
  number = Eigen::VectorXd::Zero(V.rows());

  // construct edge-matrix
  igl::edges(F, E);

  // check for neighbours
  for(int i=0; i < V.rows(); i++){
    int count = 0;
    for(int j=0; j < E.rows(); j++){
      if(E(j,0) == i){
        // store index of neighbour
        N(i,count) = E(j,1);
        count++;
      }
      else if(E(j,1) == i){
        // store index of neighbour
        N(i,count) = E(j,0);
        count++;
      }
    }
    // store number of Neighbours
    number(i) = count;
  }

  return N;
}

void loadMesh() {
  igl::readOFF(MODEL_PATH "/" + offModel + ".off", V, F);
  ComputeAreasAndAngles();
   // compute number of neighbours and store neighbours per vertex
  N = getNeighbours(V, F, E, numberNeighbours);

  viewer.data().clear();
  viewer.data().set_mesh(V, F);
  C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
  viewer.data().set_colors(C);
  viewer.core.align_camera_center(V, F);

  boundMin = V.colwise().minCoeff();
  boundMax = V.colwise().maxCoeff();
  boundDiag = (boundMax - boundMin).norm();
  handle << 0, 0, 0;

  if (firstRun) {
    origV = V;
    firstRun = !firstRun;
  }
}

bool isContained(Eigen::RowVector3d point) {
  bool inX = false;
  bool inY = false;
  bool inZ = false;

  // Confirm x is more than box lefts and less than box rights
  if (
    point(0) > boxV(0, 0) && point(0) > boxV(3, 0) && point(0) > boxV(4, 0) && point(0) > boxV(7, 0) &&
    point(0) < boxV(1, 0) && point(0) < boxV(2, 0) && point(0) < boxV(5, 0) && point(0) < boxV(6, 0)
    )
      inX = true;

  // Confirm y is more than box bottoms and less than box tops
  if (
    point(1) > boxV(0, 1) && point(1) > boxV(1, 1) && point(1) > boxV(4, 1) && point(1) > boxV(5, 1) &&
    point(1) < boxV(2, 1) && point(1) < boxV(3, 1) && point(1) < boxV(6, 1) && point(1) < boxV(7, 1)
    )
      inY = true;

  // Confirm z is more than box backs and less than box fronts
  if (
    point(2) > boxV(0, 2) && point(2) > boxV(1, 2) && point(2) > boxV(2, 2) && point(2) > boxV(3, 2) &&
    point(2) < boxV(4, 2) && point(2) < boxV(5, 2) && point(2) < boxV(6, 2) && point(2) < boxV(7, 2)
    )
      inZ = true;

  return (inX && inY && inZ);
}

void drawBoundBox(Eigen::RowVector3d boxMax, Eigen::RowVector3d boxMin, Eigen::RowVector3d posVec, float scaleL, float scaleW, float scaleH) {
  // TODO: Scale variables scaleL and scaleW should be used to resize the bounding box

  // If no handle has been fixed on the mesh, then return
  if (handle.norm() == 0.0) {
    return;
  }

  // Corners of the bounding box; using RH coordinate frame
  // Ensure that the handle is within the prototyped bound; x, y, z -> z, y, x
  boxV <<
  (boxMin(0) * scaleL) + posVec(0), (boxMin(1) * scaleW) + posVec(1), (boxMin(2) * scaleH) + posVec(2), // 0 - back bottom left
  (boxMax(0) * scaleL) + posVec(0), (boxMin(1) * scaleW) + posVec(1), (boxMin(2) * scaleH) + posVec(2), // 1 - back bottom right
  (boxMax(0) * scaleL) + posVec(0), (boxMax(1) * scaleW) + posVec(1), (boxMin(2) * scaleH) + posVec(2), // 2 - back top right
  (boxMin(0) * scaleL) + posVec(0), (boxMax(1) * scaleW) + posVec(1), (boxMin(2) * scaleH) + posVec(2), // 3 - back top left
  (boxMin(0) * scaleL) + posVec(0), (boxMin(1) * scaleW) + posVec(1), (boxMax(2) * scaleH) + posVec(2), // 4 - front bottom left
  (boxMax(0) * scaleL) + posVec(0), (boxMin(1) * scaleW) + posVec(1), (boxMax(2) * scaleH) + posVec(2), // 5 - front bottom right
  (boxMax(0) * scaleL) + posVec(0), (boxMax(1) * scaleW) + posVec(1), (boxMax(2) * scaleH) + posVec(2), // 6 - front top right
  (boxMin(0) * scaleL) + posVec(0), (boxMax(1) * scaleW) + posVec(1), (boxMax(2) * scaleH) + posVec(2), // 7 - front top left

  // Edges of the bounding box
  boxE <<
    0, 1,
    1, 2,
    2, 3,
    3, 0,
    4, 5,
    5, 6,
    6, 7,
    7, 4,
    0, 4,
    1, 5,
    2, 6,
    7, 3;

  viewer.data().set_points(handle, GOLD);
  viewer.data().add_points(boxV, RED);
  viewer.data().set_edges(boxV, boxE, GREY);
}

bool keyDownCallback(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers) {
  switch(key)
  {
    case 'C':
    {
      std::cout << "No Mode." << std::endl;
      mode = NO_MODE;
      return false;
    }
  }
  return true;
}

bool mouseDownCallback(igl::opengl::glfw::Viewer& viewer, int, int) {
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  float minDist;

  if (igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view * viewer.core.model,
    viewer.core.proj, viewer.core.viewport, V, F, fid, baryC) && mode != NO_MODE)
  {
    C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);

    // Compute nearest vertex to barycenter point baryC
    // BUG: The camera align function centers the frame such that the nearest vertex is not regular across the mesh.
    // std::cout << baryC << std::endl << V.row(F.row(fid)[0]) << " " << (baryC - V.row(F.row(fid)[0])).norm() << std::endl << V.row(F.row(fid)[1])  << " " << (baryC - V.row(F.row(fid)[1])).norm() << std::endl << V.row(F.row(fid)[2]) << " " << (baryC - V.row(F.row(fid)[2])).norm() << std::endl;
    for (int c = 0; c < V.cols(); c++) {
      float bcDist = (baryC - V.row(F.row(fid)[c])).norm();

      if (c == 0) {
        minDist = bcDist;
      }
      else if (bcDist < minDist)
      {
        minDist = bcDist;
        handleInd = F.row(fid)[c];
        handle << V.row(handleInd);
      }
    }

    if (!firstRun) {
      // Paint previously selected face WHITE
      C.row(prevFid) = WHITE;
    }

    // Paint currently selected face GREY
    C.row(fid) = GREY;
    viewer.data().set_colors(C);
    viewer.data().set_points(handle, GOLD);

    prevFid = fid;
    prevHandle = handle;

    return true;
  }
  return false;
}

void computeMeanCurvature(Eigen::SparseMatrix<double> & L, Eigen::MatrixXd & H) {
  //    Compute Mean Curvature using the formula H = ||UL*x||/2
  Eigen::MatrixXd Lx(V.rows(), V.cols());

  for (int c = 0; c < V.cols(); c++) {
      Lx.col(c) = L * V.col(c);
  }

  H = Lx.rowwise().norm()/2;
}

void computeGaussianCurvature(bool & normalized, Eigen::MatrixXd & K) {
  K = Eigen::MatrixXd::Zero(V.rows(), 1);

  // Check if normalization flag is set
  // Then compute Gaussian curvature appropriately
  if (normalized == true) {
    for (int r = 0; r < V.rows(); r++) {
      K(r) = (2 * M_PI - ANK(r))/AA(r);
    }
  } else {
    for (int r = 0; r < V.rows(); r++) {
      K(r) = (2 * M_PI - ANK(r));
    }
  }
}

void computeLaplace(Eigen::MatrixXd &V, Eigen::SparseMatrix<double> &Laplace, Eigen::VectorXd &numberNeighbours, Eigen::MatrixXi &N, bool cotangentMode){

    if (!cotangentMode) {
      std::vector<T> tripletList;

      for(int i=0; i < V.rows(); i++) {
          // fill diagonal of Laplacian
          tripletList.push_back(T(i, i, 1));
          // fill rows of Laplacian according to number of neighbours
          for(int j=0; j < numberNeighbours(i); j++){
              tripletList.push_back(T(i, N(i,j), -1/numberNeighbours(i)));
          }
      }
      // set sparse matrix
      Laplace.setFromTriplets(tripletList.begin(), tripletList.end());
    } else {
      Eigen::SparseMatrix<double> LM(V.rows(), V.rows()), LC(V.rows(), V.rows()), LM_inv(V.rows(), V.rows());

      ComputeCotanAngles();
      LM.setFromTriplets(m_coeffs.begin(), m_coeffs.end());
      LC.setFromTriplets(c_coeffs.begin(), c_coeffs.end());
      igl::invert_diag(LM, LM_inv);
      Laplace = LM_inv * LC;
    }
}

void laplacianEdit(Eigen::RowVector3d newHandle, Eigen::RowVector3d oldHandle, int handleVInd, std::vector<int> vROI, Eigen::MatrixXd &V) {
  // std::cout << vROI.size() << std::endl << handleVInd << std::endl << oldHandle << std::endl << newHandle << std::endl;

  // Homogeneous coordinates of V
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(V.rows());
  Eigen::MatrixXd Vhom(V.rows(), V.cols()+1);
  Vhom << V, ones;

  // compute uniform Laplace
  Eigen::SparseMatrix<double> Laplace(V.rows(), V.rows());
  computeLaplace(V, Laplace, numberNeighbours, N, cotanMode);

  // compute Laplacian Coordinates
  Eigen::SparseMatrix<double> Delta;

  Eigen::SparseMatrix<double> Vsparse = V.sparseView();
  Delta = Laplace * Vsparse;

  std::cout << Delta.row(0) << std::endl;
  std::cout << Delta.row(1) << std::endl;

  Eigen::SparseMatrix<double> Ones = ones.sparseView();
  Eigen::SparseMatrix<double> HomDelta;
  igl::cat(2, Delta, Ones, HomDelta);

  // Compute positions VPrime
  Eigen::Vector3d handleMovement;
  handleMovement = newHandle - oldHandle;
  // Apply handleMovement to each vertex in ROI
  Eigen::VectorXd simulatedIndices(vROI.size());
  Eigen::MatrixXd Vprime(V.rows(), V.cols());
  Vprime = V;
  for(int i = 0; i < vROI.size(); i++){
    Eigen::Vector3d currentRow = Vprime.row(vROI[i]).transpose();
    currentRow += handleMovement;
    Vprime.row(vROI[i]) = currentRow.transpose();
    simulatedIndices(i) = vROI[i];
  }

  Eigen::SparseMatrix<double> DeltaPrime;
  Eigen::SparseMatrix<double> VPrimeSparse = Vprime.sparseView();
  DeltaPrime = Laplace * VPrimeSparse;

  Eigen::SparseMatrix<double> HomDeltaPrime;
  igl::cat(2, DeltaPrime, Ones, HomDeltaPrime);

  double error = 1;
  double pre_error = 0;
  //while(error > pre_error) {
  for(int h = 0; h < 1; h++){
      pre_error = error;
      error = 0;
      for (int i = 0; i < simulatedIndices.rows(); i++) {
          int count = 0;
          Eigen::VectorXd b_0 = Eigen::VectorXd::Zero((numberNeighbours(simulatedIndices(i)) + 1) * 3);
          b_0(count) = Vprime(simulatedIndices(i), 0);
          count++;
          b_0(count) = Vprime(simulatedIndices(i), 1);
          count++;
          b_0(count) = Vprime(simulatedIndices(i), 2);
          count++;
          for (int j = 0; j < numberNeighbours(simulatedIndices(i)); j++) {
              b_0(count) = Vprime(N(simulatedIndices(i), j), 0);
              count++;
              b_0(count) = Vprime(N(simulatedIndices(i), j), 1);
              count++;
              b_0(count) = Vprime(N(simulatedIndices(i), j), 2);
              count++;
          }

          //build A
          Eigen::MatrixXd A_0 = Eigen::MatrixXd::Zero((numberNeighbours(simulatedIndices(i)) + 1) * 3, 7);
          count = 0;
          //add vertix i
          A_0(count, 0) = V(simulatedIndices(i), 0);
          A_0(count, 2) = V(simulatedIndices(i), 2);
          A_0(count, 3) = -V(simulatedIndices(i), 1);
          A_0(count, 4) = 1;
          count++;
          A_0(count, 0) = V(simulatedIndices(i), 1);
          A_0(count, 1) = -V(simulatedIndices(i), 2);
          A_0(count, 3) = V(simulatedIndices(i), 0);
          A_0(count, 5) = 1;
          count++;
          A_0(count, 0) = V(simulatedIndices(i), 2);
          A_0(count, 1) = V(simulatedIndices(i), 1);
          A_0(count, 2) = -V(simulatedIndices(i), 0);
          A_0(count, 6) = 1;
          count++;

          //add neighbours
          for (int j = 0; j < numberNeighbours(simulatedIndices(i)); j++) {
              A_0(count, 0) = V(N(simulatedIndices(i), j), 0);
              A_0(count, 2) = V(N(simulatedIndices(i), j), 2);
              A_0(count, 3) = -V(N(simulatedIndices(i), j), 1);
              A_0(count, 4) = 1;
              count++;
              A_0(count, 0) = V(N(simulatedIndices(i), j), 1);
              A_0(count, 1) = -V(N(simulatedIndices(i), j), 2);
              A_0(count, 3) = V(N(simulatedIndices(i), j), 0);
              A_0(count, 5) = 1;
              count++;
              A_0(count, 0) = V(N(simulatedIndices(i), j), 2);
              A_0(count, 1) = V(N(simulatedIndices(i), j), 1);
              A_0(count, 2) = -V(N(simulatedIndices(i), j), 0);
              A_0(count, 6) = 1;
              count++;
          }

          Eigen::MatrixXd Coeff;
          Coeff = A_0.transpose() * A_0;
          Coeff = Coeff.inverse();
          Coeff = Coeff * A_0.transpose();
          Coeff = Coeff * b_0;
          //cout << "Coeff: " << Coeff.transpose() << endl;

          // build T
          Eigen::MatrixXd T(4, 4);
          T << Coeff(0, 0), -Coeff(3, 0), Coeff(2, 0), Coeff(4, 0), Coeff(3, 0), Coeff(0, 0), -Coeff(1, 0),
                  Coeff(5, 0), -Coeff(2, 0), Coeff(1, 0), Coeff(0, 0), Coeff(6, 0), 0, 0, 0, 1;

          // Compute error
          // Eigen::VectorXd currentError = T * HomDelta.row(0).transpose();
          // currentError -= HomDeltaPrime.row(0).transpose();
          // error += currentError.squaredNorm();

          // update V
          Eigen::VectorXd currentRow(4);
          currentRow(0) = V(simulatedIndices(i),0);
          currentRow(1) = V(simulatedIndices(i),1);
          currentRow(2) = V(simulatedIndices(i),2);
          currentRow(3) = 1;
          currentRow = T * currentRow;

          V(simulatedIndices(i),0) = currentRow(0);
          V(simulatedIndices(i),1) = currentRow(1);
          V(simulatedIndices(i),2) = currentRow(2);
      }
      std::cout << "error: " << error << std::endl;
  }
  viewer.data().set_mesh(V, F);
}

void initGUI() {

  // Attach a menu plugin
  viewer.plugins.push_back(& menu);

  // Add content to the default menu window
  menu.draw_viewer_menu_func = [&]() {
    // Draw parent menu content
    // menu.draw_viewer_menu();

    static float lScale = 1.0, wScale = 1.0, hScale = 1.0, handleX = 0.00f, handleY = 0.00f, handleZ = 0.00f;

     // Add new group
    if (ImGui::CollapsingHeader("CW3", ImGuiTreeNodeFlags_DefaultOpen))
    {
      // Display mesh properties and messages
      ImGui::Text("Bounding diagonal: %f", boundDiag);
      ImGui::Text("Vertices in Mesh: %i", V.rows());
      ImGui::Text("Faces in Mesh: %i", F.rows());
      ImGui::Text("Vertices in ROI: %i", v_roi.size());
      if (handle.norm() == 0)
      {
        ImGui::TextDisabled("No Handle placed on mesh!");
      }
      ImGui::Separator();

      // Handle Config
      ImGui::Text("Handle Config");
      if (ImGui::Button("Position Handle", ImVec2(-1,0))){
        mode = HANDLE_MODE;
      }
      if (ImGui::Button("Fix Handle", ImVec2(-1,0))) {
        mode = NO_MODE;
        handleX = float(handle(0));
        handleY = float(handle(1));
        handleZ = float(handle(2));
        posROI = handle/2.0;
        drawBoundBox(boundMax, boundMin, posROI, 1.0, 1.0, 1.0);
      }
      ImGui::Separator();

      // ROI Controls
      ImGui::Text("ROI Controls");

      ImGui::SliderFloat("Length", &lScale, 0.1f, 1.0f);
      ImGui::SliderFloat("Width", &wScale, 0.1f, 1.0f);
      ImGui::SliderFloat("Height", &hScale, 0.1f, 1.0f);

      if (ImGui::Button("Increase X")) {
        posROI(0) = posROI(0) + (0.01 * boundDiag);
        drawBoundBox(boundMax, boundMin, posROI, lScale, wScale, hScale);
      }
      ImGui::SameLine();
      if (ImGui::Button("Decrease X")) {
        posROI(0) = posROI(0) - (0.01 * boundDiag);
        drawBoundBox(boundMax, boundMin, posROI, lScale, wScale, hScale);
      }

      if (ImGui::Button("Increase Y")) {
        posROI(1) = posROI(1) + (0.01 * boundDiag);
        drawBoundBox(boundMax, boundMin, posROI, lScale, wScale, hScale);
      }
      ImGui::SameLine();
      if (ImGui::Button("Decrease Y")) {
        posROI(1) = posROI(1) - (0.01 * boundDiag);
        drawBoundBox(boundMax, boundMin, posROI, lScale, wScale, hScale);
      }

      if (ImGui::Button("Increase Z")) {
        posROI(2) = posROI(2) + (0.01 * boundDiag);
        drawBoundBox(boundMax, boundMin, posROI, lScale, wScale, hScale);
      }
      ImGui::SameLine();
      if (ImGui::Button("Decrease Z")) {
        posROI(2) = posROI(2) - (0.01 * boundDiag);
        drawBoundBox(boundMax, boundMin, posROI, lScale, wScale, hScale);
      }
      ImGui::Separator();

      if (ImGui::Button("Refresh ROI Bound", ImVec2(-1,0))) {
        // loadMesh();
        drawBoundBox(boundMax, boundMin, posROI, lScale, wScale, hScale);
      }

      if (ImGui::Button("Compute ROI", ImVec2(-1,0))) {
        v_roi.clear();
        for (unsigned int r = 0; r < V.rows(); r++) {
          if (isContained(V.row(r))) {
            v_roi.push_back(r);
          }
        }
      }
      // Shade ROI for Visual Inspection
      if (ImGui::Button("Shade ROI", ImVec2(-1,0))) {
        Eigen::MatrixXd VC;
        VC = Eigen::MatrixXd::Constant(V.rows(), 3, 1);

        for (int i=0; i < v_roi.size(); i++) {
          VC.row(v_roi[i]) = GREEN;
        }

        viewer.data().set_colors(VC);
      }
      ImGui::Separator();

      // Manipulate Handle Controls
      ImGui::Text("Handle Controls");
      ImGui::DragFloat("X-Axis", & handleX, 0.001f);
      ImGui::DragFloat("Y-Axis", & handleY, 0.001f);
      ImGui::DragFloat("Z-Axis", & handleZ, 0.001f);
      ImGui::Separator();

      // Laplacian Edit Controls
      ImGui::Text("Laplacian Editing");
      ImGui::Checkbox("Use Cotan", &cotanMode);
      ImGui::Checkbox("Simple Cost", &simpleCostMode);

      if (ImGui::Button("Execute", ImVec2(-1,0))) {
        laplacianEdit(Eigen::RowVector3d(double(handleX), double(handleY), double(handleZ)), handle, handleInd, v_roi, V);
      }
      ImGui::Separator();

      if (ImGui::Button("Apply Uniform Mean Curvature", ImVec2(-1,0))) {
        Eigen::SparseMatrix<double> L(V.rows(), V.rows());
        Eigen::MatrixXd COLOR;

        computeLaplace(V, L, numberNeighbours, N, cotanMode);
        computeMeanCurvature(L, H);
        igl::jet(H, true, COLOR);
        viewer.data().set_colors(COLOR);
        std::cout << "Applied Mean Curvature" << std::endl;
      }
      if (ImGui::Button("Apply Uniform Gaussian Curvature", ImVec2(-1,0))) {
        Eigen::MatrixXd COLOR;
        bool normalized = true;

        computeGaussianCurvature(normalized, K);
        igl::jet(K, true, COLOR);
        viewer.data().set_colors(COLOR);
        std::cout << "Applied Gaussian Curvature" << std::endl;
      }
      ImGui::Separator();

      // Settings and Config
      if (ImGui::CollapsingHeader("Settings", true)) {
        // Load new mesh
        ImGui::InputText("Enter mesh name", offModel);
        if (ImGui::Button("Reload", ImVec2(-1,0))) {
            firstRun = true;
            loadMesh();
        }
        // Reset to the original mesh
        if (ImGui::Button("Reset", ImVec2(-1,0))) {
            V = origV;
            viewer.data().set_mesh(V, F);
        }
      }

      // Exit application
      if (ImGui::Button("Exit", ImVec2(-1,0))) {
          exit(0);
      }
    }
  };
}

int main(int argc, char *argv[]) {
  loadMesh();

  viewer.callback_mouse_down = & mouseDownCallback;
  viewer.callback_key_down = & keyDownCallback;

  initGUI();

  viewer.data().point_size = 15;
  viewer.launch();

  return 0;
}
