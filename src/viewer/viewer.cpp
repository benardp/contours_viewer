#include "viewer.h"
#include "svg.h"

#include <algorithm>
#include <geogram_gfx/third_party/ImGui/imgui.h>
#include <geogram_gfx/third_party/ImGui/imgui_internal.h>

Viewer::Viewer(int argc, char **argv, const std::string &usage)
    : SimpleMeshApplication(argc, argv, usage) {
  name_ = "Contour viewer";
  lighting_ = true;
  left_pane_visible_ = true;
  right_pane_visible_ = true;
  console_visible_ = false;
  show_mesh_ = true;
  show_surface_sides_ = true;
  algo_pane_visible_ = true;
  surface_color_ = GEO::vec4f(frontFaceColor.x(), frontFaceColor.y(),
                              frontFaceColor.z(), 1.0f);
  surface_color_2_ =
      GEO::vec4f(backFaceColor.x(), backFaceColor.y(), backFaceColor.z(), 1.0f);
  mesh_color_ = GEO::vec4f(wireframeColor.x(), wireframeColor.y(),
                           wireframeColor.z(), 1.0f);
  background_color_1_ = GEO::vec4f(1.0f, 1.0f, 1.0f, 1.0f);
  background_color_2_ = GEO::vec4f(1.0f, 1.0f, 1.0f, 1.0f);

  save_dialog_.set_default_filename("out.svg");

  m_mesh = m_backup_mesh = nullptr;

  m_normals = false;
  m_surface = true;
  m_contours = true;
  m_draw_camera = false;
  m_2D_intersections = true;
  m_3D_intersections = true;
  m_curtain_folds = true;
  m_compute_visibility = false;
  m_build_view_graph = true;
  m_visible_only = false;
  m_boundaries = true;
  m_surfintersect = true;
  m_ray_intersections = false;
  m_surface_intersections = true;
  m_debug_points = true;
  m_compute_chains = true;

  m_ContourMode = 0;
  m_prev_ContourMode = -1;
  m_CameraMode = 0;
  m_prev_CameraMode = 0;

#ifdef GEO_OS_EMSCRIPTEN
  m_linewidth = 2.0f;
#else
  m_linewidth = 1.0f;
#endif
  m_pointSize = 8.0f;
  m_normal_length = 0.01f;
}

bool Viewer::load(const std::string &filename) {
  if (!SimpleMeshApplication::load(filename))
    return false;

  current_file_ = filename.substr(0, filename.length() - 3) + "svg";

  assert(mesh_.vertices.nb() > 0);
  assert(mesh_.facets.nb() > 0);
  mesh_.facets.triangulate();

  mesh_.vertices.set_double_precision();

  if (m_mesh)
    delete m_mesh;
  m_mesh = nullptr;
  if (m_backup_mesh)
    delete m_backup_mesh;
  m_backup_mesh = new Mesh();

  std::vector<Vertex> vertices;
  vertices.resize(mesh_.vertices.nb());
  for (GEO::index_t i = 0; i < mesh_.vertices.nb(); i++) {
    const GEO::vec3 &p = mesh_.vertices.point(i);
    vertices[i] = m_backup_mesh->add_vertex(Vector3f(p.x, p.y, p.z));
    assert(m_backup_mesh->is_valid(vertices[i]));
  }

  for (GEO::index_t f = 0; f < mesh_.facets.nb(); f++) {
    auto face =
        m_backup_mesh->add_triangle(vertices[mesh_.facets.vertex(f, 0)],
                                    vertices[mesh_.facets.vertex(f, 1)],
                                    vertices[mesh_.facets.vertex(f, 2)]);
    if (!m_backup_mesh->is_valid(face)) {
      GEO::Logger::out("I/O") << "Unsupported mesh" << std::endl;
      delete m_backup_mesh;
      m_backup_mesh = nullptr;
      return false;
    }
  }

  m_backup_mesh->init();

  m_backup_mesh->tagConcaveEdges();

  mesh_.vertices.set_single_precision();
  return true;
}

std::string Viewer::supported_write_file_extensions() { return "svg"; }

bool Viewer::save(const std::string &filename) {
  Vector2i mSize;
  glup_viewer_get_screen_size(&mSize(1), &mSize(0));
  SVG svgWriter(filename, mSize);
  int start_index = 0;
  Matrix3Xf contours = m_mesh->get_contours(VISIBLE_AND_INVISIBLE);
  int offset = (m_ContourMode == INTERPOLATED_CONTOUR ? 0 : 1);
  std::vector<Vector2f> polyline;
  for (size_t i = 0; i < m_mesh->get_chain_lengths().size(); i++) {
    polyline.clear();
    projectToViewport(contours.block(0, start_index, 3,
                                     m_mesh->get_chain_lengths()[i] + offset),
                      polyline, m_projection, m_modelview, mSize);
    svgWriter.writePolyline(polyline, 1, alphabetColors[i % 26], false);
    start_index += m_mesh->get_chain_lengths()[i] + offset;
  }
  if (m_mesh->get_boundaries_lengths().size() > 0) {
    Matrix3Xf boundaries = m_mesh->get_boundaries(false);
    start_index = 0;
    for (size_t i = 0; i < m_mesh->get_boundaries_lengths().size(); i++) {
      polyline.clear();
      projectToViewport(boundaries.block(0, start_index, 3,
                                         m_mesh->get_boundaries_lengths()[i]),
                        polyline, m_projection, m_modelview, mSize);
      svgWriter.writePolyline(polyline, 1, boundariesColor, false);
      start_index += m_mesh->get_boundaries_lengths()[i];
    }
  }
  if (m_mesh->get_surface_intersections_lengths().size() > 0) {
    Matrix3Xf surface_intersections = m_mesh->get_surface_intersections();
    start_index = 0;
    for (size_t i = 0; i < m_mesh->get_surface_intersections_lengths().size();
         i++) {
      polyline.clear();
      projectToViewport(surface_intersections.block(
                            0, start_index, 3,
                            m_mesh->get_surface_intersections_lengths()[i]),
                        polyline, m_projection, m_modelview, mSize);
      svgWriter.writePolyline(polyline, 1, surfIntersectionsColor, false);
      start_index += m_mesh->get_surface_intersections_lengths()[i];
    }
  }
  return true;
}

void Viewer::init_graphics() {
  SimpleMeshApplication::init_graphics();
  glup_viewer_add_toggle('n', &m_normals, "show normals");
  glup_viewer_add_toggle('o', &m_contours, "show contours");
  glup_viewer_add_toggle('d', &m_debug_points, "show debug points");
  glup_viewer_add_toggle('v', &m_visible_only, "visible curves");

  setStyle();

  for (unsigned int i = 0; i < 2; i++) {
    m_trans[i] = Eigen::Map<Eigen::Matrix<float, 3, 1>>(
        glup_viewer_get_scene_translation());
    m_rot[i] = Eigen::Map<Eigen::Matrix<float, 4, 1>>(
        glup_viewer_get_scene_quaternion());
    m_zoom[i] = glup_viewer_get_float(GLUP_VIEWER_ZOOM);
  }
}

void Viewer::setStyle() {
  ImGuiStyle &style = ImGui::GetStyle();
  ImVec4 *colors = style.Colors;

  colors[ImGuiCol_Text] = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
  colors[ImGuiCol_TextDisabled] = ImVec4(0.40f, 0.40f, 0.40f, 1.00f);
  colors[ImGuiCol_ChildBg] = ImVec4(0.25f, 0.25f, 0.25f, 1.00f);
  colors[ImGuiCol_WindowBg] = ImVec4(0.25f, 0.25f, 0.25f, 1.00f);
  colors[ImGuiCol_PopupBg] = ImVec4(0.25f, 0.25f, 0.25f, 1.00f);
  colors[ImGuiCol_Border] = ImVec4(0.12f, 0.12f, 0.12f, 0.71f);
  colors[ImGuiCol_BorderShadow] = ImVec4(1.00f, 1.00f, 1.00f, 0.06f);
  colors[ImGuiCol_FrameBg] = ImVec4(0.42f, 0.42f, 0.42f, 0.54f);
  colors[ImGuiCol_FrameBgHovered] = ImVec4(0.42f, 0.42f, 0.42f, 0.40f);
  colors[ImGuiCol_FrameBgActive] = ImVec4(0.56f, 0.56f, 0.56f, 0.67f);
  colors[ImGuiCol_TitleBg] = ImVec4(0.19f, 0.19f, 0.19f, 1.00f);
  colors[ImGuiCol_TitleBgActive] = ImVec4(0.22f, 0.22f, 0.22f, 1.00f);
  colors[ImGuiCol_TitleBgCollapsed] = ImVec4(0.17f, 0.17f, 0.17f, 0.90f);
  colors[ImGuiCol_MenuBarBg] = ImVec4(0.335f, 0.335f, 0.335f, 1.000f);
  colors[ImGuiCol_ScrollbarBg] = ImVec4(0.24f, 0.24f, 0.24f, 0.53f);
  colors[ImGuiCol_ScrollbarGrab] = ImVec4(0.41f, 0.41f, 0.41f, 1.00f);
  colors[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.52f, 0.52f, 0.52f, 1.00f);
  colors[ImGuiCol_ScrollbarGrabActive] = ImVec4(0.76f, 0.76f, 0.76f, 1.00f);
  colors[ImGuiCol_CheckMark] = ImVec4(0.65f, 0.65f, 0.65f, 1.00f);
  colors[ImGuiCol_SliderGrab] = ImVec4(0.52f, 0.52f, 0.52f, 1.00f);
  colors[ImGuiCol_SliderGrabActive] = ImVec4(0.64f, 0.64f, 0.64f, 1.00f);
  colors[ImGuiCol_Button] = ImVec4(0.54f, 0.54f, 0.54f, 0.35f);
  colors[ImGuiCol_ButtonHovered] = ImVec4(0.52f, 0.52f, 0.52f, 0.59f);
  colors[ImGuiCol_ButtonActive] = ImVec4(0.76f, 0.76f, 0.76f, 1.00f);
  colors[ImGuiCol_Header] = ImVec4(0.38f, 0.38f, 0.38f, 1.00f);
  colors[ImGuiCol_HeaderHovered] = ImVec4(0.47f, 0.47f, 0.47f, 1.00f);
  colors[ImGuiCol_HeaderActive] = ImVec4(0.76f, 0.76f, 0.76f, 0.77f);
  colors[ImGuiCol_Separator] = ImVec4(0.000f, 0.000f, 0.000f, 0.137f);
  colors[ImGuiCol_SeparatorHovered] = ImVec4(0.700f, 0.671f, 0.600f, 0.290f);
  colors[ImGuiCol_SeparatorActive] = ImVec4(0.702f, 0.671f, 0.600f, 0.674f);
  colors[ImGuiCol_ResizeGrip] = ImVec4(0.26f, 0.59f, 0.98f, 0.25f);
  colors[ImGuiCol_ResizeGripHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.67f);
  colors[ImGuiCol_ResizeGripActive] = ImVec4(0.26f, 0.59f, 0.98f, 0.95f);
  colors[ImGuiCol_PlotLines] = ImVec4(0.61f, 0.61f, 0.61f, 1.00f);
  colors[ImGuiCol_PlotLinesHovered] = ImVec4(1.00f, 0.43f, 0.35f, 1.00f);
  colors[ImGuiCol_PlotHistogram] = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
  colors[ImGuiCol_PlotHistogramHovered] = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
  colors[ImGuiCol_TextSelectedBg] = ImVec4(0.73f, 0.73f, 0.73f, 0.35f);
  colors[ImGuiCol_ModalWindowDimBg] = ImVec4(0.80f, 0.80f, 0.80f, 0.35f);
  colors[ImGuiCol_DragDropTarget] = ImVec4(1.00f, 1.00f, 0.00f, 0.90f);
  colors[ImGuiCol_NavHighlight] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  colors[ImGuiCol_NavWindowingHighlight] = ImVec4(1.00f, 1.00f, 1.00f, 0.70f);
  colors[ImGuiCol_NavWindowingDimBg] = ImVec4(0.80f, 0.80f, 0.80f, 0.20f);

  style.PopupRounding = 3;

  style.WindowPadding = ImVec2(4, 4);
  style.FramePadding = ImVec2(6, 4);
  style.ItemSpacing = ImVec2(6, 2);

  style.ScrollbarSize = 18;

  style.WindowBorderSize = 1;
  style.ChildBorderSize = 1;
  style.PopupBorderSize = 1;
  style.FrameBorderSize = 1;

  style.WindowRounding = 3;
  style.ChildRounding = 3;
  style.FrameRounding = 3;
  style.ScrollbarRounding = 2;
  style.GrabRounding = 3;
}

void Viewer::draw_left_pane() {
  int w, h;
  glup_viewer_get_screen_size(&w, &h);
  if (status_bar_->active()) {
    h -= (STATUS_HEIGHT() + 1);
  }
  if (console_visible_) {
    h -= (CONSOLE_HEIGHT() + 1);
  }
  h -= MENU_HEIGHT();
  if (algo_pane_visible_) {
    h /= 2;
  }

  if (fixed_layout_) {
    ImGui::SetNextWindowPos(ImVec2(0.0f, float(MENU_HEIGHT())),
                            ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(float(PANE_WIDTH()), float(h)),
                             ImGuiCond_Always);
  } else {
    ImGui::SetNextWindowPos(
        ImVec2(float(MENU_HEIGHT()), 2.0f * float(MENU_HEIGHT())),
        ImGuiCond_Once);
    ImGui::SetNextWindowSize(ImVec2(float(PANE_WIDTH()), float(h) / 2.0f),
                             ImGuiCond_Once);
  }

  SimpleMeshApplication::draw_viewer_properties_window();

  if (algo_pane_visible_) {
    ImGui::SetNextWindowPos(ImVec2(0.0f, float(MENU_HEIGHT() + h + 1)),
                            ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(float(PANE_WIDTH()), float(h - 1)),
                             ImGuiCond_Always);

    ImGui::Begin("Contours", fixed_layout_ ? &algo_pane_visible_ : nullptr,
                 fixed_layout_
                     ? (ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                        ImGuiWindowFlags_NoCollapse)
                     : 0);

    ImGui::Text("Camera");
    ImGui::Combo(" ", (int *)&m_CameraMode, "main\0side\0\0");

    ImGui::Separator();
    ImGui::Text("Algorithm");
    ImGui::Combo("", (int *)&m_ContourMode, "mesh\0interpolated\0\0");

    ImGui::Checkbox("surf. intersect.", &m_surface_intersections);

    if (m_compute_chains) {
      ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
      ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    }
    ImGui::Checkbox("edge visibility", &m_compute_visibility);
    if (m_compute_chains) {
      ImGui::PopItemFlag();
      ImGui::PopStyleVar();
    }

    if (m_compute_visibility) {
      ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
      ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    }
    ImGui::Checkbox("chaining", &m_compute_chains);
    if (m_compute_visibility) {
      ImGui::PopItemFlag();
      ImGui::PopStyleVar();
    }

    if (!m_compute_chains) {
      ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
      ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    }
    ImGui::Checkbox("view graph", &m_build_view_graph);
    if (!m_compute_chains) {
      ImGui::PopItemFlag();
      ImGui::PopStyleVar();
    }

    ImGui::End();
  }
}

void Viewer::draw_windows_menu() {
  if (ImGui::BeginMenu("Windows")) {
    ImGui::MenuItem("object properties", nullptr, &right_pane_visible_);
    ImGui::MenuItem("viewer properties", nullptr, &left_pane_visible_);
    ImGui::MenuItem("contour options", nullptr, &algo_pane_visible_);
    ImGui::Separator();
    if (ImGui::MenuItem(
            "show/hide GUI [T]", nullptr,
            (bool *)glup_viewer_is_enabled_ptr(GLUP_VIEWER_TWEAKBARS))) {
      glup_viewer_post_redisplay();
    }
    ImGui::Separator();
    if (ImGui::MenuItem("Big text", nullptr, &retina_mode_)) {
      scaling_ = retina_mode_ ? 2.0f : 1.0f;
      glup_viewer_post_redisplay();
    }
    ImGui::EndMenu();
  }
}

void Viewer::draw_object_properties() {
  SimpleMeshApplication::draw_object_properties();

  ImGui::Separator();
  ImGui::Text("Display");
  ImGui::Checkbox("normals [n]", &m_normals);
  ImGui::Checkbox("contours", &m_contours);
  ImGui::Checkbox("boundaries", &m_boundaries);
  ImGui::Checkbox("surf. intersection", &m_surfintersect);
  ImGui::Checkbox("visible only [v]", &m_visible_only);
  ImGui::Separator();
  ImGui::Checkbox("points [d]", &m_debug_points);
  if (m_debug_points) {
    ImGui::Checkbox("3D intersections", &m_3D_intersections);
    ImGui::Checkbox("2D intersections", &m_2D_intersections);
    ImGui::Checkbox("curtain folds", &m_curtain_folds);
    ImGui::Checkbox("ray intersections", &m_ray_intersections);
  }

  ImGui::Separator();
#ifdef GEO_OS_EMSCRIPTEN
  ImGui::SliderFloat("line w.", &m_linewidth, 0.1f, 10.0f, "%.1f");
#endif
  ImGui::SliderFloat("pt. sz.", &m_pointSize, 1.f, 20.0f, "%.1f");
  if (m_normals)
    ImGui::SliderFloat("n. wid", &m_normal_length, 0.001f, 1.f, "%.3f");
}

void Viewer::draw_about() {
  ImGui::Separator();
  if (ImGui::BeginMenu("About...")) {
    ImGui::Text("This viewer has been developed in support of the tutorial:");
    ImGui::Text("\"Line Drawings from 3D Models\"");
    ImGui::Text("by Pierre Bénard and Aaron Hertzmann");
    ImGui::Text("https://arxiv.org/abs/1810.01175");
    ImGui::Separator();
    ImGui::Text("The code is based on GEOGRAM:");
    ImGui::Text("http://alice.loria.fr/software/geogram");
    ImGui::EndMenu();
  }
}

void Viewer::draw_scene() {
  if (mesh_gfx_.mesh() == nullptr) {
    return;
  }

  // Switch cameras
  if (m_CameraMode != m_prev_CameraMode) {
    glup_viewer_set_scene_translation(m_trans[m_CameraMode].data());
    glup_viewer_set_scene_quaternion(m_rot[m_CameraMode].data());
    glup_viewer_set_float(GLUP_VIEWER_ZOOM, m_zoom[m_CameraMode]);
    m_prev_CameraMode = m_CameraMode;
    draw_scene();
    return;
  } else {
    m_trans[m_CameraMode] = Eigen::Map<Eigen::Matrix<float, 3, 1>>(
        glup_viewer_get_scene_translation());
    m_rot[m_CameraMode] = Eigen::Map<Eigen::Matrix<float, 4, 1>>(
        glup_viewer_get_scene_quaternion());
    m_zoom[m_CameraMode] = glup_viewer_get_float(GLUP_VIEWER_ZOOM);
  }

  SimpleMeshApplication::draw_scene();

  if(m_backup_mesh == nullptr)
    return;

  // Update contours
  if (m_mesh == nullptr || m_ContourMode != m_prev_ContourMode) {
    if (m_CameraMode == 0) {
      m_modelview =
          Eigen::Map<Matrix4f>(glup_viewer_get_saved_modelview_matrix());
      m_projection =
          Eigen::Map<Matrix4f>(glup_viewer_get_saved_projection_matrix());
    }
    extractContours(m_modelview.inverse().topRightCorner(3, 1));
    m_prev_ContourMode = m_ContourMode;
  } else if (m_CameraMode == 0) {
    Eigen::Map<Matrix4f> modelview =
        Eigen::Map<Matrix4f>(glup_viewer_get_saved_modelview_matrix());
    Eigen::Map<Matrix4f> projection =
        Eigen::Map<Matrix4f>(glup_viewer_get_saved_projection_matrix());
    if (!modelview.isApprox(m_modelview) ||
        !projection.isApprox(m_projection)) {
      extractContours(modelview.inverse().topRightCorner(3, 1));
      m_modelview = modelview;
      m_projection = projection;
    }
  }

  glLineWidth(m_linewidth);
  glupEnable(GLUP_VERTEX_COLORS);
  glupDisable(GLUP_LIGHTING);

  // Draw its normals
  if (m_normals) {
    MapXf positions = m_mesh->get_positions();
    MapXf normals = m_mesh->get_normals();
    glupBegin(GLUP_LINES);
    for (size_t i = 0; i < m_mesh->n_vertices(); ++i) {
      glupColor3f(normalsColor.x(), normalsColor.y(), normalsColor.z());

      glupVertex4f(positions.col(i).x(), positions.col(i).y(),
                   positions.col(i).z(), 1.0f);

      Vector3f p = positions.col(i) + m_normal_length * normals.col(i);
      glupVertex4f(p.x(), p.y(), p.z(), 1.0f);
    }
    glupEnd();
  }

  auto drawCurve = [&](const Matrix3Xf &curve, int start, int end, int step,
                       const Vector4f &color) {
    glupBegin(GLUP_LINES);
    glupColor3f(color.x(), color.y(), color.z());
    for (int j = start; j < end; j += step) {
      glupVertex4f(curve.col(j).x(), curve.col(j).y(), curve.col(j).z(), 1.0f);

      glupVertex4f(curve.col(j + 1).x(), curve.col(j + 1).y(),
                   curve.col(j + 1).z(), 1.0f);
    }
    glupEnd();
  };

  // Draw contours
  if (m_contours) {
    Matrix3Xf contours;
    if (!m_compute_chains) {
      if (m_compute_visibility) { // EDGES
        contours = m_mesh->get_contours(VISIBLE);
        drawCurve(contours, 0, contours.cols() - 1, 2, contourColor);
        if (!m_visible_only) {
          contours = m_mesh->get_contours(INVISIBLE);
          drawCurve(contours, 0, contours.cols() - 1, 2, hiddenContourColor);
        }
      } else {
        contours = m_mesh->get_contours(VISIBLE_AND_INVISIBLE);
        drawCurve(contours, 0, contours.cols() - 1, 2, contourColor);
      }
    } else { // CHAINS
      contours = m_mesh->get_chains(ContourMode(m_ContourMode));
      int offset = (ContourMode(m_ContourMode) == INTERPOLATED_CONTOUR ? 0 : 1);
      int start_index = 0;
      for (size_t i = 0; i < m_mesh->get_chain_lengths().size(); i++) {
        drawCurve(contours, start_index,
                  start_index + m_mesh->get_chain_lengths()[i] + offset - 1, 1,
                  alphabetColors[i % 26]);
        start_index += m_mesh->get_chain_lengths()[i] + offset;
      }
    }
  }

  // Draw open boundaries
  if (m_boundaries) {
    Matrix3Xf boundaries = m_mesh->get_boundaries(false);
    int start_index = 0;
    for (size_t i = 0; i < m_mesh->get_boundaries_lengths().size(); i++) {
      drawCurve(boundaries, start_index,
                start_index + m_mesh->get_boundaries_lengths()[i] - 1, 1,
                boundariesColor);
      start_index += m_mesh->get_boundaries_lengths()[i];
    }
  }

  // Draw surface-surface intersections
  if (m_surfintersect) {
    Matrix3Xf &surf_intersect = m_mesh->get_surface_intersections();
    int start_index = 0;
    for (size_t i = 0; i < m_mesh->get_surface_intersections_lengths().size();
         i++) {
      drawCurve(surf_intersect, start_index,
                start_index + m_mesh->get_surface_intersections_lengths()[i] -
                    1,
                1, surfIntersectionsColor);
      start_index += m_mesh->get_surface_intersections_lengths()[i];
    }
  }

  // Draw camera in side-view mode
  if (m_CameraMode == 1) {
    Matrix4f inverse_modelview = m_modelview.inverse();
    int mVpWidth, mVpHeight;
    glup_viewer_get_screen_size(&mVpWidth, &mVpHeight);
    float ym = tan(25.f * (float)M_PI / 360.0f);
    float xm = float(mVpWidth) * (ym * 1.0 / float(mVpHeight));
    float zm = 0.25f;

    glupPushMatrix();
    glupMultMatrixd(inverse_modelview.data());

    glupBegin(GLUP_LINES);
    glupColor3f(wireframeColor.x(), wireframeColor.y(), wireframeColor.z());

    glupVertex4f(0.f, 0.f, 0.f, 1.f);
    glupVertex4f(xm, ym, -zm, 1.f);
    glupVertex4f(xm, ym, -zm, 1.f);
    glupVertex4f(xm, -ym, -zm, 1.f);
    glupVertex4f(xm, -ym, -zm, 1.f);
    glupVertex4f(0.f, 0.f, 0.f, 1.f);
    glupVertex4f(0.f, 0.f, 0.f, 1.f);
    glupVertex4f(-xm, ym, -zm, 1.f);
    glupVertex4f(-xm, ym, -zm, 1.f);
    glupVertex4f(-xm, -ym, -zm, 1.f);
    glupVertex4f(-xm, -ym, -zm, 1.f);
    glupVertex4f(0.f, 0.f, 0.f, 1.f);
    glupVertex4f(0.f, 0.f, 0.f, 1.f);
    glupVertex4f(xm, ym, -zm, 1.f);
    glupVertex4f(xm, ym, -zm, 1.f);
    glupVertex4f(-xm, ym, -zm, 1.f);
    glupVertex4f(-xm, ym, -zm, 1.f);
    glupVertex4f(0.f, 0.f, 0.f, 1.f);
    glupVertex4f(0.f, 0.f, 0.f, 1.f);
    glupVertex4f(xm, -ym, -zm, 1.f);
    glupVertex4f(xm, -ym, -zm, 1.f);
    glupVertex4f(-xm, -ym, -zm, 1.f);
    glupVertex4f(-xm, -ym, -zm, 1.f);
    glupVertex4f(0.f, 0.f, 0.f, 1.f);
    glupVertex4f(-xm, ym, -zm, 1.f);
    glupEnd();

    glupBegin(GLUP_TRIANGLES);
    glupColor3f(wireframeColor.x(), wireframeColor.y(), wireframeColor.z());
    glupVertex4f(-0.65 * xm, 1.1f * ym, -zm, 1.f);
    glupVertex4f(0.f, 1.5f * ym, -zm, 1.f);
    glupVertex4f(0.65 * xm, 1.1f * ym, -zm, 1.f);
    glupEnd();

    glupPopMatrix();
  }

  // Draw debug points
  if (m_debug_points && (m_3D_intersections || m_2D_intersections ||
                         m_curtain_folds || m_ray_intersections)) {

    nature_t nature = VertexNature::S_VERTEX;
    if (m_build_view_graph) {
      if (m_3D_intersections)
        nature = VertexNature::INTERSECTION_3D | nature;
      if (m_2D_intersections)
        nature = VertexNature::INTERSECTION_2D | nature;
      if (m_curtain_folds)
        nature = VertexNature::CURTAIN_FOLD | nature;
    }
    if (m_ray_intersections)
      nature = VertexNature::RAY_INTERSECTION | nature;

    std::vector<Vector3f> &debug_points = m_mesh->get_debug_points(nature),
                          &debug_points_color = m_mesh->get_point_colors();

    glupSetPointSize(m_pointSize);
    glupBegin(GLUP_POINTS);
    for (size_t i = 0; i < debug_points.size(); ++i) {
      glupColor3f(debug_points_color[i].x(), debug_points_color[i].y(),
                  debug_points_color[i].z());

      glupVertex3f(debug_points[i].x(), debug_points[i].y(),
                   debug_points[i].z());
    }
    glupEnd();

    // 2D intersection lines
    if (m_mesh->num2Dintersections().first >= 0) {
      debug_points = m_mesh->get_debug_points(VertexNature::INTERSECTION_2D);
      glupBegin(GLUP_LINES);
      for (size_t i = 0; i < debug_points.size(); i += 2) {
        glupColor3f(intersection2DColor.x(), intersection2DColor.y(),
                    intersection2DColor.z());
        glupVertex3f(debug_points[i].x(), debug_points[i].y(),
                     debug_points[i].z());
        glupVertex3f(debug_points[i + 1].x(), debug_points[i + 1].y(),
                     debug_points[i + 1].z());
      }
      glupEnd();
    }
  }

  glupEnable(GLUP_LIGHTING);
}

void Viewer::extractContours(const Vector3f &view_pos) {
  if (m_backup_mesh == nullptr)
    return;

  try {
    ContourMode mode = ContourMode(m_ContourMode);

    bool extract_boundaries = true;
    bool extract_surfaceIntersections = m_surface_intersections;
    if (m_mesh) {
      extract_boundaries = m_mesh->hasBoundaries();
      extract_surfaceIntersections = m_mesh->hasSurfaceIntersections();
      delete m_mesh;
    }
    m_mesh = new Mesh(*m_backup_mesh);

    if (extract_boundaries) {
      m_mesh->extractBoundaries();
      m_mesh->extractBoundaryCurtainFolds(view_pos);
    }

    if (extract_surfaceIntersections) {
      m_mesh->extractSurfaceIntersections();
    }

    m_mesh->extractContours(mode, view_pos);

    if (m_compute_chains) {
      m_mesh->computeContourChains(mode);
    }

    if (m_build_view_graph) {
      Vector2i mSize;
      glup_viewer_get_screen_size(&mSize(1), &mSize(0));
      m_mesh->projectSVertices(m_modelview.matrix(), m_projection, mSize);
      m_mesh->computeSweepLineIntersections(mode, m_projection, m_modelview,
                                            mSize);
    }
    if (m_compute_visibility) {
      m_mesh->computeContourVisibility(view_pos, mode);
    }
  } catch (const std::runtime_error &e) {
    std::string error_msg =
        std::string("Caught a fatal error: ") + std::string(e.what());
    std::cerr << error_msg << endl;
    std::abort();
  }
}

Viewer::~Viewer() {
  if (m_mesh)
    delete m_mesh;
  if (m_backup_mesh)
    delete m_backup_mesh;
}
