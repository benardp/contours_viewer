#ifndef VIEWER_H
#define VIEWER_H

#include "bvh.h"
#include "mesh.h"

#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>

class Viewer : public GEO::SimpleMeshApplication {
public:
  Viewer(int argc, char **argv, const std::string &usage);
  ~Viewer();

  virtual void init_graphics();
  virtual void draw_scene();

  virtual void draw_windows_menu();
  virtual void draw_left_pane();
  virtual void draw_object_properties();
  virtual void draw_about();

  virtual bool load(const std::string &filename);
  virtual bool save(const std::string &filename);
  virtual std::string supported_write_file_extensions();

  void loadMesh(const std::string &filename);

  void saveSVG(const std::string &filename, ContourMode mode,
               bool visible_only);

private:
  void extractContours(const Vector3f &view_pos);
  void setStyle();

  bool m_normals, m_surface, m_contours, m_draw_camera, m_2D_intersections,
      m_3D_intersections, m_curtain_folds, m_compute_visibility,
      m_build_view_graph, m_visible_only, m_boundaries, m_surfintersect,
      m_surface_intersections, m_debug_points, algo_pane_visible_;

  int m_ContourMode, m_prev_ContourMode, m_CameraMode;

  float m_linewidth, m_pointSize, m_normal_length;

  Mesh *m_mesh, *m_backup_mesh;

  Matrix4f m_modelview, m_projection;
};

#endif // _VIEWER_H
