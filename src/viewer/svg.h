#ifndef SVG_H_
#define SVG_H_

#include "common.h"
#include <fstream>

class SVG {
public:
  SVG(const std::string &filename, const Vector2i &size);
  ~SVG();

  void openAnimGroup();
  void closeAnimGroup(int begin, int end, real_t framerate = 1.f / 12.f);

  void writePolyline(const std::vector<Vector2f> &polyline,
                     int stroke_width, const Color &color, bool segments);

private:
  void header();
  void footer();

  Vector2i m_size;
  std::ofstream m_file;
};

#endif
