#include "svg.h"

using namespace std;

SVG::SVG(const std::string &filename, const Vector2i &size)
    : m_size(size) {

  m_file.open(filename, ofstream::out);
  if (!m_file.is_open()) {
    std::cout << "ERROR: Unable to open SVG output file " << filename
              << std::endl;
  }

  header();
}

SVG::~SVG() {
  footer();
  m_file.close();
}

void SVG::header() {
  m_file << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>"
         << endl
         << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl
         << " \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
         << "<svg viewBox=\"0 0 " << m_size.x() << " " << m_size.y()
         << "\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">" << endl;
}

void SVG::footer() { m_file << "</svg>" << endl; }

void SVG::openAnimGroup() { m_file << "<g visibility=\"hidden\">"; }

void SVG::closeAnimGroup(int begin, int end, real_t framerate) {
  m_file << "<animate id=\"" << begin
         << "\" attributeName=\"visibility\" attributeType=\"XML\" "
            "begin=\"";
  if (begin > 0) {
    m_file << (begin - 1) << ".end\"";
  } else {
    m_file << "0s;" << end << ".end\"";
  }
  m_file << " dur=\"" << framerate << "s\" to=\"visible\"/>" << endl;
  m_file << "</g>" << endl;
}

void SVG::writePolyline(const std::vector<Vector2f> &polyline,
                        int stroke_width, const Color &color, bool segments) {

  m_file << "<g fill=\"none\" stroke=\"rgb(" << int(color.x() * 255) << ","
         << int(color.y() * 255) << "," << int(color.z() * 255)
         << ")\" stroke-width=\"" << stroke_width << "\">";

  m_file << "<path d=\"";
  bool first = true;
  for (size_t i = 0; i < polyline.size(); ++i) {
    const Vector2f &p = polyline.at(i);
    m_file << (first ? "M" : "L") << p.x() << " " << m_size.y() - p.y() << " ";
    if (segments) {
      i++;
      const Vector2f &p = polyline.at(i);
      m_file << "L" << p.x() << " " << m_size.y() - p.y();
      if (i < polyline.size() - 1) {
        m_file << "\"/>" << endl << "<path d=\"";
      }
    } else {
      first = false;
    }
  }
  m_file << "\"/>" << endl;
  m_file << "</g>" << endl;
}