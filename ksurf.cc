#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>

//#include <surface-generalized-coons.hh>
#include <surface-midpoint.hh>

#include "kcurve.hh"
#include "ksurf.hh"

KSurf::KSurf(std::string filename) : Object(filename) {
    reload();
}

KSurf::~KSurf() {
}

void KSurf::draw(const Visualization &vis) const {
    Object::draw(vis);
    if (vis.show_control_points) {
        glDisable(GL_LIGHTING);
        glLineWidth(3.0);
        glColor3d(0.3, 0.3, 1.0);
        for (auto e : cage.edges()) {
            glBegin(GL_LINE_STRIP);
            glVertex3dv(cage.point(e.v0()).data());
            glVertex3dv(cage.point(e.v1()).data());
            glEnd();
        }
        glLineWidth(1.0);
        glPointSize(8.0);
        glColor3d(1.0, 0.0, 1.0);
        glBegin(GL_POINTS);
        for (auto v : cage.vertices())
            glVertex3dv(cage.point(v).data());
        glEnd();
        glPointSize(1.0);
        glEnable(GL_LIGHTING);
    }
}

void KSurf::drawWithNames(const Visualization &vis) const {
    if (!vis.show_control_points)
        return;
    for (auto v : cage.vertices()) {
        glPushName(v.idx());
        glRasterPos3dv(cage.point(v).data());
        glPopName();
    }
}

Vector KSurf::postSelection(int selected) {
    return cage.point(Cage::VertexHandle(selected));
}

void KSurf::movement(int selected, const Vector &pos) {
    cage.set_point(Cage::VertexHandle(selected), pos);
}

static std::shared_ptr<Transfinite::Curve> createCurve(
        const Vector &p1, const Vector &n1,
        const Vector &p2, const Vector &n2) {
#if 0
    // Simple cubic Bezier version
    auto d1 = p1 + (p2 - p1) / 3, d2 = p2 + (p1 - p2) / 3;
    d1 -= n1 * ((d1 - p1) | n1); d2 -= n2 * ((d2 - p2) | n2);
    Transfinite::Point3D cp1(p1.data()), cp2(d1.data()), cp3(d2.data()), cp4(p2.data());
    Transfinite::BSplineCurve curve(Transfinite::BSCurve({ cp1, cp2, cp3, cp4 }));
    return std::make_shared<Transfinite::BSplineCurve>(curve);
#else
    // Kappa-curve version
    return std::make_shared<KappaCurve>(
            Transfinite::Point3D(p1.data()),
            Transfinite::Vector3D(n1.data()),
            Transfinite::Point3D(p2.data()),
            Transfinite::Vector3D(n2.data()));
#endif
}

//#define CURVE_OUTPUT
#ifdef CURVE_OUTPUT
Transfinite::CurveVector all_curves;
#endif

std::unique_ptr<Transfinite::Surface> KSurf::createPatch(Cage::FaceHandle f) const {
    Transfinite::CurveVector curves;
    for (auto he : cage.fh_range(f)) {
        auto v1 = he.from(), v2 = he.to();
        auto p1 = cage.point(v1), p2 = cage.point(v2);
        auto n1 = cage.normal(v1), n2 = cage.normal(v2);
        curves.push_back(createCurve(p1, n1, p2, n2));
    } 
#ifdef CURVE_OUTPUT
    all_curves.insert(all_curves.end(), curves.begin(), curves.end());
#endif
    std::unique_ptr<Transfinite::Surface> patch =
        //std::make_unique<Transfinite::SurfaceGeneralizedCoons>();
        std::make_unique<Transfinite::SurfaceMidpoint>();
    patch->setCurves(curves);
    patch->setupLoop();
    patch->update();
    return patch;
}

void KSurf::updateBaseMesh() {
    size_t resolution = 15;
    double tolerance = 1e-5;
    Transfinite::TriMesh m;
    mesh.clear();
#ifdef CURVE_OUTPUT
    all_curves.clear();
#endif
    for (auto f : cage.faces()) {
        auto patch = createPatch(f);
        m.insert(patch->eval(resolution), tolerance);
    }
#ifdef CURVE_OUTPUT
    std::ofstream f("/tmp/curves.obj");
    for (size_t i = 0; i < all_curves.size(); ++i) {
        for (size_t j = 0; j < 100; ++j) {
            double u = j / 99.0;
            f << "v " << all_curves[i]->eval(u) << std::endl;
        }
        f << 'l';
        for (size_t j = 0; j < 100; ++j)
            f << ' ' << i * 100 + j + 1;
        f << std::endl;
    }
#endif
    std::vector<BaseMesh::VertexHandle> handles;
    for (const auto &p : m.points())
        handles.push_back(mesh.add_vertex({p[0], p[1], p[2]}));
    for (const auto &t : m.triangles())
        mesh.add_face({handles[t[0]], handles[t[1]], handles[t[2]]}); 
    Object::updateBaseMesh(false, false);
}

bool KSurf::reload() {
    if (!OpenMesh::IO::read_mesh(cage, filename))
        return false;
    cage.request_face_normals();
    cage.request_vertex_normals();
    cage.update_face_normals();
    cage.update_vertex_normals();
    updateBaseMesh();
    return true;
}
