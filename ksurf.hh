#pragma once

#include "object.hh"

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

namespace Transfinite { class Curve; class Surface; }

class KSurf : public Object {
public:
    KSurf(std::string filename);
    virtual ~KSurf();
    virtual void draw(const Visualization &vis) const override;
    virtual void drawWithNames(const Visualization &vis) const override;
    virtual Vector postSelection(int selected) override;
    virtual void movement(int selected, const Vector &pos) override;
    virtual void updateBaseMesh() override;
    virtual bool reload() override;
private:
    struct CageTraits : public OpenMesh::DefaultTraits {
        using Point  = OpenMesh::Vec3d;
        using Normal = OpenMesh::Vec3d;
    };
    using Cage = OpenMesh::PolyMesh_ArrayKernelT<CageTraits>;
    std::unique_ptr<Transfinite::Surface> createPatch(Cage::FaceHandle f) const;
    std::map<Cage::EdgeHandle, std::shared_ptr<Transfinite::Curve>> boundaries;
    Cage cage;
};

