#ifndef TGEOTET_H
#define TGEOTET_H

#include "TGeoBBox.h"
#include "TGeoVolume.h"
#include "TVector3.h"

class TGeoTet : public TGeoBBox
{
  public:
    TGeoTet(); // Do not use this. Can't = delete because of dictionary generator
    TGeoTet(TVector3 anchor, TVector3 p2, TVector3 p3, TVector3 p4);

    Double_t Capacity() const override { return fCubicVolume; }

    void ComputeBBox() override;

    void ComputeNormal(const Double_t* point, const Double_t* dir, Double_t* norm) override;

    Bool_t Contains(const Double_t* point) const override;

    // virtual Bool_t CouldBeCrossed(const Double_t* point, const Double_t* dir) const; // Use BBox
    // Impl

    Int_t DistancetoPrimitive(Int_t px, Int_t py) override;

    Double_t DistFromInside(const Double_t* point,
                            const Double_t* dir,
                            Int_t iact,
                            Double_t step,
                            Double_t* safe) const override;

    Double_t DistFromOutside(const Double_t* point,
                             const Double_t* dir,
                             Int_t iact,
                             Double_t step,
                             Double_t* safe) const override;

    TGeoVolume* Divide(TGeoVolume* voldiv, const char* divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step)
        override;

    // virtual Double_t GetAxisRange(Int_t iaxis, Double_t& xlo, Double_t& xhi) const; // Use BBox
    // Impl

    void GetBoundingCylinder(Double_t* param) const override;

    const TBuffer3D& GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;

    Int_t GetByteCount() const override { return sizeof(*this); } // Note: huh?

    Double_t GetFacetArea(Int_t index) const override;

    Bool_t GetPointsOnFacet(Int_t index, Int_t npoints, Double_t* array) const override;

    Bool_t GetPointsOnSegments(Int_t npoints, Double_t* array) const override;

    Int_t GetFittingBox(const TGeoBBox* parambox,
                        TGeoMatrix* mat,
                        Double_t& dx,
                        Double_t& dy,
                        Double_t& dz) const override;

    TGeoShape* GetMakeRuntimeShape(TGeoShape* mother, TGeoMatrix* mat) const override;

    void GetMeshNumbers(Int_t& nvert, Int_t& nsegs, Int_t& npols) const override;

    Int_t GetNmeshVertices() const override { return 4; }

    void InspectShape() const override;

    TBuffer3D* MakeBuffer3D() const override;

    Double_t Safety(const Double_t* point, Bool_t in) const override;

    void SavePrimitive(std::ostream& out, Option_t* option) override;

    void SetDimensions(Double_t* param) override;

    void SetPoints(Double_t* points) const override;

    void SetPoints(Float_t* points) const override;

    void SetSegsAndPols(TBuffer3D& buffer) const override;

    // From G4Tet
    TVector3 SurfaceNormal(const TVector3& p) const;
    Bool_t Inside(const TVector3& p) const;
    Double_t DistanceToIn(const TVector3& p, const TVector3& v) const;
    Double_t DistanceToIn(const TVector3& p) const;
    Double_t DistanceToOut(const TVector3& p,
                           const TVector3& v,
                           const Bool_t calcNorm = false,
                           Bool_t* validNorm = nullptr,
                           TVector3* n = nullptr) const;
    Double_t DistanceToOut(const TVector3& p) const;

  protected:
    const Double_t kInfinity;
    Double_t kCarTolerance;

    TVector3 fAnchor, fP2, fP3, fP4, fMiddle;

    TVector3 fNormal123, fNormal142, fNormal134, fNormal234;
    Double_t fCdotN123, fCdotN142, fCdotN134, fCdotN234;
    Double_t fXMin, fXMax, fYMin, fYMax, fZMin, fZMax;

    Double_t fCubicVolume, fTol, fMaxSize;

    ClassDefOverride(TGeoTet, 1)
};

#endif // TGEOTET_H
