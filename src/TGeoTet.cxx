#include "TGeoTet.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

TGeoVolume* TGeoTet::Divide(TGeoVolume* /*voldiv*/,
                            const char* /*divname*/,
                            Int_t /*iaxis*/,
                            Int_t /*ndiv*/,
                            Double_t /*start*/,
                            Double_t /*step*/)
{
    throw std::runtime_error("TGeoTet::Divide Not yet implemented");
    return nullptr;
}

void TGeoTet::GetBoundingCylinder(Double_t* /*param*/) const { throw std::runtime_error("Not yet implemented"); }

Double_t TGeoTet::GetFacetArea(Int_t /*index*/) const
{
    throw std::runtime_error("TGeoTet::GetBoundingCylinder Not yet implemented");
    return 0.;
}

Bool_t TGeoTet::GetPointsOnFacet(Int_t /*index*/, Int_t /*npoints*/, Double_t* /*array*/) const
{
    throw std::runtime_error("TGeoTet::GetPointsOnFacet Not yet implemented");
    return false;
}

Int_t TGeoTet::GetFittingBox(const TGeoBBox* /*parambox*/,
                             TGeoMatrix* /*mat*/,
                             Double_t& /*dx*/,
                             Double_t& /*dy*/,
                             Double_t& /*dz*/) const
{
    throw std::runtime_error("TGeoTet::GetFittingBox Not yet implemented");
    return 0;
}

void TGeoTet::SetDimensions(Double_t* /*param*/)
{
    throw std::runtime_error("TGeoTet::SetDimensions Not yet implemented");
}

TGeoTet::TGeoTet()
    : TGeoBBox()
    , kInfinity(TGeoShape::Big())
{
    throw std::runtime_error("TGeoTet::TGeoTet() Not yet implemented");
}

// From TGeoTet constructor
TGeoTet::TGeoTet(TVector3 anchor, TVector3 p2, TVector3 p3, TVector3 p4)
    : TGeoBBox()
    , kInfinity(TGeoShape::Big())
{
    // fV<x><y> is vector from vertex <y> to vertex <x>
    //
    TVector3 fV21 = p2 - anchor;
    TVector3 fV31 = p3 - anchor;
    TVector3 fV41 = p4 - anchor;

    // make sure this is a correctly oriented set of points for the tetrahedron
    //
    Double_t signed_vol = fV21.Cross(fV31).Dot(fV41);
    if (signed_vol < 0.0)
    {
        std::swap(p3, p4);
        std::swap(fV31, fV41);
    }
    fCubicVolume = std::fabs(signed_vol) / 6.;

    TVector3 fV24 = p2 - p4;
    TVector3 fV43 = p4 - p3;
    TVector3 fV32 = p3 - p2;

    fXMin = std::min(std::min(std::min(anchor.x(), p2.x()), p3.x()), p4.x());
    fXMax = std::max(std::max(std::max(anchor.x(), p2.x()), p3.x()), p4.x());
    fYMin = std::min(std::min(std::min(anchor.y(), p2.y()), p3.y()), p4.y());
    fYMax = std::max(std::max(std::max(anchor.y(), p2.y()), p3.y()), p4.y());
    fZMin = std::min(std::min(std::min(anchor.z(), p2.z()), p3.z()), p4.z());
    fZMax = std::max(std::max(std::max(anchor.z(), p2.z()), p3.z()), p4.z());

    fDX = (fXMax - fXMin) * 0.5;
    fDY = (fYMax - fYMin) * 0.5;
    fDZ = (fZMax - fZMin) * 0.5;

    fMiddle = TVector3(fXMax + fXMin, fYMax + fYMin, fZMax + fZMin) * 0.5;
    fMaxSize = std::max(std::max(std::max((anchor - fMiddle).Mag(), (p2 - fMiddle).Mag()), (p3 - fMiddle).Mag()),
                        (p4 - fMiddle).Mag());

    fTol = 1e-9 * (std::fabs(fXMin) + std::fabs(fXMax) + std::fabs(fYMin) + std::fabs(fYMax) + std::fabs(fZMin) +
                   std::fabs(fZMax));

    fAnchor = anchor;
    fP2 = p2;
    fP3 = p3;
    fP4 = p4;

    TVector3 fCenter123 = (anchor + p2 + p3) * (1.0 / 3.0); // face center
    TVector3 fCenter134 = (anchor + p4 + p3) * (1.0 / 3.0);
    TVector3 fCenter142 = (anchor + p4 + p2) * (1.0 / 3.0);
    TVector3 fCenter234 = (p2 + p3 + p4) * (1.0 / 3.0);

    // compute area of each triangular face by cross product
    // and sum for total surface area

    TVector3 normal123 = fV31.Cross(fV21);
    TVector3 normal134 = fV41.Cross(fV31);
    TVector3 normal142 = fV21.Cross(fV41);
    TVector3 normal234 = fV32.Cross(fV43);

    fNormal123 = normal123.Unit();
    fNormal134 = normal134.Unit();
    fNormal142 = normal142.Unit();
    fNormal234 = normal234.Unit();

    fCdotN123 = fCenter123.Dot(fNormal123);
    fCdotN134 = fCenter134.Dot(fNormal134);
    fCdotN142 = fCenter142.Dot(fNormal142);
    fCdotN234 = fCenter234.Dot(fNormal234);
}

void TGeoTet::ComputeBBox()
{
    // already calculated
}

// Did someone say outdated legacy interface?
void TGeoTet::ComputeNormal(const Double_t* point, const Double_t* dir, Double_t* norm)
{
    const auto n = SurfaceNormal({ point[0], point[1], point[2] });

    // TODO: Not sure if this is corrected and the intended behavior
    if (n.X() * dir[0] + n.Y() * dir[1] + n.Z() * dir[2] < 0)
    {
        norm[0] = -n.X();
        norm[1] = -n.Y();
        norm[2] = -n.Z();
    }
    else
    {
        norm[0] = +n.X();
        norm[1] = +n.Y();
        norm[2] = +n.Z();
    }
}

// From TGeoTet::SurfaceNormal
TVector3 TGeoTet::SurfaceNormal(const TVector3& p) const
{
    Double_t r123 = std::fabs(p.Dot(fNormal123) - fCdotN123);
    Double_t r134 = std::fabs(p.Dot(fNormal134) - fCdotN134);
    Double_t r142 = std::fabs(p.Dot(fNormal142) - fCdotN142);
    Double_t r234 = std::fabs(p.Dot(fNormal234) - fCdotN234);

    const Double_t delta = 0.5 * kCarTolerance;
    TVector3 sumnorm(0., 0., 0.);
    Int_t noSurfaces = 0;

    if (r123 <= delta)
    {
        noSurfaces++;
        sumnorm = fNormal123;
    }

    if (r134 <= delta)
    {
        noSurfaces++;
        sumnorm += fNormal134;
    }

    if (r142 <= delta)
    {
        noSurfaces++;
        sumnorm += fNormal142;
    }
    if (r234 <= delta)
    {
        noSurfaces++;
        sumnorm += fNormal234;
    }

    if (noSurfaces > 0)
    {
        if (noSurfaces == 1)
        {
            return sumnorm;
        }
        else
        {
            return sumnorm.Unit();
        }
    }
    else // Approximative Surface Normal
    {
        if ((r123 <= r134) && (r123 <= r142) && (r123 <= r234))
        {
            return fNormal123;
        }
        else if ((r134 <= r142) && (r134 <= r234))
        {
            return fNormal134;
        }
        else if (r142 <= r234)
        {
            return fNormal142;
        }
        return fNormal234;
    }
}

// this method returns true if POINT is actually inside the shape. The point has to be defined in
// the local shape reference.
Bool_t TGeoTet::Contains(const Double_t* point) const { return Inside(TVector3(point[0], point[1], point[2])); }

// From G4Tet::Inside
Bool_t TGeoTet::Inside(const TVector3& p) const
{
    // Note: Original has a too close to surface check
    return !((p.Dot(fNormal123) - fCdotN123) > fTol || (p.Dot(fNormal134) - fCdotN134) > fTol ||
             (p.Dot(fNormal142) - fCdotN142) > fTol || (p.Dot(fNormal234) - fCdotN234) > fTol);
}

// Calculate distance to box from an outside point
// - return kInfinity if no intersection.
// All this is very unrolled, for speed.
Double_t TGeoTet::DistanceToIn(const TVector3& p, const TVector3& v) const
{
    TVector3 vu(v.Unit()), hp;
    Double_t vdotn, t, tmin = kInfinity;

    Double_t extraDistance = 10.0 * fTol; // a little ways into the solid

    vdotn = -vu.Dot(fNormal123);
    if (vdotn > 1e-12)
    {                                                // this is a candidate face, since it is pointing at us
        t = (p.Dot(fNormal123) - fCdotN123) / vdotn; // #  distance to intersection
        if ((t >= -fTol) && (t < tmin))
        {                                      // if not true, we're going away from this face or it's not close
            hp = p + vu * (t + extraDistance); // a little beyond point of intersection
            if ((hp.Dot(fNormal134) - fCdotN134 < 0.0) && (hp.Dot(fNormal142) - fCdotN142 < 0.0) &&
                (hp.Dot(fNormal234) - fCdotN234 < 0.0))
            {
                tmin = t;
            }
        }
    }

    vdotn = -vu.Dot(fNormal134);
    if (vdotn > 1e-12)
    {                                                // # this is a candidate face, since it is pointing at us
        t = (p.Dot(fNormal134) - fCdotN134) / vdotn; // #  distance to intersection
        if ((t >= -fTol) && (t < tmin))
        {                                      // if not true, we're going away from this face
            hp = p + vu * (t + extraDistance); // a little beyond point of intersection
            if ((hp.Dot(fNormal123) - fCdotN123 < 0.0) && (hp.Dot(fNormal142) - fCdotN142 < 0.0) &&
                (hp.Dot(fNormal234) - fCdotN234 < 0.0))
            {
                tmin = t;
            }
        }
    }

    vdotn = -vu.Dot(fNormal142);
    if (vdotn > 1e-12)
    {                                                // # this is a candidate face, since it is pointing at us
        t = (p.Dot(fNormal142) - fCdotN142) / vdotn; // #  distance to intersection
        if ((t >= -fTol) && (t < tmin))
        {                                      // if not true, we're going away from this face
            hp = p + vu * (t + extraDistance); // a little beyond point of intersection
            if ((hp.Dot(fNormal123) - fCdotN123 < 0.0) && (hp.Dot(fNormal134) - fCdotN134 < 0.0) &&
                (hp.Dot(fNormal234) - fCdotN234 < 0.0))
            {
                tmin = t;
            }
        }
    }

    vdotn = -vu.Dot(fNormal234);
    if (vdotn > 1e-12)
    {                                                // # this is a candidate face, since it is pointing at us
        t = (p.Dot(fNormal234) - fCdotN234) / vdotn; // #  distance to intersection
        if ((t >= -fTol) && (t < tmin))
        {                                      // if not true, we're going away from this face
            hp = p + vu * (t + extraDistance); // a little beyond point of intersection
            if ((hp.Dot(fNormal123) - fCdotN123 < 0.0) && (hp.Dot(fNormal134) - fCdotN134 < 0.0) &&
                (hp.Dot(fNormal142) - fCdotN142 < 0.0))
            {
                tmin = t;
            }
        }
    }

    return std::max(0.0, tmin);
}

// Approximate distance to tet.
// returns distance to sphere centered on bounding box
// - If inside return 0
Double_t TGeoTet::DistanceToIn(const TVector3& p) const
{
    Double_t dd = (p - fMiddle).Mag() - fMaxSize - fTol;
    return std::max(0.0, dd);
}

// Calcluate distance to surface of box from inside
// by calculating distances to box's x/y/z planes.
// Smallest distance is exact distance to exiting.
Double_t TGeoTet::DistanceToOut(const TVector3& p,
                                const TVector3& v,
                                const Bool_t calcNorm,
                                Bool_t* validNorm,
                                TVector3* n) const
{
    TVector3 vu(v.Unit());
    Double_t t1 = kInfinity, t2 = kInfinity, t3 = kInfinity, t4 = kInfinity, vdotn, tt;

    vdotn = vu.Dot(fNormal123);
    if (vdotn > 1e-12) // #we're heading towards this face, so it is a candidate
    {
        t1 = (fCdotN123 - p.Dot(fNormal123)) / vdotn; // #  distance to intersection
    }

    vdotn = vu.Dot(fNormal134);
    if (vdotn > 1e-12) // #we're heading towards this face, so it is a candidate
    {
        t2 = (fCdotN134 - p.Dot(fNormal134)) / vdotn; // #  distance to intersection
    }

    vdotn = vu.Dot(fNormal142);
    if (vdotn > 1e-12) // #we're heading towards this face, so it is a candidate
    {
        t3 = (fCdotN142 - p.Dot(fNormal142)) / vdotn; // #  distance to intersection
    }

    vdotn = vu.Dot(fNormal234);
    if (vdotn > 1e-12) // #we're heading towards this face, so it is a candidate
    {
        t4 = (fCdotN234 - p.Dot(fNormal234)) / vdotn; // #  distance to intersection
    }

    tt = std::min(std::min(std::min(t1, t2), t3), t4);

    if (calcNorm && n)
    {
        TVector3 normal;
        if (tt == t1)
        {
            normal = fNormal123;
        }
        else if (tt == t2)
        {
            normal = fNormal134;
        }
        else if (tt == t3)
        {
            normal = fNormal142;
        }
        else if (tt == t4)
        {
            normal = fNormal234;
        }
        *n = normal;
        if (validNorm)
        {
            *validNorm = true;
        }
    }

    return std::max(tt, 0.0); // avoid tt<0.0 by a tiny bit
    // if we are right on a face
}

// Calculate exact shortest distance to any boundary from inside
// - If outside return 0
Double_t TGeoTet::DistanceToOut(const TVector3& p) const
{
    Double_t t1, t2, t3, t4;
    t1 = fCdotN123 - p.Dot(fNormal123); //  distance to plane, positive if inside
    t2 = fCdotN134 - p.Dot(fNormal134); //  distance to plane
    t3 = fCdotN142 - p.Dot(fNormal142); //  distance to plane
    t4 = fCdotN234 - p.Dot(fNormal234); //  distance to plane

    // if any one of these is negative, we are outside,
    // so return zero in that case

    Double_t tmin = std::min(std::min(std::min(t1, t2), t3), t4);
    return (tmin < fTol) ? 0 : tmin;
}

// computes the distance to exiting a shape from a given point INSIDE, along a given direction. The
// direction is given by its director cosines with respect to the local shape coordinate system.
// This method provides additional information according the value of IACT input parameter. (NOTE:
// W.T.F. This is why we can't have nice things.)
Double_t TGeoTet::DistFromInside(const Double_t* point,
                                 const Double_t* dir,
                                 Int_t iact,
                                 Double_t step,
                                 Double_t* safe) const
{
    TVector3 vPoint(point[0], point[1], point[2]);
    TVector3 vDir(dir[0], dir[1], dir[2]);

    switch (iact)
    {
        case 0:
            // IACT = 0 => compute only safe distance and fill it at the location given by SAFE
            if (safe != nullptr)
            {
                *safe = DistanceToOut(vPoint);
            }
            return kInfinity;
        case 1:
            // IACT = 1 => a proposed STEP is supplied. The safe distance is computed first. If this is
            // bigger than STEP than the proposed step is approved and returned by the method since it
            // does not cross the shape boundaries. Otherwise, the distance to exiting the shape is
            // computed and returned.
            if (safe != nullptr)
            {
                *safe = DistanceToOut(vPoint);
                if (*safe > step)
                {
                    // Note: This does not actually return step, behavior as in TGeoBBox
                    return kInfinity;
                }
            }
            return DistanceToOut(vPoint, vDir);
        case 2:
            // IACT = 2 => compute both safe distance and distance to exiting, ignoring the proposed
            // step.
            if (safe != nullptr)
            {
                *safe = DistanceToOut(vPoint);
            }
            return DistanceToOut(vPoint, vDir);
        default:
            // IACT > 2 => compute only the distance to exiting, ignoring anything else
            return DistanceToOut(vPoint, vDir);
    }
}

// computes the distance to entering a shape from a given point OUTSIDE. Acts in the same way as B).
Double_t TGeoTet::DistFromOutside(const Double_t* point,
                                  const Double_t* dir,
                                  Int_t iact,
                                  Double_t step,
                                  Double_t* safe) const
{
    TVector3 vPoint(point[0], point[1], point[2]);
    TVector3 vDir(dir[0], dir[1], dir[2]);

    switch (iact)
    {
        case 0:
            // IACT = 0 => compute only safe distance and fill it at the location given by SAFE
            if (safe != nullptr)
            {
                *safe = DistanceToIn(vPoint);
            }
            return kInfinity;
        case 1:
            // IACT = 1 => a proposed STEP is supplied. The safe distance is computed first. If this is
            // bigger than STEP than the proposed step is approved and returned by the method since it
            // does not cross the shape boundaries. Otherwise, the distance to exiting the shape is
            // computed and returned.
            if (safe != nullptr)
            {
                *safe = DistanceToIn(vPoint);
                if (*safe > step)
                {
                    // Note: This does not actually return step, behavior as in TGeoBBox
                    return kInfinity;
                }
            }
            return DistanceToIn(vPoint, vDir);
        case 2:
            // IACT = 2 => compute both safe distance and distance to exiting, ignoring the proposed
            // step.
            if (safe != nullptr)
            {
                *safe = DistanceToIn(vPoint);
            }
            return DistanceToIn(vPoint, vDir);
        default:
            // IACT > 2 => compute only the distance to exiting, ignoring anything else
            return DistanceToIn(vPoint, vDir);
    }
}

// compute maximum shift of a point in any direction that does not change its INSIDE/OUTSIDE state
// (does not cross shape boundaries). The state of the point have to be properly supplied.
Double_t TGeoTet::Safety(const Double_t* point, Bool_t in) const
{
    TVector3 vPoint(point[0], point[1], point[2]);
    return in ? DistanceToOut(vPoint) : DistanceToIn(vPoint);
}

Int_t TGeoTet::DistancetoPrimitive(Int_t px, Int_t py)
{
    const Int_t numPoints = 4;
    return ShapeDistancetoPrimitive(numPoints, px, py);
}

// Fills array with n random points located on the line segments of the shape mesh.
// The output array must be provided with a length of minimum 3*npoints. Returns
// true if operation is implemented.
Bool_t TGeoTet::GetPointsOnSegments(Int_t npoints, Double_t* array) const
{
    return false; // Todo
}

// in case shape has some negative parameters, these has to be computed
// in order to fit the mother
TGeoShape* TGeoTet::GetMakeRuntimeShape(TGeoShape* mother, TGeoMatrix* mat) const
{
    // TODO: Not sure what this is supposed to do
    return nullptr;
}

// Returns numbers of vertices, segments and polygons composing the shape mesh.
void TGeoTet::GetMeshNumbers(Int_t& nvert, Int_t& nsegs, Int_t& npols) const
{
    nvert = 4;
    nsegs = 6;
    npols = 4;
}

void TGeoTet::InspectShape() const { throw std::runtime_error("TGeoTet::InspectShape Not yet implemented"); }

TBuffer3D* TGeoTet::MakeBuffer3D() const
{
    TBuffer3D* buff = new TBuffer3D(TBuffer3DTypes::kGeneric, 4, 4 * 3, 6, 6 * 3, 4, 4 * 5);
    if (buff)
    {
        SetPoints(buff->fPnts);
        SetSegsAndPols(*buff);
    }
    return buff;
}

void TGeoTet::SavePrimitive(std::ostream& out, Option_t* option)
{
    throw std::runtime_error("TGeoTet::SavePrimitive Not yet implemented");
}

void TGeoTet::SetPoints(Double_t* points) const
{
    points[0] = fAnchor.X();
    points[1] = fAnchor.Y();
    points[2] = fAnchor.Z();
    points[3] = fP2.X();
    points[4] = fP2.Y();
    points[5] = fP2.Z();
    points[6] = fP3.X();
    points[7] = fP3.Y();
    points[8] = fP3.Z();
    points[9] = fP4.X();
    points[10] = fP4.Y();
    points[11] = fP4.Z();
}

void TGeoTet::SetPoints(Float_t* points) const
{
    points[0] = static_cast<Float_t>(fAnchor.X());
    points[1] = static_cast<Float_t>(fAnchor.Y());
    points[2] = static_cast<Float_t>(fAnchor.Z());
    points[3] = static_cast<Float_t>(fP2.X());
    points[4] = static_cast<Float_t>(fP2.Y());
    points[5] = static_cast<Float_t>(fP2.Z());
    points[6] = static_cast<Float_t>(fP3.X());
    points[7] = static_cast<Float_t>(fP3.Y());
    points[8] = static_cast<Float_t>(fP3.Z());
    points[9] = static_cast<Float_t>(fP4.X());
    points[10] = static_cast<Float_t>(fP4.Y());
    points[11] = static_cast<Float_t>(fP4.Z());
}

// Fills TBuffer3D structure for segments and polygons.
// TODO: No idea what I'm doing here
void TGeoTet::SetSegsAndPols(TBuffer3D& buff) const
{
    const Int_t c = GetBasicColor();

    if (buff.NbSegs() != 6 || buff.NbPols() != 4)
    {
        buff.SetRawSizes(4, 4 * 3, 6, 6 * 3, 4, 4 * 5);
    }

    // Segs = Connections between Points (vertices)
    // Connection from point with id=0 to point with id=1. This connection gets id=0
    buff.fSegs[0] = c;
    buff.fSegs[1] = 0;
    buff.fSegs[2] = 1;

    buff.fSegs[3] = c + 1;
    buff.fSegs[4] = 1;
    buff.fSegs[5] = 2;

    buff.fSegs[6] = c + 2;
    buff.fSegs[7] = 0;
    buff.fSegs[8] = 2;

    buff.fSegs[9] = c + 3;
    buff.fSegs[10] = 0;
    buff.fSegs[11] = 3;

    buff.fSegs[12] = c + 4;
    buff.fSegs[13] = 1;
    buff.fSegs[14] = 3;

    buff.fSegs[15] = c + 6;
    buff.fSegs[16] = 2;
    buff.fSegs[17] = 3;

    // Pols = Polygons making up the faces, formed by connections
    // Faces here have always three lines
    buff.fPols[0] = c;
    buff.fPols[1] = 3; // number of lines
    buff.fPols[2] = 0;
    buff.fPols[3] = 1;
    buff.fPols[4] = 2;

    buff.fPols[5] = c; // + 1;
    buff.fPols[6] = 3;
    buff.fPols[7] = 0;
    buff.fPols[8] = 3;
    buff.fPols[9] = 4;

    buff.fPols[10] = c; // + 2;
    buff.fPols[11] = 3;
    buff.fPols[12] = 1;
    buff.fPols[13] = 4;
    buff.fPols[14] = 5;

    buff.fPols[15] = c;
    buff.fPols[16] = 3;
    buff.fPols[17] = 3;
    buff.fPols[18] = 2;
    buff.fPols[19] = 5;
}

const TBuffer3D& TGeoTet::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
    static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

    TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

    if (reqSections & TBuffer3D::kRawSizes)
    {
        if (buffer.SetRawSizes(4, 4 * 3, 6, 6 * 3, 4, 4 * 5))
        {
            buffer.SetSectionsValid(TBuffer3D::kRawSizes);
        }
    }

    if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes))
    {
        SetPoints(buffer.fPnts);
        if (!buffer.fLocalFrame)
        {
            TransformPoints(buffer.fPnts, buffer.NbPnts());
        }

        SetSegsAndPols(buffer);
        buffer.SetSectionsValid(TBuffer3D::kRaw);
    }

    return buffer;
}

ClassImp(TGeoTet);
