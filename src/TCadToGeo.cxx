#include "TCadToGeo.h"
#include "TGeoManager.h"
#include "TGeoTet.h"
#include "tetgen.h"
#include <iostream>

TCadToGeo::TCadToGeo(const TString& fileName, const Double_t scale)
    : fScale(scale)
    , fIn(std::make_shared<tetgenio>())
    , fOut(std::make_shared<tetgenio>())
{
    fIn->load_ply(const_cast<char*>(fileName.Data()));
    Tetrahedralize();
}

TGeoVolumeAssembly* TCadToGeo::MakeVolumeAssembly(const TString& name, const TGeoMedium* med) const
{
    if (gGeoManager == nullptr)
    {
        throw std::runtime_error("TCadToGeo::MakeVolumeAssembly: No gGeoManager.");
    }

    auto assembly = gGeoManager->MakeVolumeAssembly(name);
    if (assembly == nullptr)
    {
        throw std::runtime_error("TCadToGeo::MakeVolumeAssembly: Could not create volume assembly.");
    }

    for (int i = 0; i < fOut->numberoftetrahedra; i++)
    {
        const int indexOffset = i * 4; // For a tetrahedron, fOut->numberofcorners == 4
        const TVector3 p1 = GetTetPoint(indexOffset);
        const TVector3 p2 = GetTetPoint(indexOffset + 1);
        const TVector3 p3 = GetTetPoint(indexOffset + 2);
        const TVector3 p4 = GetTetPoint(indexOffset + 3);

        auto tet = new TGeoTet(p1, p2, p3, p4);
        const TString nameTetVol = TString::Format("%s_tet_%d", name.Data(), i);
        auto tetVol = new TGeoVolume(nameTetVol, tet, med);
        assembly->AddNode(tetVol, 0);
    }
    return assembly;
}

TVector3 TCadToGeo::GetTetPoint(const int i) const
{
    return TVector3(fOut->pointlist[fOut->tetrahedronlist[i] * 3 + 0] * fScale,
                    fOut->pointlist[fOut->tetrahedronlist[i] * 3 + 1] * fScale,
                    fOut->pointlist[fOut->tetrahedronlist[i] * 3 + 2] * fScale);
}

void TCadToGeo::Tetrahedralize()
{
    auto behavior = new tetgenbehavior();
    behavior->nobisect = 1;
    behavior->plc = 1;
    // behavior->quality = 1;
    // behavior->coarsen = 0;
    // behavior->optlevel = 2;

    tetrahedralize(behavior, fIn.get(), fOut.get());
}

ClassImp(TCadToGeo)
