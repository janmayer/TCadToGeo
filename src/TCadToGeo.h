#ifndef TCADTOGEO_H
#define TCADTOGEO_H

#include "TGeoVolume.h"
#include "TVector3.h"
#include <memory>

class tetgenio;

class TCadToGeo
{
  public:
    TCadToGeo(const TString& fileName, const Double_t scale = 1);
    virtual ~TCadToGeo() = default;

    TGeoVolumeAssembly* MakeVolumeAssembly(const TString& name, const TGeoMedium* med) const;

  private:
    Double_t fScale;
    std::shared_ptr<tetgenio> fIn;
    std::shared_ptr<tetgenio> fOut;

    TVector3 GetTetPoint(const int i) const;
    void Tetrahedralize();

    ClassDef(TCadToGeo, 1);
};

#endif // TCADTOGEO_H
