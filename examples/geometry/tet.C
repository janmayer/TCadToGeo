// Demonstrate basic display
void tet()
{
    auto a = new TGeoTet({ 0, 0, 0 }, { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 });

    gStyle->SetCanvasPreferGL(kTRUE);
    a->Draw("ogl");
}
