// Demonstrate Basic Geometry Definition
void tet2()
{
    auto geom = new TGeoManager("World", "IES");
    auto matVac = new TGeoMaterial("Vacuum", 0, 0, 0);
    auto medVac = new TGeoMedium("Vacuum", 0, matVac);

    auto world = geom->MakeBox("World", medVac, 200., 200., 200.);
    world->SetLineColor(kMagenta);
    geom->SetTopVolume(world);

    auto Sha1 = new TGeoTet({ 0, 0, 0 }, { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 });
    auto Vol1 = new TGeoVolume("tetrahedron", Sha1, medVac);
    Vol1->SetLineColor(kBlue);
    world->AddNode(Vol1, 0);

    auto Sha2 = new TGeoTet({ 0, 0, 0 }, { -1, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 });
    auto Vol2 = new TGeoVolume("tetrahedron2", Sha2, medVac);
    Vol2->SetLineColor(kGreen);
    world->AddNode(Vol2, 0);

    geom->CloseGeometry();
    geom->SetVisLevel(4);

    gStyle->SetCanvasPreferGL(kTRUE);
    world->Draw("ogl");
    auto v = (TGLViewer*)gPad->GetViewer3D();
    v->CurrentCamera().RotateRad(0, 0.1);
    v->SetStyle(TGLRnrCtx::kOutline);
    v->RequestDraw();
    v->SavePicture("tet2.png");

    geom->Export("tet2.root");
}
