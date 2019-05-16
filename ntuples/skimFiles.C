void skimFiles(TString filename = "", TString cut = "")
{
    TFile *of = new TFile(filename);
    TTree *ot = (TTree*)of->Get("PDsecondTree");
    TString filename_( filename(0, filename.Index(".root")) );
    TFile *nf = new TFile(filename_ + "_skim.root", "RECREATE");
    TTree *nt = ot->CopyTree(cut);
    nt->Write();
    nf->Close();
}