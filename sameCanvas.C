// simple_plot.C
// Simple macro to plot 4 histograms from 4 root files in same canvas

void Ncl_sameCanvas()
{
    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "4 Histograms", 1000, 800);
    gStyle->SetOptStat(0);
    
    // Define file names and histogram names
    TString files[4] = {"/home/msahil/Desktop/Education/PhD_Doctor_of_Philisophy/research/digitization_data/clustring/dch_updates_2025_11_16/Proton_data/protondch_digi_alg_N1000_debug.root",
    "/home/msahil/Desktop/Education/PhD_Doctor_of_Philisophy/research/digitization_data/clustring/dch_updates_2025_11_16/Muon_data/muondch_digi_alg_N1000_debug.root",
    "/home/msahil/Desktop/Education/PhD_Doctor_of_Philisophy/research/digitization_data/clustring/dch_updates_2025_11_16/Pion_data/piondch_digi_alg_N1000_debug.root",
    "/home/msahil/Desktop/Education/PhD_Doctor_of_Philisophy/research/digitization_data/clustring/dch_updates_2025_11_16/Kaon_data/kaondch_digi_alg_N1000_debug.root"};
    
    TString hist_Ncl[4] = {"hNcl_perStep", "hNcl_perStep", "hNcl_perStep", "hNcl_perStep"}; // Change if different names
    
    TString legends[4] = {"Ncl Proton", "Ncl Muon", "Ncl Pion", "Ncl Kaon"};
    int colors[4] = {kRed, kBlue, kGreen, kBlack};
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);

    // Loop over files
    for (int i = 0; i < 4; i++) {
        // Open file
        TFile *file = TFile::Open(files[i]);
        if (!file)
        {
            cout << "Error: Cannot open " << files[i] << endl;
            continue;
        }

        TH1F *hist = (TH1F*)file->Get(hist_Ncl[i]);
        if (!hist)
        {
            cout << "Error: Cannot find histogram in " << files[i] << endl;
            continue;
        }

        hist->SetLineColor(colors[i]);
        hist->SetLineWidth(3);

        if (i == 0)
        {
            hist->Draw("HIST");
            hist->SetTitle("Ncl Distribution");
            hist->GetXaxis()->SetTitle("Ncl");
            hist->GetYaxis()->SetTitle("Counts");
            
            // OPTIONAL: AXIS LABEL SIZE BADA SAKTE HO
            hist->GetXaxis()->SetTitleSize(0.04);
            hist->GetYaxis()->SetTitleSize(0.04);
        } 
        else {
            hist->Draw("HIST SAME");
        }
        
        // Add to legend
        leg->AddEntry(hist, legends[i], "l");
    }
    

    leg->Draw();
    c1->Update();
    
    c1->SaveAs("../analysis_plots/particle_Ncl_comparison.png");
}
