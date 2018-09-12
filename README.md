# splots

'heppy3.C' and 'heppy3.h' : Selection cuts and input variables for BDT (output : 'MuonEG_49.root' for MC and 'TTbar_14.root').

'tmva.h' and 'r_tmva.h' : BDT output for leading, sub-leading and others jet (output : 'TMVA_output.root' for MC and 'r_TMVA_output.root' for data)

'splots_hist_mc.C' and 'splots_hist_r.C' : histograms of BDTs of jets (output : Updated files 'TMVA_output.root' for MC and 'r_TMVA_output.root' for data)

'splots_tree_mc.C' and 'splots_tree_r.C' : pt_binned MC and data BDT templates and vectors->doubles tuple conversion (From this stage we lose the count of events.) (output : 'TMVA_output_splots.root' for MC and 'r_TMVA_output_splots.root' for data)

'splots_pt_mc.C' and ''splots_pt_r.C' : plots pt_binned fits to MC and data BDT using MC templates (output : folders 'bdt_fit_mc' and 'bdt_fit_r').

'splots_tree_rho.C' : compute rho_weight and add a branch (output : Update 'TMVA_output_splots.root' file)

'splots.C' : splots for pt_inclusive sample (ouptut : 'qreg_file.root')

'splots_*_**.C' : splots for pt bin from * to ** . (output : 'qreg_file_*_**.root')

'sig_bg_plots.C' : plots of data/MC comparision,sig/bg comparision of i/p variables of BDT and BDT itself. Also plots of rho_reweighting and sweights. This code gives all the plots for my thesis initial chapters. (folders : sig_bg_KIN, data_MC_KIN, rho_reweight, sweights)

IMP: In initial days, I have many times drawn a histogram in ROOT, saved it and latter imported it in RooFit. This is rubbish. Use createHistogram method in RooFit directly.
