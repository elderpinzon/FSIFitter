///////////////////////////////////////////////////////////
// Simple macro to read the summary csv file into the root
// file that will be used later on by the FSIFitter
// Note: the input and output file names are hard-coded
///////////////////////////////////////////////////////////

void read_summary_file() {

  // input text file
  string fname = "scan_all_c_o_al_fe_cu_pb.dat";
  string fname_root = "scan_all_c_o_al_fe_cu_pb.root";

  TFile *f = new TFile(fname_root.data(),"RECREATE");
  TTree *tree = new TTree("t","Scan of FSI parameters");
  tree->ReadFile(fname.data(),"target/I:mom/D:pid/I:FEFQE/D:FEFABS/D:FEFCX/D:FEFINEL/D:FEFQEH/D:FEFCXH/D:FEFALL/D:nall/I:nreac/I:nqe/I:nabs/I:ncx/I:ndcx/I:nhadr/I:xreac/D:xqe/D:xabs/D:xcx/D:xdcx/D:xhadr/D",' ');
  f->Write();
}
