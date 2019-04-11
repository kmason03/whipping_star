void hand_make_matrix(){

TFile *f = new TFile("unit1b_matrix.root","recreate");
TMatrixD *m = new TMatrixD(6,6);
m->Zero();

(*m)(0,0) = 0.2;
(*m)(1,1) = 0.2;

(*m)(2,2) = 0.15;
(*m)(3,3) = 0.15;

(*m)(4,4) = 0.1;
(*m)(5,5) = 0.1;

(*m)(0,2) = 0.21;
(*m)(2,0) = 0.21;
(*m)(1,3) = 0.21;
(*m)(3,1) = 0.21;

(*m)(2,4) = 0.08;
(*m)(4,2) = 0.08;
(*m)(3,5) = 0.08;
(*m)(5,3) = 0.08;

(*m)(0,4) = 0.12;
(*m)(4,0) = 0.12;
(*m)(1,5) = 0.12;
(*m)(5,1) = 0.12;


m->Print();
m->Write();
f->Close();

/*
TFile *f = new TFile("unit1a_matrix.root","recreate");
TMatrixD *m = new TMatrixD(6,6);
m->Zero();

(*m)(0,0) = 0.2;
(*m)(1,1) = 0.15;
(*m)(2,2) = 0.1;

(*m)(0,1) = 0.21;
(*m)(1,0) = 0.21;

(*m)(0,2) = 0.12;
(*m)(2,0) = 0.12;

(*m)(1,2) = 0.08;
(*m)(2,1) = 0.08;

(*m)(3,3) = 0.2;
(*m)(4,4) = 0.15;
(*m)(5,5) = 0.1;

(*m)(3,4) = 0.21;
(*m)(4,3) = 0.21;

(*m)(3,5) = 0.12;
(*m)(5,3) = 0.12;

(*m)(4,5) = 0.08;
(*m)(5,4) = 0.08;

m->Print();
m->Write();
f->Close();

*/
}