// -*- mode: c++ -*-
/*
 * Groom.h
 *
 * Copyright 2018 Andrii Verbytskyi <andriish@mppmu.mpg.de>
 * Max-Planck Institut für Physik
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#ifndef GROOM_H
#define GROOM_H

template <class T> double MOMX(T &A) { return A.px(); }
template <class T> double MOMY(T &A) { return A.py(); }
template <class T> double MOMZ(T &A) { return A.pz(); }
template <class T> double MOME(T &A) { return A.e(); }

// #ifdef ROOT_TLorentzVector
template <> double MOMX(TLorentzVector  &A) { return A.Px(); }
template <> double MOMY(TLorentzVector  &A) { return A.Py(); }
template <> double MOMZ(TLorentzVector  &A) { return A.Pz(); }
template <> double MOME(TLorentzVector  &A) { return A.E(); }
// #endif

/* #ifdef HEPMC_GENPARTICLE_H */
/* template <> double MOMX(HepMC::GenParticlePtr  &A) { return A->momentum().px(); }; */
/* template <> double MOMY(HepMC::GenParticlePtr  &A) { return A->momentum().py(); }; */
/* template <> double MOMZ(HepMC::GenParticlePtr  &A) { return A->momentum().pz(); }; */
/* template <> double MOME(HepMC::GenParticlePtr  &A) { return A->momentum().e(); }; */
/* #endif */


template <class T>
bool GroomInput( std::vector<T> A, const double beta, const double zcut,
		 std::vector<fastjet::PseudoJet> & C1,
		 std::vector<fastjet::PseudoJet> & C2,
		 fastjet::PseudoJet & AX1,
		 fastjet::PseudoJet & AX2 ) {
  std::vector<fastjet::PseudoJet> fjInputsA;
  for( size_t i=0; i<A.size(); i++ ) {
    fastjet::PseudoJet part_momA;
    part_momA.reset( MOMX<T>( A[i] ), MOMY<T>( A[i] ), MOMZ<T>( A[i] ), MOME<T>( A[i] ) );
    fjInputsA.push_back( part_momA );
  }
  fastjet::JetDefinition jet_def( fastjet::ee_genkt_algorithm, 3*M_PI/2, 0.0 );
  fastjet::ClusterSequence clust_seq_A( fjInputsA, jet_def );
  std::vector<fastjet::PseudoJet> jets_A= clust_seq_A.exclusive_jets(2);
  double R0 = 1.0;
  fastjet::contrib::SoftDrop
    fSoftDrop( beta, zcut,
	       fastjet::contrib::RecursiveSymmetryCutBase::theta_E,R0,
	       std::numeric_limits<double>::infinity(),
	       fastjet::contrib::RecursiveSymmetryCutBase::larger_E );
  fastjet::PseudoJet T1= fSoftDrop( jets_A[0] );
  fastjet::PseudoJet T2= fSoftDrop( jets_A[1] );
  if( !T1.has_constituents() ) return false;
  if( !T2.has_constituents() ) return false;
  C1= T1.constituents();
  C2= T2.constituents();
  AX1= jets_A[0];
  AX2= jets_A[1];
  return true;
}


template <class T>
std::map<std::string,double> GetGroomedObservables( std::vector<T> C1,
						    std::vector<T> C2,
						    const T AX1,
						    const T AX2,
						    const std::string& prefix="" ) {
  /* This is a main function for all observables*/
  std::map<std::string, double> P;
  P["Ctot"]=0;
  P["T1"]=0;
  P["Product1"]=0;
  P["T2"]=0;
  P["Product2"]=0;
  for( size_t i=0; i<C1.size(); i++ ) {
    P["T1"]+= std::abs(AX1.px()*C1[i].px()+AX1.py()*C1[i].py()+AX1.pz()*C1[i].pz());
    double p1i= std::sqrt(C1[i].px()*C1[i].px()+C1[i].py()*C1[i].py()+C1[i].pz()*C1[i].pz());
    P["Product1"]+= p1i;
  }
  for( size_t i=0; i<C2.size(); i++ ) {
    P["T2"]+= std::abs(AX2.px()*C2[i].px()+AX2.py()*C2[i].py()+AX2.pz()*C2[i].pz());
    double p2i= std::sqrt(C2[i].px()*C2[i].px()+C2[i].py()*C2[i].py()+C2[i].pz()*C2[i].pz());
    P["Product2"]+= p2i;
  }
  std::vector<T > Ctot;
  Ctot.insert(Ctot.end(),C1.begin(),C1.end());
  Ctot.insert(Ctot.end(),C2.begin(),C2.end());

  for( size_t i=0; i<Ctot.size(); i++ ) {
    double pi= std::sqrt(Ctot[i].px()*Ctot[i].px()+Ctot[i].py()*Ctot[i].py()+Ctot[i].pz()*Ctot[i].pz());
    for( size_t j=0; j<Ctot.size(); j++ ) {
      double pj= std::sqrt(Ctot[j].px()*Ctot[j].px()+Ctot[j].py()*Ctot[j].py()+Ctot[j].pz()*Ctot[j].pz());
      P["Ctot"]+= 3.0/2.0*pi*pj*(1.0-  std::pow((Ctot[i].px()*Ctot[j].px()+Ctot[i].py()*Ctot[j].py()+Ctot[i].pz()*Ctot[j].pz())/(pi*pj),2));
    }
  }
  std::map<std::string,double> res;
  /* Set status to other value in case something goes wrong*/
  res["status"]=0;
  P["norm1"]=std::sqrt(AX1.px()*AX1.px()+AX1.py()*AX1.py()+AX1.pz()*AX1.pz());
  P["norm2"]=std::sqrt(AX2.px()*AX2.px()+AX2.py()*AX2.py()+AX2.pz()*AX2.pz());
  res[prefix+"T"]=(P["T1"]/P["norm1"]+P["T2"]/P["norm2"])/(P["Product1"]+P["Product2"]);
  res[prefix+"C"]=P["Ctot"]/std::pow(P["Product1"]+P["Product2"],2.0);
  // printf("%f %f %f %f %f\n",res["C"], P["T1"],P["T2"],P["Product1"],P["Product2"]);
  return res;
}
#endif

