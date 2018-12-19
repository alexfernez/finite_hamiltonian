#include "Matrix_Hamiltonian.h"

// importantly, all units will be MeV and fm
const double C_LIGHT = 2.99792458e23;
const double H_BAR = 6.5821e-22;
const double H_BAR_C = 197.33;

using namespace std;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix_Double_Dynamic;

//////////////////////////////////////////////// Some useful helper functions ////////////////////////////////////////////////////////

int roundUp(double x) { //rounds a double up to the nearest integer
  int x_int = int(x);
  if (x_int < x) return x_int + 1; //this should always be the case, but just for safety
  else return x_int;
}

int C3n(int n) { //returns the number of ways to sum three squared non-negative integers to n. Implemented as O(n^(3/2)) algorithm
  int count = 0;
  int root_n_int = roundUp(sqrt(n));
  for (int nx=0; nx<=root_n_int; nx++) {
    for (int ny=0; ny<=root_n_int; ny++) {
      for (int nz=0; nz<=root_n_int; nz++) {
        if ((nx*nx + ny*ny + nz*nz) == n) count++;
      }
    }
  }
  return count;
} // as a note, I would imagine this function could be improved to be O(1) using a generating function (or some other method) to count C3n


//////////////////////////////////////////////// Class implementations ////////////////////////////////////////////////////////

bare_state::bare_state(double m) {
  if (m < 0) {
    cout << "Provide a positive mass for the bare state" << endl;
    exit(1);
  }
  mass = m;
}

double bare_state::M() {return mass;}

bool bare_state::operator ==(bare_state bs) {return (mass == bs.mass);}

/////////

bare_state_list::bare_state_list(bare_state bs) {bare_states.push_back(bs);}

void bare_state_list::Add(bare_state bs) {bare_states.push_back(bs);}

int bare_state_list::Num_BareStates() {return bare_states.size();}

bare_state bare_state_list::operator [](int i) {
  if (i < 0 || i >= bare_states.size()) {
    cout << "Bare state index out of range" << endl;
    exit(1);
  }
  return bare_states[i];
}

bare_state bare_state_list::Get_BareState(int i) {
  if (i < 0 || i >= bare_states.size()) {
    cout << "Bare state index out of range" << endl;
    exit(1);
  }
  return bare_states[i];
}

bool bare_state_list::operator ==(bare_state_list bsl) {
  if (bsl.Num_BareStates() != bare_states.size()) return false;
  for (int i=0; i<bsl.Num_BareStates(); i++) {
    if (!(bsl.Get_BareState(i) == bare_states[i])) return false;
  }
  return true;
}

/////////

class discrete_state {
  // Discrete state is defined by two particles' mass and momentum (equal magnitude in CM frame, also quantized since we are on a lattice)
  // The user does not need to see the implementation of the discrete states
  // Though the paper mentions using vectors k and n, the useful quantities in the states only ever involve |k| and |n|, so I will only have
  // a field for |n| (where k = 2pi/L * sqrt(n))
  double m1; // assumed to be in MeV/c^2
  double m2;
  int n; // n = nx^2 + ny^2 + nz^2

  public:
    discrete_state(double m1=0.0, double m2=0.0, int n=0);
    double M_1();
    double M_2();
    double k(double L); // k = p/h_bar: the same for each particle, since we are in CM frame
    double k2(double L);
    double E_1(double L); // energy of state in MeV: (potentially) differs between particles because of different mass values
    double E_2(double L);
    //bool operator ==(const discrete_state& ds); // I don't ever use this
    ~discrete_state() {}
};

discrete_state::discrete_state(double m1val, double m2val, int nval) {
  if (m1val < 0 || m2val < 0 || nval < 0) {
    cout << "Entered something wrong in discrete state parameters" << endl;
    exit(1);
  }
  m1 = m1val;
  m2 = m2val;
  n = nval;
}

double discrete_state::M_1() {return m1;}

double discrete_state::M_2() {return m2;}

double discrete_state::k(double L) {return (2.0*M_PI/L)*sqrt(n);}

double discrete_state::k2(double L) {return k(L)*k(L);}

double discrete_state::E_1(double L) {return sqrt(m1*m1 + H_BAR_C*H_BAR_C*k2(L));}

double discrete_state::E_2(double L) {return sqrt(m2*m2 + H_BAR_C*H_BAR_C*k2(L));}

//bool discrete_state::operator ==(const discrete_state& ds) {return (m1==ds.m1 && m2==ds.m2 && n==ds.n);}

//////////

two_particle_channel::two_particle_channel(double m1, double m2) {
  if (m1 < 0.0 || m2 < 0.0) {
    cout << "Provide positive masses for two particle channel" << endl;
    exit(1);
  }
  mass_1 = m1; mass_2 = m2;
}

double two_particle_channel::M_1() {return mass_1;}

double two_particle_channel::M_2() {return mass_2;}

bool two_particle_channel::operator ==(two_particle_channel tpc) {
  return (mass_1 == tpc.mass_1 && mass_2 == tpc.mass_2) || (mass_1 == tpc.mass_2 && mass_2 == tpc.mass_1);
}

//////////

two_particle_channel_list::two_particle_channel_list(two_particle_channel tpc) {two_particle_channels.push_back(tpc);}

void two_particle_channel_list::Add(two_particle_channel tpc) {two_particle_channels.push_back(tpc);}

int two_particle_channel_list::Num_Channels() {return two_particle_channels.size();}

two_particle_channel two_particle_channel_list::operator [](int i) {
  if (i < 0 || i >= two_particle_channels.size()) {
    cout << "Two particle channel index out of range" << endl;
    exit(1);
  }
  return two_particle_channels[i];
}

two_particle_channel two_particle_channel_list::Get_TwoParticleChannel(int i) {
  if (i < 0 || i >= two_particle_channels.size()) {
    cout << "Two particle channel index out of range" << endl;
    exit(1);
  }
  return two_particle_channels[i];
}

bool two_particle_channel_list::operator ==(two_particle_channel_list tpcl) {
  if (tpcl.Num_Channels() != two_particle_channels.size()) return false;
  for (int i=0; i<tpcl.Num_Channels(); i++) {
    if (!(tpcl.Get_TwoParticleChannel(i) == two_particle_channels[i])) return false;
  }
  return true;
}

//////////

Hamiltonian::Hamiltonian(int n_val, double L_val, bare_state_list bsl_val, two_particle_channel_list tpcl_val, g_func g_val, v_func v_val) {
// Sets the Hamiltonian with given parameters according to paper by Wu, Lee, Thomas, and Young
// This is really the most important function in the implementation

  //check all the input values are okay, and then set the fields
  if (n_val <= 0) {
    cout << "Provide a positive number of momenta (n>0)" << endl;
    exit(1);
  }
  n = n_val;
  if (L_val <= 0) {
    cout << "Provide a positive lattice spacing (L>0.0)" << endl;
    exit(1);
  }
  L = L_val;
  if (bsl_val.Num_BareStates() == 0) {
    cout << "Provide a non-zero number of bare states" << endl;
    exit(1);
  }
  bsl = bsl_val;
  if (tpcl_val.Num_Channels() == 0) {
    cout << "Provide a non-zero number of two particle channels" << endl;
    exit(1);
  }
  tpcl = tpcl_val;
  if (g_val != nullptr) g = g_val;
  else g = default_g_potential;
  if (v_val != nullptr) v = v_val;
  else v = default_v_potential;

  // first, get the size of the Hamiltonian and initialize the entries
  int num_bare = bsl.Num_BareStates();
  int num_chan = tpcl.Num_Channels();
  int size = n*num_chan + num_bare;
  for (int i=0; i<size; i++) {
    vector<double> row;
    for (int j=0; j<size; j++) row.push_back(0); // ensures that Hamiltonian is a square matrix
    elements.push_back(row);
  }


  // now, add in the bare states (part of H_0)
  for (int i=0; i<num_bare; i++) {
    elements[i][i] += bsl[i].M();
  }

  // Note: when entering in a discrete state entry in elements, the index will be offset by num_bare

  // add in the discrete state energies (part of H_0)
  for (int ni=0; ni<n; ni++) {
    for (int alpha=0; alpha<num_chan; alpha++) {
      discrete_state ds(tpcl.Get_TwoParticleChannel(alpha).M_1(), tpcl.Get_TwoParticleChannel(alpha).M_2(), ni);
      elements[num_bare + (num_chan*ni)+alpha][num_bare + (num_chan*ni)+alpha] += ds.E_1(L) + ds.E_2(L);
    }
  }

  // add in the vertex interaction (part of H_I)
  for (int i=0; i<num_bare; i++) {
    for (int ni=0; ni<n; ni++) {
      for (int alpha=0; alpha<num_chan; alpha++) {
        double k = (2.0*M_PI/L)*sqrt(ni); // could make a discrete state to get this value, but this is faster
        elements[num_bare + (num_chan*ni)+alpha][i] += (sqrt(C3n(ni)/(4.0*M_PI)) * pow(2.0*M_PI/L, 1.5) * g(k, bsl[i], tpcl[alpha]));
        elements[i][num_bare + (num_chan*ni)+alpha] += (sqrt(C3n(ni)/(4.0*M_PI)) * pow(2.0*M_PI/L, 1.5) * g(k, bsl[i], tpcl[alpha]));
      }
    }
  }

  // add in two particle two particle interaction (part of H_I)
  for (int ni=0; ni<n; ni++) {
    for (int nj=0; nj<n; nj++) {
      for (int alpha1=0; alpha1<num_chan; alpha1++) {
        for (int alpha2=0; alpha2<num_chan; alpha2++) {
          double ki = (2.0*M_PI/L)*sqrt(ni);
          double kj = (2.0*M_PI/L)*sqrt(nj);
          double v_elem = sqrt(C3n(ni)*C3n(nj)/(4.0*M_PI)) * pow(2.0*M_PI/L, 3.0) * v(ki, kj, tpcl[alpha1], tpcl[alpha2]);
          elements[num_bare + (num_chan*ni)+alpha1][num_bare + (num_chan*nj)+alpha2] += v_elem;
        }
      }
    }
  }

  // solve eigenvalue equation and set evecs and evals to eigenvectors and eigenvalues
  Matrix_Double_Dynamic MatHam(size, size);
  for (int i=0; i<size; i++) {
    for (int j=0; j<size; j++) {
      MatHam(i,j) = elements[i][j];
    }
  }
  // When n becomes very large, the next line of code, where the eigenvalues and eigenvectors are solved for, will likely take the most time
  Eigen::EigenSolver<Matrix_Double_Dynamic> eigensolver(MatHam, true); // true tells Eigensolver to compute eigenvectors

  for (int i=0; i<size; i++) {
    evals.push_back(eigensolver.eigenvalues()(i).real()); // note that eigenvalues must be real, anyway, but Eigen code stores them as complex
    vector<double> evec;
    for (int j=0; j<size; j++) evec.push_back(eigensolver.eigenvectors()(i,j).real());
    evecs.push_back(evec);
  }
}

int Hamiltonian::Get_n() {return n;}

double Hamiltonian::Get_L() {return L;}

bare_state_list Hamiltonian::Get_BareStateList() {return bsl;}

two_particle_channel_list Hamiltonian::Get_TwoParticleChannelList() {return tpcl;}

vector<double> Hamiltonian::Get_Eigenvector(int i) {return evecs[i];}

double Hamiltonian::Get_Eigenvalue(int i) {return evals[i];}

vector< vector<double> > Hamiltonian::Get_Elements() {return elements;}

void Hamiltonian::Set_n(int n_new) {
  if (n_new == n) {
    cout << "Hamiltonian already uses this number of momenta" << endl;
    return;
  }

  //else, recall the constructor with new n value and set current object to it
  Hamiltonian H_new(n_new, L, bsl, tpcl, g, v);
  *this = H_new;
}

void Hamiltonian::Set_L(double L_new) {
  if (L == L_new) {
    cout << "Hamiltonian already uses this lattice spacing" << endl;
    return;
  }

  //else, recall the constructor with the new L value
  Hamiltonian H_new(n, L_new, bsl, tpcl, g, v);
  *this = H_new;
}

void Hamiltonian::Set_BareStates(bare_state_list bsl_new) {
  if (bsl == bsl_new) {
    cout << "Hamiltonian already uses this bare state list" << endl;
    return;
  }
  Hamiltonian H_new(n, L, bsl_new, tpcl, g, v);
  *this = H_new;
}

void Hamiltonian::Set_TwoParticleChannels(two_particle_channel_list tpcl_new) {
  if (tpcl == tpcl_new) {
    cout << "Hamiltonian already uses this two particle channel list" << endl;
    return;
  }
  Hamiltonian H_new(n, L, bsl, tpcl_new, g, v);
  *this = H_new;
}

void Hamiltonian::Set_g(g_func g_new) {
  Hamiltonian H_new(n, L, bsl, tpcl, g_new, v);
  *this = H_new;
}

void Hamiltonian::Set_v(v_func v_new) {
  Hamiltonian H_new(n, L, bsl, tpcl, g, v_new);
  *this = H_new;
}

double Hamiltonian::Eval_g(double k, bare_state bs, two_particle_channel tpc) {return (*g)(k, bs, tpc);}

double Hamiltonian::Eval_v(double k1, double k2, two_particle_channel tpc1, two_particle_channel tpc2) {return (*v)(k1, k2, tpc1, tpc2);}

int Hamiltonian::Get_Dimension() {return elements.size();}

vector<double> Hamiltonian::operator [](int i) {
  if (i < 0 || i >= elements.size()) {
    cout << "Hamiltonian row index out of range" << endl;
    exit(1);
  }
  return elements[i];
}

double Hamiltonian::Get_MatrixElement(int i, int j) {
  if (i < 0 || i >= elements.size()) {
    cout << "Hamiltonian row index out of range" << endl;
  }
  if (j < 0 || j >= elements.size()) { // uses the fact that H is a square matrix
    cout << "Hamiltonian col index out of range" << endl;
  }
  return elements[i][j];
}

// prints out the upper left hand square block of matrix Hamiltonian, up to row max and col max
void Hamiltonian::Print_SquareBlock(int max) {
  cout << "Matrix Hamiltonian up to row " << max << " and column " << max << ":" << endl;
  for (int i=0; i<max; i++) {
    for (int j=0; j<max; j++) {
      cout << elements[i][j] << "  ";
    }
    cout << endl;
  }
}

void Hamiltonian::Sort_Eigenvalues() { // first makes a copy of evals and then sorts that, later will swap elements in evecs and evals
  vector<double> evals_sorted(evals);
  sort(evals_sorted.begin(), evals_sorted.end());

  //check which elements have changed in evals, and then change them in evecs: efficiency of algorithm relies on evals being mostly sorted
  for (int i=0; i<evals_sorted.size(); i++) {
    if (evals[i] != evals_sorted[i]) {
      for (int j=0; j<evals.size(); j++) {
        if (evals[j] == evals_sorted[i]) {
          vector<double> temp_evec(evecs[i]);
          double temp_eval(evals[i]);
          evecs[i] = evecs[j];
          evecs[j] = temp_evec;
          evals[i] = evals[j];
          evals[j] = temp_eval;
          break; // if you've already found the right thing to switch with in evals_sorted, no need to keep looking through that list
        }
      }
    }
  }
}

void Hamiltonian::Print_Eigen() { // prints out eval and evec pairs, in current order (may not be sorted unless Sort_Eigenvalues() called)
  int size = evals.size();
  for (int i=0; i<size; i++) {
    cout << "Eigenvalue " << i << ": " << evals[i] << endl;
    cout << "Eigenvector " << i << ": (";
    for (int j=0; j<size; j++) {
      if (j != size-1) cout << evecs[i][j] << ", ";
      else cout << evecs[i][j] << ")" << endl;
    }
  }
}

void Hamiltonian::Print_Eigenvalues() { // prints sorted eigenvalues, but doesn't sort evals or evecs fields
  vector<double> evals_sorted(evals);
  sort(evals_sorted.begin(), evals_sorted.end());

  cout << "Eigenvalues of current Matrix Hamiltonian:" << endl;
  for (int i=0; i<evals_sorted.size(); i++) {
    cout << evals_sorted[i] << endl;
  }
}

// allows you to print out Hamiltonian using cout-- note that columns won't be aligned perfectly, however
ostream& operator <<(ostream& out, Hamiltonian& H) {
  out << "Matrix Hamiltonian:" << endl;
  for (int i=0; i<H.Get_Dimension(); i++) {
    for (int j=0; j<H.Get_Dimension(); j++) {
      out << H[i][j] << "  ";
    }
    out << endl;
  }
  return out;
}


///////////////////////////////////////////// User function implementation  //////////////////////////////////////////////

double default_g_potential(double k, bare_state bs, two_particle_channel tpc) {
// g according to Wu, Lee, Thomas, Young paper eqn 20, constants from table 1. The constants are given for some bare state that the table asserts
// has mass 700 MeV/c^2. In order to be able to input different bare states than this value, the function will have to be edited, or another
// function will have to be written.

  bare_state bs_compare(700.0);
  if (!(bs == bs_compare)) {
    cout << "Constants not defined for any bare state other than bare state with mass 700 MeV/c^2: Try using different g potential..." << endl;
    exit(1);
  }

  // the paper also only defines constants for two channels: pipi and kkbar. Make sure the input channel is one of these, and use the approrpriate
  // constant depending on what tpc is
  two_particle_channel pipi_tpc(134.9767, 134.9767);
  two_particle_channel kkbar_tpc(497.611, 497.611);
  two_particle_channel_list tpcl_compare;
  tpcl_compare.Add(pipi_tpc);
  tpcl_compare.Add(kkbar_tpc);
  bool tpc_okay = false;
  for (int alpha=0; alpha<tpcl_compare.Num_Channels(); alpha++) if (tpcl_compare[alpha] == tpc) tpc_okay = true;
  if (!tpc_okay) {
    cout << "Constants not defined for input two particle channel: Try using different g potential..." << endl;
    exit(1);
  }

  // based on the input bs and tpc, decide what the constants used in the potential should be
  double g_const;
  double c_const;
  if (tpc == pipi_tpc) {
    g_const = 2.0;
    c_const = 0.6722;
  }
  if (tpc == kkbar_tpc) {
    g_const = 0.6451;
    c_const = 1.0398;
  }

  double returnval = g_const/sqrt(M_PI);
  returnval *= pow(1 + pow(c_const*k, 2.0), -1.0);

  return returnval;
}

double default_v_potential(double k1, double k2, two_particle_channel tpc1, two_particle_channel tpc2) {
// v according to Wu, Lee, Thomas, Young paper eqn 21, constants from table 1

  // the paper only defines constants two channels; make sure both tpc1 and tpc2 are one of these channels
  two_particle_channel pipi_tpc(134.9767, 134.9767);
  two_particle_channel kkbar_tpc(497.611, 497.611);
  two_particle_channel_list tpcl_compare;
  tpcl_compare.Add(pipi_tpc);
  tpcl_compare.Add(kkbar_tpc);
  bool tpc1_okay = false;
  bool tpc2_okay = false;
  for (int alpha=0; alpha<tpcl_compare.Num_Channels(); alpha++) if (tpcl_compare[alpha] == tpc1) tpc1_okay = true;
  for (int alpha=0; alpha<tpcl_compare.Num_Channels(); alpha++) if (tpcl_compare[alpha] == tpc2) tpc2_okay = true;
  if (!(tpc1_okay && tpc2_okay)) {
    cout << "Constants not defined for input two particle channel: Try using different g potential..." << endl;
    exit(1);
  }

  // based on the two input tpc, decide what the constants should be: this should be done a little more cleverly if there are a lot of
  // two particle channels to consider
  double G_const;
  double d1_const;
  double d2_const;
  if (tpc1 == pipi_tpc && tpc2 == pipi_tpc) {
    G_const = 2.4998;
    d1_const = 0.244;
    d2_const = 0.244;
  }
  if (tpc1 == kkbar_tpc && tpc2 == kkbar_tpc) {
    G_const = 0.02;
    d1_const = 0.1;
    d2_const = 0.1;
  }
  if (tpc1 == pipi_tpc && tpc2 == kkbar_tpc) {
    G_const = 0.35;
    d1_const = 0.244;
    d2_const = 0.1;
  }
  if (tpc1 == kkbar_tpc && tpc2 == pipi_tpc) {
    G_const = 0.35;
    d1_const = 0.1;
    d2_const = 0.244;
  }

  double returnval = G_const/pow(134.9767, 2); // 134.9767 is mass of neutral pion; this is because the entire paper focuses on pipi scattering
  returnval *= pow(1 + pow(d1_const*k1, 2.0), -1.0);
  returnval *= pow(1 + pow(d2_const*k2, 2.0), -1.0);

  return returnval;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
  /*
  // A test
  bare_state bs(700.0);
  bare_state_list bsl(bs);
  two_particle_channel pipi_tpc(134.9767, 134.9767);
  two_particle_channel kkbar_tpc(497.611, 497.611);
  two_particle_channel_list tpcl;
  tpcl.Add(pipi_tpc);
  tpcl.Add(kkbar_tpc);
  int n = 6; //number of momenta considered: Hamiltonian will be size (#channels)*n + (#bare states)
  double L = 4.0; // in fm
  Hamiltonian H(n, L, bsl, tpcl);

  H.Sort_Eigenvalues();
  H.Print_Eigen();
  H.Print_SquareBlock(8);
  */

  return 0;
}
