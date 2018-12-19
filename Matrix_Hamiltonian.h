#ifndef Matrix_Hamiltonian_h
#define Matrix_Hamiltonian_h

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <vector>
#include "Eigen/Eigenvalues"

//////////////////////////////// User Classes

class bare_state {
  // May add more fields, if desired. Though, for program's usage, the mass is all that matters
  double mass; // assumed to be in MeV/c^2

  public:
    bare_state(double m=0.0);
    double M();
    bool operator ==(bare_state bs);
    ~bare_state() {}
};

class bare_state_list {
  std::vector<bare_state> bare_states;

  public:
    bare_state_list() {}
    bare_state_list(bare_state bs);
    void Add(bare_state bs); //adds to back of bare_states vector
    int Num_BareStates();
    bare_state operator [](int i);
    bare_state Get_BareState(int i); //the same as above, though more convenient if using a bare_state_list* object and can't make assignments
    bool operator ==(bare_state_list bsl);
    ~bare_state_list() {}
};

class two_particle_channel {
  // Like for the bare state, you can define more fields if you desire
  double mass_1; // Masses assumed to be in MeV/c^2
  double mass_2;

  public:
    two_particle_channel(double m1=0.0, double m2=0.0);
    double M_1();
    double M_2();
    bool operator ==(two_particle_channel tpc);
    ~two_particle_channel() {}
};

class two_particle_channel_list {
  std::vector<two_particle_channel> two_particle_channels;

  public:
    two_particle_channel_list() {}
    two_particle_channel_list(two_particle_channel tpc);
    void Add(two_particle_channel tpc);
    int Num_Channels();
    two_particle_channel operator [](int i);
    two_particle_channel Get_TwoParticleChannel(int i); //more convenient for two_particle_channel_list* object, but can't make assignments
    bool operator ==(two_particle_channel_list tpcl);
    ~two_particle_channel_list() {}
};

typedef double (*g_func)(double, bare_state, two_particle_channel);
typedef double (*v_func)(double, double, two_particle_channel, two_particle_channel);

class Hamiltonian {
  // Matrix Hamiltonian, as proposed by the Australian lattice group's papers. The size of the Hamiltonian will be n*size(tpcl) + size(bsl)
  int n; // number of momenta considered
  double L; // assumed to be in fm
  bare_state_list bsl;
  two_particle_channel_list tpcl;
  std::vector<double> evals;
  std::vector< std::vector<double> > evecs;
  g_func g;
  v_func v;
  std::vector< std::vector<double> > elements;

  public:
    Hamiltonian(int n, double L, bare_state_list bsl, two_particle_channel_list tpcl, g_func g=nullptr, v_func v=nullptr);
    int Get_n();
    double Get_L();
    bare_state_list Get_BareStateList();
    two_particle_channel_list Get_TwoParticleChannelList();
    std::vector<double> Get_Eigenvector(int i);
    double Get_Eigenvalue(int i);
    std::vector< std::vector<double> > Get_Elements();
    void Set_n(int n);
    void Set_L(double L);
    void Set_BareStates(bare_state_list bsl);
    void Set_TwoParticleChannels(two_particle_channel_list tpcl);
    void Set_g(g_func g); // can pass in a function g as g or &g
    void Set_v(v_func v); // can pass in a function v as v or &v
    double Eval_g(double k, bare_state bs, two_particle_channel tpc);
    double Eval_v(double k1, double k2, two_particle_channel tpc1, two_particle_channel tpc2);
    int Get_Dimension(); // returns the number of rows (= number of columns) of Hamiltonian object
    std::vector<double> operator [](int i); // returns a row of the Hamiltonian matrix-- ensures that the user can use H[][]
    double Get_MatrixElement(int i, int j); // returns entry in row i, column j (starting index from 0, naturally)
    void Print_SquareBlock(int i); // prints (to std out) out the upper left square block of matrix, up to row i and column i
    void Sort_Eigenvalues(); // sorts eigenvalues by least to greatest, and also changes the order of vectors in evecs to correspond to new evals list
    void Print_Eigen(); // prints out eigenvalues and their associated eigenvectors
    void Print_Eigenvalues(); // prints out just the list of eigenvalues (eigenvalues will be sorted, but this won't change the order of evals above)
};

//////////////////////////////// User Functions

double default_g_potential(double k, bare_state bs, two_particle_channel tpc);

double default_v_potential(double k1, double k2, two_particle_channel tpc1, two_particle_channel tpc2);

/////////////////////////////////////////////////////////


#endif
