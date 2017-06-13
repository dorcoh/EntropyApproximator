#include <sys/time.h>
#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <ilcp/cpext.h>
#include <vector>
#include <set>
#include <bitset>
#include <iterator>

// use ILOG's STL namespace
ILOSTLBEGIN
static const int DEFAULT_FILTER_THRESHOLD = 4;
static const int WHEN_GRAPH_FILTERING_IS_BETTER = 5;

IloInt parity_number          = 0;
IloInt parity_minlength       = -1;
IloInt parity_maxlength       = -1;
unsigned long parity_seed;
IloBool parity_use_given_seed = false;
IloInt parity_filterLevel     = 2; 
   // parity_filterLevel = 0: binary representation, individual xors
   // parity_filterLevel = 1: binary representation, Gaussian elimination
   // parity_filterLevel = 2: CP variable representation, individual xors
IloInt parity_filterThreshold = DEFAULT_FILTER_THRESHOLD; 

bool PARITY_DONT_HANDLE_RANDOM_SEED = false;


bool yannakis =true;
bool jaroslow = false;
bool wainr = false;
long short_xor_max_length  = 10;
bool use_pairwise_subs = false;

// global parameters, etc.
unsigned long seed;
bool          use_given_seed = false;
IloInt        timelimit      = -1;
bool          use_tb2        = false;
char          instanceName[1024];

unsigned long get_seed(void) {
  struct timeval tv;
  struct timezone tzp;
  gettimeofday(&tv,&tzp);
  return (( tv.tv_sec & 0177 ) * 1000000) + tv.tv_usec;
}

////////////////////////////////

void parseParityArgs(int & argc, char **argv)
{
  // this method eats up all arguments that are relevant for the
  // parity constraint, and returns the rest in argc, argv

  char residualArgv[argc][64];
  strcpy(residualArgv[0], argv[0]);
  int residualArgc = 1;

  for (int argIndex=1; argIndex < argc; ++argIndex) {
    if ( !strcmp(argv[argIndex], "-paritylevel") ) {
      argIndex++;
      parity_filterLevel = atol( argv[argIndex] );
    }
    else if ( !strcmp(argv[argIndex], "-paritythreshold") ) {
      argIndex++;
      parity_filterThreshold = atol( argv[argIndex] );
      if (parity_filterThreshold < 2) {
        cerr << "ERROR: paritythreshold must be at least 2." << endl;
        exit(1);
      }
    }
    else if ( !strcmp(argv[argIndex], "-number") ) {
      argIndex++;
      parity_number = (IloInt)atol( argv[argIndex] );
    }
    else if ( !strcmp(argv[argIndex], "-minlength") ) {
      argIndex++;
      parity_minlength = (IloInt)atol( argv[argIndex] );
    }
	   else if ( !strcmp(argv[argIndex], "-feldman") ) {
		argIndex++;
     short_xor_max_length = atol( argv[argIndex] );
		wainr = true;
	  }
	  	   else if ( !strcmp(argv[argIndex], "-jaroslow") ) {
		argIndex++;
     short_xor_max_length = atol( argv[argIndex] );
		jaroslow = true;
	  }
	  else if ( !strcmp(argv[argIndex], "-yannakis") ) {
		yannakis = true;
	  }
	  	  else if ( !strcmp(argv[argIndex], "-pairwisesubs") ) {
		use_pairwise_subs = true;
	  }
    else if ( !strcmp(argv[argIndex], "-maxlength") ) {
      argIndex++;
      parity_maxlength = (IloInt)atol( argv[argIndex] );
    }
    else if ( !strcmp(argv[argIndex], "-seed") && !PARITY_DONT_HANDLE_RANDOM_SEED ) {
      argIndex++;
      parity_seed =  atol( argv[argIndex] );
      parity_use_given_seed = true;
    }
    else {
      // save this option to be returned back
      strcpy(residualArgv[residualArgc++], argv[argIndex]);
    }
  }

  argc = residualArgc;
  for (int i=1; i<argc; ++i) {
    // free(argv[i]);
    argv[i] = new char[strlen(residualArgv[i])+1];
    strcpy(argv[i], residualArgv[i]);
  }
}

void printParityUsage(ostream & os = cout) {
  os << "ParityConstraint options:" << endl
     << "   -paritylevel        0: binary individual Xors," << endl 
     << "                       1: binary Gaussian elimination, " << endl 
     << "                       2: non-binary individual Xors (default)" << endl
     << "   -paritythreshold    >= 2, for individual Xors (default: 3)" << endl 
     << "   -number             Number of random XORs (default: 0)" << endl
     << "   -minlength          Minlength of XORs (default: nvars/2)" << endl
     << "   -maxlength          Maxlength of XORs (default: nvars/20)" << endl;
  if (!PARITY_DONT_HANDLE_RANDOM_SEED)
    os << "   -seed               Random seed" << endl;
  else {
    os << "   -seed               Disabled; to enable, remove from the model" << endl
       << "                         PARITY_DONT_HANDLE_RANDOM_SEED = true" << endl;
  }
  os << endl;
}


void parseArgs(int argc, char **argv) 
{
  // one argument must be the instance filename
  if (argc <= 1) {
    cerr << "ERROR: instance name must be specified" << endl
         << "       See usage (WH_cpo_uai -h)" << endl;
    exit(1);
  }

  for (int argIndex=1; argIndex < argc; ++argIndex) {
    if ( !strcmp(argv[argIndex], "-tb2") ) {
      use_tb2 = true;
    }
    else if ( !strcmp(argv[argIndex], "-timelimit") ) {
      argIndex++;
      timelimit = atol(argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-seed") ) {
      argIndex++;
      seed =  atol( argv[argIndex] );
      use_given_seed = true;
    }
    else if ( !strcmp(argv[argIndex], "-verbosity") ) {
      argIndex++;
    }
    else if ( !strcmp(argv[argIndex], "-h") || !strcmp(argv[argIndex], "-help") ) {
      cout << endl
           << "USAGE: iloglue_uai [options] instance.uai" << endl
           << endl
           << "   -timelimit          Timelimit in seconds (default None)" << endl
           << "   -seed               Random seed" << endl
           << endl;
      // print parity constraint options usage
      //printParityUsage(cout);
      exit(0);
    }
    else if (argv[argIndex][0] != '-') {
     // must be the instance name
     strcpy(instanceName, argv[argIndex]);
    }
    else {
      cerr << "ERROR: Unexpected option: " << argv[argIndex] << endl
           << "       See usage (iloglue_uai -h)" << endl;
      exit(1);
    }
  }
}

double median(vector<double> &v)
{
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}

void print_matrix (vector <vector <bool> > A)
{
for (unsigned int i =0;i<A.size();i++)
	{
		for (unsigned int j =0;j<A[i].size();j++)					// last column is for coefficients b
			cout << A[i][j] << ",";
		cout << endl;
	}
}

vector <vector <bool> > generate_matrix(int m, int n)
{
	vector <vector <bool> > A;
	A.resize(m);
	for (int i =0;i<m;i++)
	{
	A[i].resize(n+1);
	for (int j =0;j<n+1;j++)					// last column is for coefficients b
		//if (rnd_uniform()<0.5)
		if (rand()%2==0)
			A[i][j] = true;
		else
			A[i][j]=false;
	}
	
	// print
	//cout << "Random matrix" <<endl;
	//print_matrix(A);
	
	return A;	
}

vector <bool>  feasiblesol;

void row_echelon(vector <vector <bool> > & A)
{
	bool solvable = true;						// is the system A x + b =0 solvable?
	
	size_t m = A.size();
	size_t n = A[0].size()-1;
	
	vector <int> indep_columns;
		vector <int> indep_columns_rindex;
		set <int> indep_vars;
		
	// put A in row echelon form
	for (size_t i = 0;i<m;i++)
	{
   //Find pivot for column k:
   bool empty_row=true;
   int j_max =0;
	for (int s = 0;s<n;s++)
		if (A[i][s]) 
			{
			empty_row=false;
			j_max = s;
			break;
			}
		if (empty_row)				// low rank
			{
				if (A[i][n])		//0=1
				{
				solvable = false;
				}
				else
				{
				// 0 =0, we can remove this row
				}			
			}
		else
			{
				indep_vars.insert(j_max);
				indep_columns.push_back(j_max);					// index of a basis of A
				indep_columns_rindex.push_back(i);				// row index of pivot
				
				for (size_t h=i+1;h<m;h++)
					if (A[h][j_max])			// if not already zero
						{
							for (int q=0;q<n+1;q++)			// sum the two rows
								A[h][q] = A[h][q] xor A[i][q];
						}
			}

	}
	
	for (size_t i = 0;i<indep_columns.size();i++)
	{
	int j_max = indep_columns[i];
	int p = indep_columns_rindex[i];
	
					for (int h=p-1;h>=0;h--)
					if (A[h][j_max])			// if not already zero
						{
							//print_matrix(A);

							for (int q=0;q<n+1;q++)			// sum the two rows
								A[h][q] = A[h][q] xor A[p][q];
							//cout<<h << " = "<<h <<  "+"<<p << endl;
							//print_matrix(A);
						}
	
	
	}
	
	
	// produce a solution
	

	vector <bool>  b;
	vector <bool>  y;
	feasiblesol.resize(n);
	y.resize(n);
	
	// initialize b to the last column of A
	b.resize(m);
	for (size_t i =0;i<m;i++)
		b[i] = A[i][n];
	
	for (size_t i =0;i<n;i++)
		y[i] = rand()%2;
		
	// sum all the dependent variables that are already fixed
	for (size_t i =0;i<n;i++)	
		{
		feasiblesol[i] = y[i];
		if ( (indep_vars.count(i)==0) && (y[i]==1))		// dependent variable, and non zero
			{				
			// b = b + x[i] * A[] [i]
				for (size_t j =0;j<m;j++)
					b[j] = b[j] xor A[j][i];
			}
		}
		
	/*
	cout << "Printing b" << endl;
	for (size_t j =0;j<m;j++)
		cout << b[j] << ",";
	cout <<  endl;
	*/
	// backsubstitute r

	for (int i =indep_columns_rindex.size()-1;i>=0;i--)
		{
			int c = indep_columns_rindex[i];		// lowest pivot
			if (b[c]==1)		// we need to add a 1
				{
				y[indep_columns[i]] = 1;
				feasiblesol[indep_columns[i]] = 1;
				for (size_t j =0;j<m;j++)
					b[j] = b[j] xor A[j][indep_columns[i]];
				}
			else
			{
				y[indep_columns[i]] = 0;
				feasiblesol[indep_columns[i]] = 0;
			}
		}
	
	
}

//

int sparsify(vector <vector <bool> > & A)
{
vector< bitset<1000> > bv;

bv.resize(A.size());
size_t m = A.size();
size_t n = A[0].size();
int saved_bits = 0;
size_t initialnbr=0;

for (size_t i = 0;i<m;i++)
	{
	bitset<1000> row;
	for (size_t s = 0;s<n;s++)
		if (A[i][s]) 
			{
			row.set(s,1);
			initialnbr++;
			}
	bv[i]=row;	
	}

cout << "Initial # of bits: " << 	initialnbr << endl;


if (m<100)
{
for (size_t i = 0;i<m;i++)
	for (size_t l = 0;l<m;l++)
		for (size_t z = 0;z<m;z++)
		for (size_t g = 0;g<m;g++)
		if (i!=l && i!=z && z!=l && i!=g && z!=g && l!=g)
		{
		int cursize = bv[i].count();
		int newsize =  (bv[i]^bv[l]^bv[z]^bv[g]).count();
		if (newsize<cursize)
			{
			saved_bits = saved_bits + cursize - newsize;
			bv[i]^=bv[l]^bv[z]^bv[g];
			for (int q=0;q<n;q++)			// sum the two rows
				A[i][q] = A[i][q] xor A[l][q] xor A[z][q] xor A[g][q];
			}
		}

}
if (m<500)
{
for (size_t i = 0;i<m;i++)
	for (size_t l = 0;l<m;l++)
		for (size_t z = 0;z<m;z++)
		if (i!=l && i!=z && z!=l)
		{
		int cursize = bv[i].count();
		int newsize =  (bv[i]^bv[l]^bv[z]).count();
		if (newsize<cursize)
			{
			saved_bits = saved_bits + cursize - newsize;
			bv[i]^=bv[l]^bv[z];
			for (int q=0;q<n;q++)			// sum the two rows
				A[i][q] = A[i][q] xor A[l][q] xor A[z][q];
			}
		
		}
}
if (m<10000)
{		
for (size_t i = 0;i<m;i++)
	for (size_t l = 0;l<m;l++)
		if (i!=l)
		{
		int cursize = bv[i].count();
		int newsize =  (bv[i]^bv[l]).count();
		if (newsize<cursize)
			{
			saved_bits = saved_bits + cursize - newsize;
			bv[i]^=bv[l];
		//	print_matrix(A);
			
			for (int q=0;q<n;q++)			// sum the two rows
				A[i][q] = A[i][q] xor A[l][q];
			
			//cout<<i << " = "<<i <<  "+"<<l << endl;
		//	print_matrix(A);
			//cout<<endl;
			}
		
		}
}
/*
for (size_t att = 0;att<10000000;att++)	
{

bitset<1000> rndcomb;
for (size_t i = 0;i<m;i++)
	rndcomb[i] = rand()%2;
	
bitset<1000> res;
for (size_t i = 0;i<m;i++)
	if (rndcomb[i]>0)
		res^=bv[i];

for (size_t i = 0;i<m;i++)
	if (rndcomb[i]>0 && bv[i].count()>res.count())		
		{
		cout << "+";
		saved_bits = saved_bits +bv[i].count() - res.count();
		for (size_t l = 0;l<m;l++)
			if (rndcomb[l]>0 && i!=l)
				{
				for (int q=0;q<n;q++)
					A[i][q] = A[i][q] xor A[l][q];
				bv[i] ^= bv[l];
				}
		break;			// ??
		}
		
}		
*/
cout << "final # of bits: " << 	(int) initialnbr- saved_bits<< endl;		
return saved_bits;	
}

void add_linear_combinations(vector <vector <bool> > & A, size_t M)
{
//size_t cursize = A.size();
//size_t newrindex = 0;

for (size_t i =0;i<M;i++)
	for (size_t k =i;k<M;k++)
		if (k!=i)
			{
			vector <bool> newrow;
			newrow.resize(A[i].size());
			for (size_t j =0;j<A[i].size();j++)
			{
			newrow[j] = A[i][j] xor  A[k][j];
			}
			A.push_back(newrow);
			cout << "adding" << endl;
			}
}

vector <vector <bool> > generate_matrix_maxlength(int m, int n, int k)
{
	vector <vector <bool> > A;
	A.resize(m);
	
	vector <size_t> index;
	index.resize(n);
	for (int j =0;j<n;j++)
		index[j] = j;
			
	
	for (int i =0;i<m;i++)
	{
	A[i].resize(n+1);
	std::random_shuffle(index.begin(), index.end());
	for (int j =0;j<k;j++)					// last column is for coefficients b
		//if (rnd_uniform()<0.5)
		A[i][index[j]] = true;
	}
	
	// fill parity bits at random
	for (int i =0;i<m;i++)
		if (rand()%2==0)
			A[i][n] = true;
		else
			A[i][n] = false;
	// print
	//cout << "Random matrix" <<endl;
	//print_matrix(A);
	
	return A;	
}

vector <vector <bool> > generate_Toeplitz_matrix(int m, int n)
{
	
	vector <vector <bool> > A;
	if (m==0)
		return A;
		
	A.resize(m);
	int i;
	for (i =0;i<m;i++)
	{
	A[i].resize(n+1);
	}
	
	// first column
	for (i =0;i<m;i++)
	{
		if (rand()%2==0)
			A[i][0] = true;
		else
			A[i][0]=false;
		for (int j =1;j<m-i;j++)
			if (j<n)
				A[i+j][j] = A[i][0];
	}
	
		// last column
	for (i =0;i<m;i++)
	{
		if (rand()%2==0)
			A[i][n] = true;
		else
			A[i][n]=false;

	}


	
	// first row
	for (int j =1;j<n;j++)
	{
		if (rand()%2==0)
			A[0][j] = true;
		else
			A[0][j]=false;
			
		for (i =1;i<m;i++)
			if (j+i<n)
				A[i][j+i] = A[0][j];
	}
	
	// print
	//cout << "Random matrix3" <<endl;
	//print_matrix(A);
	
	return A;	
}

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;
 
powerset_type powerset2(set_type const& set)
{
  typedef set_type::const_iterator set_iter;
  typedef std::vector<set_iter> vec;
  typedef vec::iterator vec_iter;
 
  struct local
  {
    static int dereference(set_iter v) { return *v; }
  };
 
  powerset_type result;
 
  vec elements;
  do
  {
    set_type tmp;
    std::transform(elements.begin(), elements.end(),
                   std::inserter(tmp, tmp.end()),
                   local::dereference);
    result.insert(tmp);
    if (!elements.empty() && ++elements.back() == set.end())
    {
      elements.pop_back();
    }
    else
    {
      set_iter iter;
      if (elements.empty())
      {
        iter = set.begin();
      }
      else
      {
        iter = elements.back();
        ++iter;
      }
      for (; iter != set.end(); ++iter)
      {
        elements.push_back(iter);
      }
    }
  } while (!elements.empty());
 
  return result;
}

set <set <int> > powerset (set <int> s) {
   set <set <int> > result;
   set <int> nullset; //  the default constructor builds a set with no elements

   /*
   	set <int> ::iterator it4;
		for (it4 = s.begin ( ); it4 != s.end (); it4++)
			cout << (*it4) << " ";
	cout << endl;
	*/		
   if (s.size( ) == 0) {
      result.insert (nullset);
      return result;
   }
   
 //  if (s.size() == k) { result.insert(s); return result; }
   else {
      set <int>::iterator it;
      for (it = s.begin(); it != s.end(); it++) {
         int elem = *it;

         //  copy the original set, and delete one element from it.
         set <int> s1 (s);
         s1.erase (elem);

         //  compute the power set of this smaller set.
         set <set <int> > p1 = powerset (s1);
        
		

		set <set <int> >::iterator it3;
		for (it3 = p1.begin ( ); it3 != p1.end (); it3++) {
		result.insert (*it3);
		}

         //  add the deleted element to each member of this power set,
         //  and insert each new set into the desired result.
         set <set <int> >::iterator iter;
         for (iter = p1.begin(); iter != p1.end(); iter++) {
            set <int>  next = *iter;
			
				next.insert (elem);

				result.insert (next);
				
         };
      };
      return result;
   }
}


vector <vector <bool> > substitute_pairwise_vars(vector <vector <bool> > & A	,  std::vector < std::vector< std::vector < std::vector< IloBoolVar> > > > Mu)
{
vector <vector <bool> > B;
size_t m = A.size();
B.resize(m);

// copy A
for (size_t i = 0;i<m;i++)
	{
	B[i].resize(A[i].size()+Mu.size()*Mu.size());
	for (size_t s = 0;s<A[i].size();s++)
		B[i][s]= A[i][s];
	}

// for each entry, check all other entries. 
for (size_t i = 0;i<m;i++)
	{
	for (size_t s = 0;s<A[i].size()-1;s++)
		for (size_t s2 = 0;s2<A[i].size()-1;s2++)
			if (B[i][s] && B[i][s2])
			{
			if (!Mu[s][s2].empty())		// we have the pairwise var already
				{
				//cout << s << "," << s2 << " -->" << A[i].size()+s*Mu.size()+s2 << endl;
				B[i][s] = false;
				B[i][s2] = false;
				B[i][A[i].size()+s*Mu.size()+s2] = true;
				// add the pairwise var
				}
			}
	}

return B;	
}

IloCP IlogSolver;

// Usage: iloglue problem_name.wcsp [verbosity]
int main(int argc, char **argv)
{
  char pbname[1024];
  int nbvar,nbval,nbconstr;
  IloEnv env;
  IloTimer timer(env);
    
  try {
    // first parse and remove parity-related command-line arguments
    PARITY_DONT_HANDLE_RANDOM_SEED = true;
    parseParityArgs(argc, argv);

    // now parse regular arguments
    parseArgs(argc, argv);

    // associate the CP solver with the environment
    IlogSolver = IloCP(env);

    // associate a model with the environment
    IloModel model(env);
        
    // open the instance file
    ifstream file(instanceName);
    if (!file) {
      cerr << "Could not open file " << instanceName << endl;
      exit(EXIT_FAILURE);
    }
  
    // stefano mod, read uai file
    // reads uai file to parse domain sizes; creates variables along the way
    cerr << "Creating variables"<< endl;
    file >> pbname;
    file >> nbvar;
    IloIntVarArray vars(env, nbvar, 0, 100);
    nbval = 0;
    int tmp;
    for (int i=0; i<nbvar; i++) {
      file >> tmp;
      if (tmp>nbval)
        nbval = tmp;
      vars[i].setBounds(0, tmp-1);									// (17)
      char *name = new char[16];
      sprintf(name, "x%d", i);
      vars[i].setName(name);
    }
    model.add(vars);
    file >> nbconstr;
    cerr << "Var:"<< nbvar <<" max dom size:" <<nbval<<" constraints:"<<nbconstr << endl;

    // define variable that captures the value of the objective function
    IloIntVar obj(env, 0, IloIntMax, "objective");


      // use a native CP Optimizer representation of the .uai file
      // read in variable scopes of CPT tables
      std::vector < std::vector< int> > scopes;
      int arity;
      for (int i=0; i<nbconstr; i++) {
        file >> arity;
		std::vector< int> scope;
		int id;
        for (int j=0; j<arity; j++) {
          file >> id;
          scope.push_back(id);
        }
		scopes.push_back(scope);
      }
      // read in values of CPT tables
      IloInt l;
      IloArray<IloNumArray> cost(env);
      int TableSize;
      for (int i=0; i<nbconstr; i++) {
        file >> TableSize;
        double prod;
        IloNum entry;
        IloNumArray table(env);
        for (int j=0; j<TableSize; j++) {
          file >> prod;				
          entry = log10(prod);		// integer vs real
          //entry = 10;		// integer vs real
        //  cout << "adding " << entry << endl;
          table.add(entry);
        }
        cost.add(table);
      }
      cout << "done reading CPTs"<< endl;
      // define cost expression
      IloNumExpr objexpr(env);
	  std::vector < std::vector< std::vector < std::vector< IloBoolVar> > > > Mu;
	  
	  Mu.resize(nbvar);
	  for (size_t q= 0; q<nbvar;q++)
			Mu[q].resize(nbvar);
			
      for (l = 0; l < nbconstr; l++) {        
        IloIntExpr pos(env);			// init to 0
	
	
	
	if (scopes[l].size()==1)
	{
	//objexpr += cost[l][0]* vars[scopes[l][0]]+cost[l][1]* (1-vars[scopes[l][0]]);
	
	if (isfinite(cost[l][0]))
		objexpr += cost[l][0]* (1-vars[scopes[l][0]]);
	else
		{
		if (isinf(cost[l][0]))
			{
			model.add(vars[scopes[l][0]]==1);
			}
		else
			{
			cout << "Cannot generate ILP"<< endl;
			exit(-1);
			}
		}
	
	if (isfinite(cost[l][1]))
		objexpr += cost[l][1]* (vars[scopes[l][0]]);
	else
		{
		if (isinf(cost[l][1]))
			{
			model.add(vars[scopes[l][0]]==0);
			}
		else
			{
			cout << "Cannot generate ILP"<< endl;
			exit(-1);
			}
		}
	
	}
	else
	{
			int i = scopes[l][0];
			int j = scopes[l][1];
			
			char *name = new char[32];
			sprintf(name, "mu_%d_%d (0,0)", (int) i, (int) j);
			IloBoolVar mu_i_j_0_0 (env, 0, 1, name);					// (18)
			model.add(mu_i_j_0_0);
			
			sprintf(name, "mu_%d_%d (0,1)", (int) i, (int) j);
			IloBoolVar mu_i_j_0_1 (env, 0, 1, name);					// (18)
			model.add(mu_i_j_0_1);
			
			sprintf(name, "mu_%d_%d (1,0)", (int) i, (int) j);
			IloBoolVar mu_i_j_1_0 (env, 0, 1, name);
			model.add(mu_i_j_1_0);
			
			sprintf(name, "mu_%d_%d (1,1)", (int) i, (int) j);
			IloBoolVar mu_i_j_1_1 (env, 0, 1, name);
			model.add(mu_i_j_1_1);
			
			model.add((mu_i_j_0_0+mu_i_j_1_0 == 1-vars[j]));
			//model.add((mu_i_j_0_0+mu_i_j_1_0 <= vars[j]));
			
			model.add((mu_i_j_0_1+mu_i_j_1_1 == vars[j]));
			//model.add((mu_i_j_0_1+mu_i_j_1_1 <= 1-vars[j]));
			
			model.add((mu_i_j_0_0+mu_i_j_0_1 == 1-vars[i]));
			//model.add((mu_i_j_0_0+mu_i_j_0_1 <= 1-vars[i]));
			
			model.add((mu_i_j_1_0+mu_i_j_1_1 == vars[i]));
			//model.add((mu_i_j_1_0+mu_i_j_1_1 <= vars[i]));
			
			
			model.add((mu_i_j_0_1+mu_i_j_1_1 <= 1));
			model.add((mu_i_j_0_0+mu_i_j_1_0 <= 1));
			model.add((mu_i_j_1_0+mu_i_j_1_1 <= 1));
			model.add((mu_i_j_0_0+mu_i_j_0_1 <= 1));
			
			
			Mu[i][j].resize(2);
			Mu[i][j][0].resize(2);
			Mu[i][j][1].resize(2);
			Mu[i][j][0][0]= mu_i_j_0_0 ;
			Mu[i][j][0][1]= mu_i_j_0_1 ;
			Mu[i][j][1][0]= mu_i_j_1_0 ;
			Mu[i][j][1][1]= mu_i_j_1_1 ;
			
		
		
		//objexpr += cost[l][0]* mu_i_j_0_0;
		
		if (isfinite(cost[l][0]))
			objexpr += cost[l][0]* mu_i_j_0_0;
		else
		{
		if (isinf(cost[l][0]))
			{
			model.add(mu_i_j_0_0==0);
			}
		else
			{
			cout << "Cannot generate ILP"<< endl;
			exit(-1);
			}
		}	
		
		//objexpr += cost[l][1]* mu_i_j_0_1;
		
		if (isfinite(cost[l][1]))
			objexpr += cost[l][1]* mu_i_j_0_1;
		else
		{
		if (isinf(cost[l][1]))
			{
			model.add(mu_i_j_0_1==0);
			}
		else
			{
			cout << "Cannot generate ILP"<< endl;
			exit(-1);
			}
		}
		
		
		//objexpr +=cost[l][2]* mu_i_j_1_0;
		
		if (isfinite(cost[l][2]))
			objexpr +=cost[l][2]* mu_i_j_1_0;
		else
		{
		if (isinf(cost[l][2]))
			{
			model.add(mu_i_j_1_0==0);
			}
		else
			{
			cout << "Cannot generate ILP"<< endl;
			exit(-1);
			}
		}
		
		//objexpr += cost[l][3]* mu_i_j_1_1;	

		if (isfinite(cost[l][3]))
			objexpr += cost[l][3]* mu_i_j_1_1;	
		else
		{
		if (isinf(cost[l][3]))
			{
			model.add(mu_i_j_1_1==0);
			}
		else
			{
			cout << "Cannot generate ILP"<< endl;
			exit(-1);
			}
		}

		
				
	//	objexpr += cost[l][0]* ((vars[scopes[l][0]]<=0)&&(vars[scopes[l][1]]<=0))+ cost[l][1]* ((vars[scopes[l][0]]<=0)&&(vars[scopes[l][1]]>=1))+cost[l][2]* ((vars[scopes[l][0]]>=1)&&(vars[scopes[l][1]]<=0))+ cost[l][3]* ((vars[scopes[l][0]]>=1)&&(vars[scopes[l][1]]>=1));
	}
     //   for (unsigned j=0; j<scopes[l].size(); j++) {
     //     pos += ((int) pow(2.0,(double) scopes[l].size()-1-j))*vars[scopes[l][j]];
	//	}
       // objexpr +=cost[l][pos];
      }
      // make obj, the variable capturing the objective, equal this cost expression
    //  model.add(obj == objexpr);
	  

model.add(IloMaximize(env, objexpr ));

if (!use_given_seed) seed = get_seed();
srand(seed);

// generate matrix of coefficients A x = b. b is the last column

vector <vector <bool> > A = generate_Toeplitz_matrix(parity_number, nbvar);
cout << "here" << endl;
//vector <vector <bool> > A = generate_matrix_maxlength(parity_number, nbvar,parity_maxlength );


if (!A.empty())
{
	//print_matrix(A);
	cout << "row echelon form:" << endl;
	row_echelon(A);
		
	print_matrix(A);
	
	cout << "Bits saved: ";
	for (int i=0; i<2; i++)
		cout << sparsify(A) << " ";
	cout << endl;	

	//add_linear_combinations(A,5);
			

/*
for (int i=0; i<2; i++)
	{
	row_echelon(A);
	sparsify(A) ;
	sparsify(A) ;
	}
*/	




if (use_pairwise_subs)
{
	vector <vector <bool> > B  = substitute_pairwise_vars(A,Mu);
	//cout<<"after pairwise subtitution: "<< endl;
	//print_matrix(B);
	A=B;
	//cout << sparsify(A) << " " << endl;
}
}

std::vector < std::set <size_t> > varAppearancesInXors;
//varAppearancesInXors.resize(nbvar+1);						// dummy parity var

if (!A.empty())
varAppearancesInXors.resize(A[0].size());						// dummy parity var

IloArray<IloArray<IloIntVarArray> > zeta_vars(env);

//IloArray<IloArray<IloNumVarArray> > zeta_vars(env);

IloArray <IloIntVarArray> alpha_vars(env);

//IloArray<IloNumVarArray> alpha_vars(env);

IloBoolVar dummy_parity (env, 0, 1, "dummy");					// (18)
model.add(dummy_parity);
model.add((dummy_parity==1));
vector <size_t> xors_length;


	
if (!A.empty())
{

	xors_length.resize(A.size());
	
	for (size_t j = 0; j<A.size();j++)
		{
		size_t f = 0;
		
		for (size_t l = 0; l<A[j].size();l++)			// last column is the parity bit b
			if (A[j][l])
				{
				f++;
				}
		xors_length[j] = f;		// save length of j-th xor

		}
	
		for (size_t j = 0; j<A.size();j++)
		{
		if (!( (wainr || jaroslow) && xors_length[j]<=short_xor_max_length ))				// use yannakis encoding for longer ones
			for (size_t l = 0; l<A[j].size();l++)			// last column is the parity bit b
				if (A[j][l])
					{
					varAppearancesInXors[l].insert(j);		// for each var, save list of xors involved	
					}
		}
	
	cout << "XOR minimum length: " << *std::min_element(xors_length.begin(),xors_length.end()) <<" . XOR maximum length: " << *std::max_element(xors_length.begin(),xors_length.end()) << endl;
	
	if (yannakis)
	{
		
		for (size_t j= 0; j<A.size();j++)
			{
			IloIntVarArray alphas(env);
			if (!( (wainr || jaroslow) && xors_length[j]<=short_xor_max_length ))				// use yannakis encoding for longer ones
			{			

			
	//		IloNumVarArray alphas(env);
			
			// compute xor length
			size_t f =xors_length[j];					
			
			// add alpha_j_k var
			IloNumExpr alpha_sum_to_one(env);
			for (size_t k = 0; k<= 2* (size_t) floor(f/2);k=k+2)
				{
				char *name = new char[32];
				sprintf(name, "alpha_%d_%d", (int) j, (int) k);
				
				IloBoolVar alpha_j_k (env, 0, 1, name);					// (18)
	//			IloNumVar alpha_j_k (env, 0.0, 1.0, ILOFLOAT, name);					// (18)
				
				alphas.add(alpha_j_k);
				alpha_sum_to_one = alpha_sum_to_one + alpha_j_k;
				}

			//alpha_vars[j] = alphas;
			
			model.add(alphas);
			model.add((alpha_sum_to_one==1));							// (15)
			
			}
			alpha_vars.add(alphas);
			}

		// add zeta_i_j_k var nbvar
		
		for (size_t i= 0; i<A[0].size();i++)
			{
				IloArray<IloIntVarArray> zet_jk (env,A.size());
	//			IloArray<IloNumVarArray> zet_jk (env,A.size());
					
				std::set< size_t > XorsContainingi =  varAppearancesInXors[i];
				for (std::set<size_t >::iterator it=XorsContainingi.begin(); it!=XorsContainingi.end(); ++it)
					{
						size_t f = xors_length[*it];

						IloIntVarArray zet(env);
	//					IloNumVarArray zet(env);
						
						IloNumExpr zeta_sum_to_f(env);
						
						for (size_t k = 0; k<= 2* (size_t) floor(f/2);k=k+2)
						{

						char *name = new char[32];
						sprintf(name, "zeta_%d_%d_%d", (int) i, (int) *it, (int) k);
						
						IloBoolVar zeta_i_j_k (env, 0, 1, name);
						
						//IloNumVar zeta_i_j_k (env, 0, 1, ILOFLOAT,name);
						
						zet.add(zeta_i_j_k);
						zeta_sum_to_f = zeta_sum_to_f + zeta_i_j_k;
						
						model.add((zeta_i_j_k<=alpha_vars[*it][k/2]));	// (19)
						}
						model.add(zet);
						
						// if (i==nbvar)
							// model.add((zeta_sum_to_f==dummy_parity));			// (14)
						// else
							// model.add((zeta_sum_to_f==vars[i]));			// (14)
							
						if (i<nbvar)
							model.add((zeta_sum_to_f==vars[i]));			// (14)							
						else
							if (i==nbvar)															// dummy
								model.add((zeta_sum_to_f==dummy_parity));			// (14)
							else			// it's a pairwise
								{
								int i1 =(i - nbvar-1)/nbvar;
								int j1 =(i - nbvar-1)%nbvar;
								//cout << i << " -->" << i1 << " " << j1 << endl;
								model.add((zeta_sum_to_f==Mu[i1][j1][0][1]+Mu[i1][j1][1][0]));			// (14)
								}
							
						zet_jk[*it]=zet;
					}
				
			zeta_vars.add(zet_jk);	
			}
		
		for (size_t j= 0; j<A.size();j++)
			if (!( (wainr || jaroslow) && xors_length[j]<=short_xor_max_length ))				// use yannakis encoding for longer ones
			{	
				size_t f = xors_length[j];
				for (size_t k = 0; k<= 2* (size_t) floor(f/2);k=k+2)
					{
						IloNumExpr zeta_sum_to_alpha_k(env);
						for (size_t l = 0; l<A[j].size();l++)			// last column is the parity bit b
							{
							if (A[j][l])
								zeta_sum_to_alpha_k = zeta_sum_to_alpha_k + zeta_vars[l][j][k/2];
							}
						model.add((zeta_sum_to_alpha_k==IloInt(k)*alpha_vars[j][k/2]));			// (16)
		
					}
			}
	}
	
}



// jaroslaw encoding and wainwright for short xors
for (size_t j= 0; j<A.size();j++)
	{
	size_t f =  xors_length[j];
	if (f<=short_xor_max_length)
		{
		set <int> variables_involved;
		vector <IloNumExpr> sum_over_S_fori;
		//sum_over_S_fori.resize(nbvar+1);
		sum_over_S_fori.resize(A[0].size());

		for (size_t l= 0; l<A[j].size();l++)		// also the parity bit, dummy var
		{
			if (A[j][l])
				{
				variables_involved.insert(l);
				}
			sum_over_S_fori[l]=IloNumExpr(env);	
		}
		
		//cout << endl << variables_involved.size() << endl;
		// generate power set, then only use odd size
		set <set <int> > subsets = powerset2 (variables_involved);
		
		//for (size_t k= 1; k<=f;k=k+2)
			//{

		
		set <set <int> >::iterator it3;
				IloNumExpr sum_w_over_S(env);
			//for (it3 = subset_size_k.begin ( ); it3 != subset_size_k.end (); it3++) 
			for (it3 = subsets.begin ( ); it3 != subsets.end (); it3++) 
				{
				IloNumExpr fi_par_check_sum(env);
				set <int> S = (*it3);
				

				int counter = 0;
				
				if (jaroslow)
				{
					if (S.size()%2>0)
					{
						set <int> ::iterator it4;
						for (it4 = S.begin ( ); it4 != S.end (); it4++)
							{
							int l = (*it4);
							//cout << (*it4) << " ";
							// if (l==nbvar)
								// fi_par_check_sum = fi_par_check_sum +dummy_parity;
							// else
								// fi_par_check_sum = fi_par_check_sum + vars[l];
								
							if (l<nbvar)
								fi_par_check_sum = fi_par_check_sum + vars[l];		
							else
								if (l==nbvar)															// dummy
								fi_par_check_sum = fi_par_check_sum +dummy_parity;
								else			// it's a pairwise
								{
								int i1 =(l - nbvar-1)/nbvar;
								int j1 =(l - nbvar-1)%nbvar;
								//cout << l << " -->" << i1 << " " << j1 << endl;
								fi_par_check_sum = fi_par_check_sum + Mu[i1][j1][0][1]+Mu[i1][j1][1][0];		
								}
								
							}
						set <int> NminusS;
						std::set<int> s1, s2;
						//cout << " --- ";
						std::set_difference(variables_involved.begin(), variables_involved.end(), S.begin(), S.end(),
							std::inserter(NminusS, NminusS.end()));
						//it4 = set_difference(variables_involved.begin(),variables_involved.end(),S.begin(),S.end(),std::inserter(NminusS, NminusS.end()));
						for (it4 = NminusS.begin ( ); it4 != NminusS.end (); it4++)
							{
							int l = (*it4);
							//cout << (*it4) << " ";
							// if (l==nbvar)
								// fi_par_check_sum = fi_par_check_sum +(1-dummy_parity);
							// else
								// fi_par_check_sum = fi_par_check_sum + (1-vars[l]);
								
							if (l<nbvar)
								fi_par_check_sum = fi_par_check_sum + (1-vars[l]);		
							else
								if (l==nbvar)															// dummy
								fi_par_check_sum = fi_par_check_sum +(1-dummy_parity);
								else			// it's a pairwise
								{
								int i1 =(l - nbvar-1)/nbvar;
								int j1 =(l - nbvar-1)%nbvar;
								//cout << l << " -->" << i1 << " " << j1 << endl;
								fi_par_check_sum = fi_par_check_sum +(1- Mu[i1][j1][0][1]-Mu[i1][j1][1][0]);		
								}								
								
							}
						model.add((fi_par_check_sum<=IloInt(f-1)));	
					}
					//cout << endl;
				}
				if (wainr)
					{
						if (S.size()%2==0)
						{
							char *name = new char[32];
							sprintf(name, "w_%d_%d", (int) j, counter++);					
							IloBoolVar w_j_S(env, 0, 1, name);				// (5)
							
							sum_w_over_S = sum_w_over_S + w_j_S;
							
							set <int> ::iterator it4;
							for (it4 = S.begin ( ); it4 != S.end (); it4++)
								{
								int l = (*it4);
								sum_over_S_fori[l] = sum_over_S_fori[l] + w_j_S;		// S contains variable l, as in constraint (7)
								}
								
						}	
					}
				}
				
			if (wainr)
				{
				set <int> ::iterator it4;
				for (it4 = variables_involved.begin ( ); it4 != variables_involved.end (); it4++)
					{
					int l = (*it4);
					if (l==nbvar)
						model.add((sum_over_S_fori[l]==dummy_parity));
					else
						model.add((sum_over_S_fori[l]==vars[l]));
					}
									
				model.add((sum_w_over_S==1));					// (6)
				}
			//}
		}
	}
	


IloCplex cplex(model);

	if (timelimit > 0)
		cplex.setParam(IloCplex::TiLim, timelimit);
	/*
	IloNumArray    ordpri(env);
	for (size_t j= 0; j<nbvar;j++)
		ordpri.add(10.0);
	cplex.setPriorities(vars,ordpri);
	*/
	
	//IlogSolver.setParameter(IloCP::LogPeriod, 1000000);
	//IlogSolver.setParameter(IloCP::LogPeriod, 1);   // for debugging
	cplex.setParam(IloCplex::Threads, 1);    // number of parallel threads

//	cplex.setParam(IloCplex::Threads, 4);    // number of parallel threads
 //    cplex.setParam(IloCplex::ParallelMode, -1);
		
//	cplex.setParam(IloCplex::Cliques, IloInt (2));
	
	//cplex.setParam(IloCplex::MIPDisplay, 5);
	//cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);
	//cplex.setParam(IloCplex::NodeAlg, IloCplex::Dual);
	
	//cplex.setParam(IloCplex::MIPEmphasis,2); //CPX_MIPEMPHASIS_BESTBOUND

	//cplex.addMIPStart(vars, new double[] {5.0, 3.0});
	
	if (!A.empty())
	{
	IloNumArray feasibleinit(env);
	//double [] feasibleinit;
	IloNumVarArray startVar(env);
	
	for (size_t l= 0; l<nbvar;l++)
		{
		startVar.add(vars[l]);
		feasibleinit.add(feasiblesol[l]);
		}
	cplex.addMIPStart(startVar, feasibleinit);
	}
	
      if ( !cplex.solve() ) {
         env.out() << "Failed to optimize LP." << endl;
         throw(-1);
      }
//cout << objexpr;

     IloNumArray vals(env);
      env.out() << "Solution status = " << cplex.getStatus() << endl;
     // env.out() << "Solution value = " << cplex.getObjValue() << endl;
      env.out() << "Solution value log10lik = " << cplex.getObjValue() << endl;

	cplex.getValues(vals, vars);
      env.out() << "Values = " << vals << endl;
	
    } catch (IloException& ex) {
    cout << "Error: " << ex << endl;

  }

	    
   env.end();
  return 0;
}

