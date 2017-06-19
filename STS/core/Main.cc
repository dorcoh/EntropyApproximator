#include <math.h>
		
		// a vector for counting #(x) and #(!x)
		std::vector<int> varSolCountPos; // counter for pos
		std::vector<int> varSolCountNeg; // counter for neg
		std::vector<double> rvPos;		 // vector for r(v)
		std::vector<double> rvNeg;
		std::vector<double> ev;			 // vector for e(v)
		std::vector<double> entropy;	 // vector for formula entropy for each sample


		// init
		varSolCountNeg.resize(var_num);
		varSolCountPos.resize(var_num);
		rvPos.resize(var_num);
		rvNeg.resize(var_num);
		ev.resize(var_num);
		for (int iter=0; iter<var_num; iter++)
		{
			varSolCountPos[iter] = 0;
			varSolCountNeg[iter] = 0;
			rvPos[iter] = 0;
			rvNeg[iter] = 0;
			ev[iter] = 0;
		}

		entropy.resize(nsamples);

						// compute #(x) and #(!x)
						if (OutputSamples[l][i] == 1)
						{
							varSolCountPos[i] += 1;
						} else {
							varSolCountNeg[i] += 1;
						}
