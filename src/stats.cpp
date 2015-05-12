#include <vector>
#include <cmath>
#include <unordered_map>
#include <string>
#include <utility>
#include <iostream>

#include "stats.hpp"
#include "typedefs.hpp"
#include "utilities.hpp"

using namespace std;

//Let's assume it wont't overflow
double computeMean(const vector<double> &vec){
    double sum = accumulate(vec.cbegin(), vec.cend(), 0.0);
    return sum/(double)vec.size();
}

double computeStandardDeviation(const vector<double> &vec, bool correction){
    double m = computeMean(vec);
    double accum = 0.0;
    for_each(vec.cbegin(), vec.cend(), [&](const double d){
    	accum += (d-m)*(d-m);
    });
    if(correction){
    	return sqrt(accum/(vec.size()-1));
    }
    else{
    	return sqrt(accum/vec.size());
    }
}

double computeZeroPercentage(const vector<double> &vec){
    return (double)count(vec.cbegin(), vec.cend(), 0.0)/(double)vec.size();
}

double computePearsonCorrelation(const vector<double> &x, const vector<double> &y, double x_mean, double x_stddev, double y_mean, double y_stddev){
    double accum = 0.0;
    for_each_two_ranges(x.cbegin(), x.cend(), y.cbegin(), [&accum, x_mean, y_mean](const double a, const double b){
    	accum += (a-x_mean)*(b-y_mean);
    });
    return accum/(x.size()*x_stddev*y_stddev);
}

double computePearsonCorrelation(const vector<double> &x, const vector<double> &y){
    double accum = 0.0;
    double x_mean = computeMean(x);
    double y_mean = computeMean(y);
    for_each_two_ranges(x.cbegin(), x.cend(), y.cbegin(), [&accum, x_mean, y_mean](const double a, const double b){
    	accum += (a-x_mean)*(b-y_mean);
    });
    double x_stddev = computeStandardDeviation(x, false);
    double y_stddev = computeStandardDeviation(y, false);
    return accum/(x.size()*x_stddev*y_stddev);
}

vector<double> computePearsonCorrelation(const vector<vector<double>> &M){
	unsigned int N = M.size();
	cout << endl << "Pearson correlation to be computed for " << N << " vectors..." << endl;
	vector<double> correlationMatrix(N*N);
	vector<double> means(N);
	vector<double> standard_deviations(N);

	cout << "Computing all means and standard deviations... " << flush;
	transform(M.cbegin(), M.cend(), means.begin(), computeMean);
	transform(M.cbegin(), M.cend(), standard_deviations.begin(), [](const vector<double> &vec){
		return computeStandardDeviation(vec, false);
	});
	cout << "Done."  << endl;

	unsigned int count = 0;

	cout << "Computing correlations for all pairs of vectors... " << endl;
	for(unsigned int i=0; i<N; ++i){
		printAdvancement(count, (N*(N-1)/2));
		correlationMatrix[i+N*i] = 1.0;
		for(unsigned int j=i+1; j<N; ++j){
			double cor = computePearsonCorrelation(M[i], M[j], means[i], standard_deviations[i], means[j], standard_deviations[j]);
			correlationMatrix[i+N*j] = cor;
			correlationMatrix[j+N*i] = cor;
		}
		count += N-i+1;
	}
	cout << "Done. " << endl << flush;

	return correlationMatrix;
}

//Assume x is sorted
void computeRankSorted(vector<double> &x){
	unsigned int current = 0;
	while(current != x.size()){
		double d = x[current];
		unsigned int next = current+1;
		int sum = current;
		int count = 1;
		while(next < x.size() && d == x[next]){
			++count;
			sum += next;
			++next;
		}

		for(unsigned int j=current; j<next; ++j){
			x[j] = (double)sum/(double)count;
		}

		current = next;
	}
}

void computeRank(vector<double> &x){
	vector<double> copyX(x);
	vector<size_t> sortedXIndexes = get_rank_increasing(copyX);
	sort(copyX.begin(), copyX.end());
	computeRankSorted(copyX);
	for(unsigned int i=0;  i != x.size(); ++i){
		x[i] = copyX[sortedXIndexes[i]];
	}
}

double computeSpearmanCorrelation(const vector<double> &x, const vector<double> &y){
	vector<double> copyX(x);
	vector<double> copyY(y);
	computeRank(copyX);
	computeRank(copyY);
	return computePearsonCorrelation(copyX, copyY);
}

vector<double> computeSpearmanCorrelation(const vector<vector<double>> &M){
	vector<vector<double>> M_copy(M);
	for_each(M_copy.begin(), M_copy.end(), computeRank);
	return computePearsonCorrelation(M_copy);
}

unordered_map<string, vector<pair<double, double>>> computeControlDistribution(const RNASeqData &controlData){

    unordered_map<string, vector<pair<double, double>>> controlDistributionParameters;

    cout << endl << "*********** Computing mean and standard deviation for each gene ************" << endl;

    for(const auto &mappedData : controlData){
        string cancerName = mappedData.first;
        if(controlData.at(cancerName).at(0).empty()){
            continue;
        }
        controlDistributionParameters.insert(make_pair(cancerName, vector<pair<double, double>>()));
        for(const vector<double> &vec : mappedData.second){
            double mean = computeMean(vec);
            double stddev = computeStandardDeviation(vec);
            controlDistributionParameters[cancerName].push_back(make_pair(mean, stddev));
        }
    }

    return controlDistributionParameters;
}

void computeZScore(RNASeqData &tumorData, const unordered_map<string, vector<pair<double, double>>> &controlDistributionParameters){

    cout << endl << "********* Computing Z Scores ********** " << endl;

    for(auto &mappedData : tumorData){
        string cancerName = mappedData.first;
        if(controlDistributionParameters.find(cancerName) == controlDistributionParameters.end()){
            cout << "Problem with Z Score for " << cancerName << endl;
            return;
        }

        int i=0;
        for(vector<double> &vec : mappedData.second){
            for(double &d : vec){
                double stddev = controlDistributionParameters.at(cancerName).at(i).second;
                if(stddev>0){
                    d = (d-controlDistributionParameters.at(cancerName).at(i).first)/stddev;
                }
                else{
                    d = 0;
                }
            }
            ++i;
        }
    }
}
