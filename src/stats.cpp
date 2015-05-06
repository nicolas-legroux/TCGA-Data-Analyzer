#include <vector>
#include "stats.hpp"
#include <cmath>
#include <unordered_map>
#include <string>
#include <utility>
#include <iostream>

using namespace std;

typedef unordered_map<string, vector<vector<double>>> RNASeqData;

double computeMean(const vector<double> &vec){
    double mean = 0.0;
    double n = vec.size();
    for(double d : vec){
        mean += d/n;
    }
    return mean;
}

double computeStandardDeviation(const vector<double> &vec){
    vector<double> vecXvec;
    for(double d : vec){
        vecXvec.push_back(d*d);
    }
    double mean = computeMean(vec);
    return sqrt(computeMean(vecXvec)-mean*mean);
}

double computeZeroPercentage(const vector<double> &vec){
    int countZero = 0;
    for(double d : vec){
        if(d == 0.0){
            countZero++;
        }
    }

    return (double)countZero / (double)vec.size();
}

double computePearsonCorrelation(const vector<double> &x, const vector<double> &y, double x_mean, double x_stddev, double y_mean, double y_stddev){
	double sum = 0.0;
	unsigned int N = x.size();
    for(unsigned int i=0; i != N; ++i){
        sum += ((x[i]-x_mean)*(y[i]-y_mean));
    }
    return sum/((double)N*(x_stddev*y_stddev));
}

double computePearsonCorrelation(const vector<double> &x, const vector<double> &y){
    vector<double> vec;
    double x_mean = computeMean(x);
    double x_stddev = computeStandardDeviation(x);
    double y_mean = computeMean(y);
    double y_stddev = computeStandardDeviation(y);
    for(unsigned int i=0; i != x.size(); ++i){
        vec.push_back((x[i]-x_mean)*(y[i]-y_mean));
    }
    return computeMean(vec)/(x_stddev*y_stddev);
}

vector<double> computePearsonCorrelation(const vector<vector<double>> &M){
	int N = M.size();
	cout << endl << "Pearson correlation to be computed for " << N << " vectors..." << endl;
	vector<double> correlationMatrix(N*N);
	vector<double> means(N);
	vector<double> standard_deviations(N);

	cout << "Computing all means and standard deviations... " << flush;
	for(int i=0; i<N; ++i){
		means[i] = computeMean(M[i]);
		standard_deviations[i] = computeStandardDeviation(M[i]);
	}
	cout << "Done."  << endl;

	cout << "Computing correlations for all pairs of vectors... " << flush;
	for(int i=0; i<N; ++i){
		correlationMatrix[i+N*i] = 1.0;
		for(int j=i+1; j<N; ++j){
			double cor = computePearsonCorrelation(M[i], M[j], means[i], standard_deviations[i], means[j], standard_deviations[j]);
			correlationMatrix[i+N*j] = cor;
			correlationMatrix[j+N*i] = cor;
		}
	}
	cout << "Done. " << endl << flush;

	return correlationMatrix;
}

//Assume x is sorted
void computeRankSorted(vector<double> &x){
	unsigned int current = 0;
	while(current != x.size()){
		double d = x[current];
		unsigned int i = current+1;
		int sum = current+1;
		int count = 1;
		while(i < x.size() && d == x[i]){
			++count;
			sum += i+1;
			++i;
		}

		for(unsigned int j=current; j<i; ++j){
			x[j] = (double)sum/(double)count;
		}

		current = i;
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

double computeSpearmanCorrelationRankedVectors(const vector<double> &rankedX, const vector<double> &rankedY){
	return computePearsonCorrelation(rankedX, rankedY);
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
	for(vector<double> &vec : M_copy){
		computeRank(vec);
	}
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
