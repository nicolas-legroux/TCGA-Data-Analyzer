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

double computeCorrelation(const vector<double> &x, const vector<double> &y){
    vector<double> vec;
    double x_mean = computeMean(x);
    double x_stddev = computeStandardDeviation(x);
    double y_mean = computeMean(y);
    double y_stddev = computeStandardDeviation(y);
    for(auto i=0; i != x.size(); ++i){
        vec.push_back((x[i]-x_mean)*(y[i]-y_mean));
    }
    return computeMean(vec)/(x_stddev*y_stddev);
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


