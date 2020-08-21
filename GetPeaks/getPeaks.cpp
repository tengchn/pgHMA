#include<cstdlib>
#include<cstring>
#include<vector>
#include<fstream>
#include<iostream>
#include<algorithm>
#include<math.h>

#define EMPTY -10.0
#define DEFAULT_K 0
#define DEFAULT_L 1
#define DEFAULT_T 0.6
#define DEFAULT_T_LIST "0.3,0.6,1.0"
#define VERSION "1.4"

using namespace std;

struct Option {
	int k;
	int l;
	double t;
	double m; // the time of lower marker
	double n; // the time of upper marker
	string t_list;
	vector<double> t_array;
};


void getStrings(string& s, vector<string>& fields) {
	string t;
	int i;
	t="";
	fields.clear();
	for (i=0; i<s.length(); i++) {
		if (s[i] == ' ')
			continue; // skip the space
		if (s[i] == ',') {
			fields.push_back(t);
			t="";
		} else {
			t.append(1,s[i]);
		}
	}
	fields.push_back(t);
}

void getNumbers(string& s, vector<double>& fields) {
	string t;
	int i;
	t="";
	fields.clear();
	for (i=0; i<s.length(); i++) {
		if (s[i] == ' ')
			continue; // skip the space
		if (s[i] == ',') {
			if (t=="")
				fields.push_back(EMPTY);
			else
				fields.push_back(atof(t.c_str()));
			t="";
		} else {
			t.append(1,s[i]);
		}
	}
	if (t=="")
		fields.push_back(EMPTY);
	else
		fields.push_back(atof(t.c_str()));
}

// take the average : data[i] = average between i-k th and i+k th items in the raw data
// k can be zero
void averaging(vector<double>& data, vector<double>& avgData, vector<int>& pos, int k) {
	int i,j;
	double sum;
	bool isValid;
	for (i=k; i<data.size()-k; i++) {
		sum = 0.0;
		isValid = true;
		for (j=i-k; j<=i+k && isValid; j++) {
			sum += data[j];
			if (data[j] == EMPTY)
				isValid = false;
		}
		if (isValid) {
			avgData.push_back(sum / (2.0 * k + 1));
			pos.push_back(i);
		}
	}
}

// data[i] is local peak if data[i-l] <= data[i-l+1] <= data[i-1+2] <= ... <= data[i] >= data[i+1] >= data[i+2] >= data[i+l]
// minimum value of l = 1
void potentialPeaks(vector<double>& avgData, vector<int>& avgDataPos, vector<int>& potentialPeakPos, int l) {
	int i,j;
	bool isPeak;
	potentialPeakPos.clear();
	for (i=l; i<avgData.size()-l; i++) {
		isPeak = true;
		for (j=i-l; j<i; j++) {
			if (avgData[j] > avgData[j+1])
				isPeak = false;
		}
		if (!isPeak)
			continue;
		for (j=i; j<i+l; j++) {
			if (avgData[j] < avgData[j+1])
				isPeak = false;
		}
		if (isPeak) {
			potentialPeakPos.push_back(avgDataPos[i]);
		}
	}
}

// if data[i] is a potential peak, the peak value = max{rawdata[i-k], rawdata[i-k+1], ... , rawdata[i+k]}, with the corresponding positions
void getPeaks(vector<double>& data, vector<int>& potentialPeakPos, int k, vector<int>& peakPos) {
	int i, j, n;
	double maxValue;
	int maxPos;
	n = data.size();
	peakPos.clear();
	for (i=0; i<potentialPeakPos.size(); i++) {
		if (potentialPeakPos[i]>=k && potentialPeakPos[i]<n-k) {
			j = potentialPeakPos[i] - k;
			maxValue = data[j];
			maxPos = j;
			for (j=potentialPeakPos[i]-k+1; j<=potentialPeakPos[i]+k; j++) {
				if (data[j] > maxValue) {
					maxValue = data[j];
					maxPos = j;
				}
			}
		}
		peakPos.push_back(maxPos);
	}
}

// if the distance between two peaks is < k, then only output the one with higher value
void removeTooClosePeaks(vector<double>& data, vector<int>& peakPos, int k) {
	vector<int> replacedByPos;
	int i,j,n;
	n = peakPos.size();
	replacedByPos.resize(n);
	for (i=0; i<n; i++) {
		replacedByPos[i] = -1;
	}
	for (i=0; i<n-1; i++) {
		if (peakPos[i+1]-peakPos[i] <= k) {
			// too close, select the higher one
			if (data[peakPos[i]] > data[peakPos[i+1]]) {
				replacedByPos[i+1] = i; // replace the i+1 by i
			} else {
				replacedByPos[i] = i+1; // replace the i by i+1
			}
		}
	}
	// remove the peakPos[i] if replacedByPos[i] != -1
	j=0;
	for (i=0; i<n; i++) {
		if (replacedByPos[i] == -1) {
			if (i>j) {
				peakPos[j] = peakPos[i];
				j++;
			}
		}
	}
	peakPos.resize(j);
}

// get the average
double getAvg(vector<double>& data) {
	int n = data.size();
	double sum = 0.0;
	int i;
	for (i=0; i<n; i++) {
		sum += data[i];
	}
	return sum/n;
}

// get the standard deviation
double getStdDev(vector<double>& data, double& avg) {
	int n = data.size();
	double stddev = 0.0;
	int i;
	double t;
	for (i=0; i<n; i++) {
		t = data[i] - avg;
		stddev += t*t;
	}
	stddev = sqrt(stddev/(double)(n-1));
	return stddev;
}

// background_noise = the value v such that b of the data are lower than v and (1-b) of the data are above v" << endl;
double backgrdValue(vector<double>& data, double r) {
	int n = data.size() * r;
	if (n >= data.size())
		n = data.size()-1;
	vector<double> sortedData(data);
	sort(sortedData.begin(), sortedData.end());
	return sortedData[n];
}

// show the peaks with value higher than the background noise
// the threshold for the background noise = avgValue + t * stdDev
void showPeak(string header, vector<double>& timeData, vector<double>& data, vector<int>& peakPos, double avgValue, double stdDev, double t) {
	int j,p;
	double thres;
	cout << "--------------" << endl;
	cout << header << endl;
	cout << "--------------" << endl;
	// cout << "Background noise: " << avgValue << endl;
	cout << "Average: " << avgValue << endl;
	cout << "Standard deviation: " << stdDev << endl;
	thres = avgValue + t * stdDev;
	cout << "Noise threshold: " << thres << endl;
	cout << "time\tvalue" << endl;
	for (j=0; j<peakPos.size(); j++) {
		p = peakPos[j];
		if (data[p] >= thres) {
			cout << timeData[p] << "\t" << data[p] << endl;
		}
	}
	cout << endl;
}

// identify the peaks
bool identifyPeak(string header, vector<double>& timeData, vector<double>& data, vector<int>& peakPos, double avgValue, double stdDev, double t, double m, double n, bool proceedEvenError) {

	vector<int> validPeaks;
	vector<int> peakTypes; // 0: unidentify; 1: homo; 2: hetero
	int lm, hm;
	double lm_dist, hm_dist;
	double d, d2;
	int j,p;
	double thres;
	double dHe, dHo;
	double dHm, dLm;
	string errMesg = "";
	
	// result =  (dHe-dHo)/(dUM - dLM )
	double result;
	int validResult;

	// identify which lower marker and higher marker	
	lm=hm=-1;
	thres = avgValue + t * stdDev;
	for (j=0; j<peakPos.size(); j++) {
		p = peakPos[j];
		if (data[p] >= thres) {
			d = fabs(timeData[p] - m);
			if (lm==-1 || d < lm_dist) {
				lm = j;
				lm_dist = d;
			}
			d = fabs(timeData[p] - n);
			if (hm==-1 || d < hm_dist) {
				hm = j;
				hm_dist = d;
			}
		}
	}
	// collect all the valid peaks excluding lower marker, higher marker and the peaks before lower marker
	for (j=lm+1; j<peakPos.size(); j++) {
		if (j==hm)
			continue;
		p = peakPos[j];
		if (data[p] >= thres) {
			validPeaks.push_back(p);
			peakTypes.push_back(0);
		}
	}

	validResult = 1;	
	if (validPeaks.size() == 1) {
		peakTypes[0] = 1;
		dHe = dHo = timeData[validPeaks[0]];
	} else if (validPeaks.size() == 2) {
		peakTypes[0] = 1;
		peakTypes[1] = 2;
		dHo = timeData[validPeaks[0]];
		dHe = timeData[validPeaks[1]];
	} else if (validPeaks.size() == 3) {
		d  = timeData[validPeaks[1]]-timeData[validPeaks[0]];
		d2 = timeData[validPeaks[2]]-timeData[validPeaks[1]];
		if (d <= d2) {
			peakTypes[0] = peakTypes[1] = 1;
			peakTypes[2] = 2;
			dHo = (timeData[validPeaks[0]] + timeData[validPeaks[1]]) / 2.0;
			dHe = timeData[validPeaks[2]];
		} else {
			peakTypes[0] = 1;
			peakTypes[1] = peakTypes[2] = 2;
			dHo = timeData[validPeaks[0]];
			dHe = (timeData[validPeaks[1]] + timeData[validPeaks[2]]) / 2.0;
		}
	} else if (validPeaks.size() == 4) {
		peakTypes[0] = peakTypes[1] = 1;
		peakTypes[2] = peakTypes[3] = 2;
		dHo = (timeData[validPeaks[0]] + timeData[validPeaks[1]]) / 2.0;
		dHe = (timeData[validPeaks[2]] + timeData[validPeaks[3]]) / 2.0;
	} else if (validPeaks.size() == 0) {
		validResult = 0;
		if (!proceedEvenError)
			return false;
		errMesg = "Error! No peak is left after excluding the lower marker and the higher marker";
	} else {
		validResult = 0;
		if (!proceedEvenError)
			return false;
		errMesg = "Error! There are more than 4 peaks after excluding lower marker and higher marker";
	}


	cout << "--------------" << endl;
	cout << header << endl;
	cout << "t = " << t << endl;
	cout << "--------------" << endl;

	// show error message if there is
	if (errMesg != "") {
		cout << errMesg << endl;
	}

	// show the result
	if (validResult) {
		dLm = timeData[peakPos[lm]];
		dHm = timeData[peakPos[hm]];
		result = (dHe - dHo) / (dHm - dLm);
		cout << "The value of M = " << result << endl;
	}

	cout << "time\tvalue\ttype" << endl;
	p = peakPos[lm];
	cout << timeData[p] << "\t" << data[p] << "\tlower marker" << endl;
	p = peakPos[hm];
	cout << timeData[p] << "\t" << data[p] << "\thigher marker" << endl;
	
	for (j=0; j<validPeaks.size(); j++) {
		p = validPeaks[j];
		cout << timeData[p] << "\t" << data[p] << "\t";
		if (peakTypes[j]==1)
			cout << "homoduplex" << endl;
		else if (peakTypes[j]==2)
			cout << "heteroduplex" << endl;
		else
			cout << "unknown" << endl;
	}
	cout << endl;
	return (errMesg == "");
}

void showComputeMenu(char* argv0) {
		// cout << "GetPeaks v" << VERSION << endl;
		// cout << endl;
		cout << "To identify the homo-duplex and hetero-duplex, assuming it is a mixture of TWO samples" << endl;
		cout << "Syntax: " << argv0 << " compute [csv file] -m [time (in double) of lower marker] -n [time (in double) of upper marker] (options)" << endl;
		cout << "Options:" << endl;
		cout << "    -k [integer value of k, default: " << DEFAULT_K << "]" << endl;
		cout << "    -l [integer value of l, default: " << DEFAULT_L << "]" << endl;
		cout << "    -tlist [list of t values, default: " << DEFAULT_T_LIST << "]" << endl;
		cout << endl;
		cout << "Procedure of homo-duplex and hetero-duplex identification:" << endl;
		cout << "1. use the above procedure (with t=t1) to identify the peaks" << endl;
		cout << "2. the peak close to the input lower marker value is regarded as lower marker" << endl;
		cout << "3. the peak close to the input higher marker value is regarded as higher marker" << endl;
		cout << "4. ignore the lower-marker peak, the higher-marker peak, and all the peaks before the lower-marker peak." << endl;
		cout << "5. if the number of remaining peaks is more than 4, then repeat step 1 to step 4 with t=t2 and" << endl; 
		cout << "   report error if the number of remaining peaks is still more than 4. Otherwise, proceed to step 6." << endl;
		cout << "6. there are 5 cases:" << endl;
		cout << "   a. only one peak (p1) is left, then p1 is homo-duplex" << endl;
		cout << "   b. only two peaks (p1, p2) are left, then p1 is homo-duplex and p2 is hetero-duplex" << endl;
		cout << "   c. three peaks (p1, p2, p3) are left and p2 is closer to p1, then p1 & p2 are homo-duplex, and p3 is hetero-duplex" << endl;
		cout << "   d. three peaks (p1, p2, p3) are left and p2 is closer to p3, then p1 is homo-duplex, and p2 & p3 are hetero-duplex" << endl;
		cout << "   e. four peaks (p1, p2, p3, p4) are left. Then p1 & p2 are homo-duplex, and p3 & p4 are hetero-duplex" << endl;
		cout << "For example:" << endl;
		cout << "   $ " << argv0 << " compute example.csv -m 11.92 -n 22.65" << endl;
}

void showHelpMenu(char* argv0) {
		cout << "GetPeaks v" << VERSION << endl;
		cout << endl;
		cout << "To get the peaks from the input data" << endl;
		cout << "Syntax: " << argv0 << " view [csv file] (options)" << endl;
		cout << endl;
		cout << "Options:" << endl;
		cout << "    -k [integer value of k, default: " << DEFAULT_K << "]" << endl;
		cout << "    -l [integer value of l, default: " << DEFAULT_L << "]" << endl;
		cout << "    -t [double value of t, default: " << DEFAULT_T << "]" << endl;
		cout << endl;
		cout << "Procedure of peak identification:" << endl;
		cout << "1. take the average : data[i] = average between i-k th and i+k th items in the raw data" << endl;
		cout << "2. data[i] is a potential peak if data[i-l] <= ... <= data[i] >= data[i+1] >= ... >= data[i+l]" << endl;
		cout << "3. if data[i] is a potential peak, the peak value = max{rawdata[i-k], ... , rawdata[i+k]}" << endl;
		cout << "4. if the distance between two peaks <= k, then only output the one with higher value" << endl;
		cout << "5. the threshold (thres) for background noise = average of all values + t * their standard deviation" << endl;
		cout << "6. show the peaks with value >= thres" << endl;
		cout << "For example:" << endl;
		cout << "   $ " << argv0 << " view example.csv" << endl;
		cout << endl;
		showComputeMenu(argv0);
}

void process_T_array(string t_array, vector<double> & t_arr) {
	string s;
	int i;
	int p;
	
	// get t_array
	t_arr.clear();
	i=0;
	while (i<t_array.length()) {
		p=i;
		while (p<t_array.length() && t_array[p]!=',')
			p++;
		if (p > i) {
			s = t_array.substr(i, p-i);
			t_arr.push_back(atof(s.c_str()));
		}
		i=p+1;
	}
	
	// check whether the size is at least 1
	if (t_arr.size() < 1) {
		cerr << "Error! There is no item inside -t!" << endl;
		exit(1);
	}
}

void getOptions(int argc, char** argv, Option& option) {
	int i;
	option.k = DEFAULT_K;
	option.l = DEFAULT_L;
	option.t = DEFAULT_T;
	option.t_list = DEFAULT_T_LIST;
	option.m = option.n = -1;
	for (i=3; i<argc; i+=2) {
		if (argc > i+1) {
			if (strcmp(argv[i],"-k")==0) {
				option.k = atoi(argv[i+1]);
			} else if (strcmp(argv[i],"-l")==0) {
				option.l = atoi(argv[i+1]);
			} else if (strcmp(argv[i],"-t")==0) {
				option.t = atof(argv[i+1]);
			} else if (strcmp(argv[i],"-m")==0) {
				option.m = atof(argv[i+1]);
			} else if (strcmp(argv[i],"-n")==0) {
				option.n = atof(argv[i+1]);
			} else if (strcmp(argv[i],"-t_list")==0) {
				option.t_list = argv[i+1];
			}
		}
	}
	process_T_array(option.t_list, option.t_array);
}

int main(int argc, char** argv) {
	vector<vector<double> > data;
	vector<string> header;
	vector<double> row;
	vector<double> avgData;
	vector<int> avgDataPos;
	vector<int> potentialPeakPos;
	vector<int> peakPos;
	ifstream fin;
	string aline;
	string action;
	double avg;
	double stddev;
	int i,j,n;
	Option option;
	bool returnStatus;
	bool proceedEvenError;
	int numErrors = 0;
	
	if (argc < 3) {
		showHelpMenu(argv[0]);
		exit(1);
	}
	
	action = argv[1];
	
	if (action!="view" && action!="compute") {
		cerr << "Error! Invalid keyword: " << action << endl << endl;
		showHelpMenu(argv[0]);
		exit(1);
	}
	
	fin.open(argv[2]);
	
	getOptions(argc, argv, option);
	
	if (action=="compute" && (option.m == -1 || option.n == -1)) {
		cerr << "Error! The value of m or n is not inputted" << endl << endl;
		showComputeMenu(argv[0]);
		exit(1);
	}

	if (action=="compute" && (option.m >= option.n)) {
		cerr << "Error! The value of m should smaller than n" << endl << endl;
		showComputeMenu(argv[0]);
		exit(1);
	}

	cout << "GetPeak v" << VERSION << endl;
	cout << "Value of l = " << option.l << endl;
	
	if (action=="compute") {
		cout << "The list of t values: " << option.t_list << endl;
		cout << "Time of lower marker = " << option.m << endl;
		cout << "Time of higher marker = " << option.n << endl;
	} else {
		cout << "Value of k = " << option.k << endl;
		cout << "Value of t = " << option.t << endl;
	}

	// get the header
	if (getline(fin,aline)) {
		getStrings(aline, header);
	} else {
		cerr << "Error! The file is empty!" << endl;
		exit(1);
	}

	// reset the data
	n = header.size();

	data.resize(n);
	for (i=0; i<n; i++) {
		data[i].clear();
	}
		
	// get the content
	j = 2;
	while (getline(fin,aline)) {
		if (aline.length()>0) {
			getNumbers(aline, row);
		}
		if (row.size() > n) {
			cerr << "Error! Line " << j << " has more than " << n << " fields" << endl;
			exit(1);
		}
		for (i=0; i<n; i++) {
			if (i < row.size()) {
				data[i].push_back(row[i]);
			} else {
				data[i].push_back(EMPTY);
			}
		}
		j++;
	}
	fin.close();

	for (i=1; i<n; i++) {
		avgData.clear();
		// show the data[i]
		// for (j=0; j<data[i].size(); j++)
		//	cout << data[i].at(j) << " ";
		// cout << endl;
		// take the average : data[i] = average between i-k th and i+k th items in the raw data
		// k can be zero
		averaging(data[i], avgData, avgDataPos, option.k);
		// data[i] is local peak if data[i-l] <= data[i-l+1] <= data[i-1+2] <= ... <= data[i] >= data[i+1] >= data[i+2] >= data[i+l]
		potentialPeaks(avgData, avgDataPos, potentialPeakPos, option.l);
		// if data[i] is a potential peak, the peak value = max{rawdata[i-k], rawdata[i-k+1], ... , rawdata[i+k]}, with the corresponding positions
		getPeaks(data[i], potentialPeakPos, option.k, peakPos);
		// if the distance between two peaks is < k, then only output the one with higher value
		removeTooClosePeaks(data[i], peakPos, option.k);
	
		// computation of background noise
		avg = getAvg(data[i]);
		stddev = getStdDev(data[i], avg);
		
		if (action=="view") {
			// show the peaks with value higher than the background noise
			// the threshold for the background noise = avgValue + t * stdDev
			showPeak(header[i], data[0], data[i], peakPos, avg, stddev, option.t);
		} else {
			// identify the peaks
			returnStatus = false;
			j = 0;
			while (!returnStatus && j<option.t_array.size()) {
				// cout << "THE VALUE OF T IS: " << option.t_array[j] << endl;
				proceedEvenError = (j==option.t_array.size()-1); // true for the last iteration
				returnStatus = identifyPeak(header[i], data[0], data[i], peakPos, avg, stddev, option.t_array[j], option.m, option.n, proceedEvenError);
				j++;
			}
			if (!returnStatus)
				numErrors++;
		}
	}
	
	// show the summary
	cout << "---------------------------------------------" << endl;
	cout << "Total number of samples: " << n-1 << endl;
	if (action != "view") {
		cout << endl;
		cout << "Number of samples where the heteroduplex was" << endl;
		cout << "successfully identified: " << n - numErrors - 1 << endl;
		cout << endl;
		cout << "Number of samples where the heteroduplex" << endl;
		cout << "could not be identified: " << numErrors << endl;
	}
	cout << "---------------------------------------------" << endl;
}
