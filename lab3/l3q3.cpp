#include <bits/stdc++.h>
using namespace std;

void u_and_d(double &u, double &d, double T, double M, double sig, double r){
   u = exp(sig*sqrt(T/M)+(r-0.5*sig*sig)*(T/M));
   d = exp(-sig*sqrt(T/M)+(r-0.5*sig*sig)*(T/M));
}

double approx(double x){
	int p=(int)(10000.0*x);
	double ans=p/10000.0;
	return ans;
	//return x;
}


double option_price(double s0, double T, double r, double sig, double M ){

    double u, d;
    u_and_d(u, d, T, M, sig, r);
    double p = (exp(r*(T/M))- d)/(u-d);

    map< pair<double,double> , unordered_set<double> >mp;

    mp[{0,0}].insert(s0);
    int ans  = 0;
    double sp;
    for(int i = 1; i<= (M); i++){

        for(int j = 0; j <= i; j++){

            sp = approx(s0*(pow(u, i-j)*pow(d, j)));
            int flag = 1;

            if(j > 0){
                for( auto k: mp[{i-1, j-1 }]){
		            k = approx(k);
                    if(sp <= k){
                        mp[{i, j}].insert(k);
                    }
                    else if(flag){
                        flag = 0;
                        mp[{i, j}].insert(sp);
                    }

                }

            }
            if(j < i){
                for( auto k: mp[{i-1, j }]){
		            k = approx(k);
                    if(sp <= k){
                        mp[{i, j}].insert(k);
                    }
                    else if(flag){
                        flag = 0;
                        mp[{i, j}].insert(sp);
                    }

                }

            }

        }
    }



    map< pair<double , double>, double>final_price;
    map< pair<double , double>, double>curr_price;
    for(int i = 0; i<= M; i++){

        for(auto j: mp[{M,i}]){

        	   double a = approx(s0*pow(u, M-i)*pow(d, i));
            double b = approx(j);
            final_price[{a, b}] = -a + b;

        }
    }

    for(int i = M-1; i>= 0; i--){
        for(int j = 0; j<=i; j++ ){
            for(auto k : mp[{i,j}]){

	            double a = approx(s0*pow(u, i-j)*pow(d, j));
               double b = approx(k);
               double c = approx(s0*pow(u, i-j)*pow(d, j)*u);
               double d1 = approx(s0*pow(u, i-j)*pow(d, j)*d);
	            curr_price[{a, b}] = (p*final_price[{c, max(c, b)}] + (1-p)*final_price[{d1, b}])*exp(-r*(T/M));
            }
        }
        final_price = curr_price;
    }

    return final_price[{s0, s0}];
}



int main(){


    double s0=100, T=1, r=0.08, sig=0.2, M;
    for(M=5; M<=25; M+=5){
        cout << "M = " << M << " : " << option_price(s0,T,r,sig,M) << "\n";
    }
    M=50;
    cout << "M = " << M << " : " << option_price(s0,T,r,sig,M) << "\n";
    return 0;
}
