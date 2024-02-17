#include <iostream>
#include <random>
#include <vector>
#include <NTL/mat_GF2.h>

using namespace std;
using namespace NTL;


float dimension(int &option, int &m, int &n){
	    int s = 1000;
	    int steps = s;
	    float totalRank=0;  
	    std::random_device rd;
	    std::mt19937 gen(rd());
	    std::uniform_int_distribution<> dis(0, 1);
	    
	if(option == 1){
	    
            while(s != 0){
			int* v = new int[n];
			
			std::vector<std::vector<int>> M;
			
			
			for (int row = 0 ; row < m ; row++){
			    
			    for (int i = 0; i < n; ++i) {
				v[i] = dis(gen);
				//cout << v[i] << " ";
			    }
			    
			    vector<int> row_vector;
			    
			    for (int i = 0; i < n; ++i) {
				row_vector.push_back(v[i]);
			    }
			
			
			    
			    M.push_back(row_vector);
			}
			
			
		   
		     
		      int vector_size = n;
		    
		    
		  //  cout<<"Full vector size : "<<vector_size<<"\n";
		    
		    Mat<GF2> M1;
			M1.SetDims(m, vector_size );
			for (int i=0; i < m ; i++){
				for(int j=0; j<vector_size; j++){
					M1[i][j] = M[i][j];
				
				}
			
			}
			
		     
			
		    int rank = gauss(M1);
		    totalRank += rank;
		 //   cout<<rank<<" "<<totalRank<<"\n";
		    s--;
		    } // end while 
		    
	}
	else if(option == 2){
		while(s != 0){
			int* v = new int[n];
			
			std::vector<std::vector<int>> M;
			
			
			for (int row = 0 ; row < m ; row++){
			    
			    for (int i = 0; i < n; ++i) {
				v[i] = dis(gen);
				//cout << v[i] << " ";
			    }
			    
			    vector<int> row_vector;
			    
			    for (int i = 0; i < n; ++i) {
				row_vector.push_back(v[i]);
			    }
			
			    for (int i = 0; i < n; ++i) {
				for (int j = i + 1; j < n; ++j) {
				    row_vector.push_back(v[i] * v[j]);
				}
			    }
			
			
		       
			    
			    M.push_back(row_vector);
			}
			
			
		    
		    
		    int vector_size = n + n*(n-1)/2;
		    
		    
		    
		  
		    
		    Mat<GF2> M1;
			M1.SetDims(m, vector_size );
			for (int i=0; i < m ; i++){
				for(int j=0; j<vector_size; j++){
					M1[i][j] = M[i][j];
				
				}
			
			}
			
		      
			
		    int rank = gauss(M1);
		    totalRank += rank;
		 //   cout<<rank<<" "<<totalRank<<"\n";
		    s--;
		    } // end while
	}
	else if(option == 3){
		while(s != 0){
			int* v = new int[n];
			
			std::vector<std::vector<int>> M;
			
			
			for (int row = 0 ; row < m ; row++){
			    
			    for (int i = 0; i < n; ++i) {
				v[i] = dis(gen);
				//cout << v[i] << " ";
			    }
			    
			    vector<int> row_vector;
			    
			    for (int i = 0; i < n; ++i) {
				row_vector.push_back(v[i]);
			    }
			
			    for (int i = 0; i < n; ++i) {
				for (int j = i + 1; j < n; ++j) {
				    row_vector.push_back(v[i] * v[j]);
				}
			    }
			
			
			    for (int j = 1; j < n; ++j) {
				for (int k = j + 1; k < n; ++k) {
				    row_vector.push_back(v[0] * v[j] * v[k]);
				}
			    }
			    
			    M.push_back(row_vector);
			}
			
			
		  
		    
		    int vector_size = n + n*(n-1)/2 + (n-1)*(n-2)/2;
		    
		   
		    
		    
		    
		  //  cout<<"Full vector size : "<<vector_size<<"\n";
		    
		    Mat<GF2> M1;
			M1.SetDims(m, vector_size );
			for (int i=0; i < m ; i++){
				for(int j=0; j<vector_size; j++){
					M1[i][j] = M[i][j];
				
				}
			
			}
			
		      
			
		    int rank = gauss(M1);
		    totalRank += rank;
		 //   cout<<rank<<" "<<totalRank<<"\n";
		    s--;
		    } // end while
	}
	else{
		printf("Invalid option!. Enter valid option 1,2 or 3");
	}
	
	float avgRank = totalRank/steps ;
	return avgRank;
}


int main(){
	int option = 0;
	printf("Enter option-1 for first order, 2 for second order, 3 for third order: ");
	int item1 = scanf("%d",&option);
	int m,n;
	printf("Enter m for number of samples: ");
	int item2 = scanf("%d",&m);
	printf("Enter n for block size: ");
	int item3 = scanf("%d",&n);
	float avgRank = dimension(option, m, n);
	cout<<"Average rank is " << avgRank<<"\n";
}


