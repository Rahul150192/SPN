#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <NTL/mat_GF2.h>
#include<vector>
#include<bitset>
#include<algorithm>
#include<map>
#include<iomanip>
#include<cmath>

using namespace std;
using namespace NTL;

#define BS 64
struct cmpBitset{
    bool operator()(const bitset<BS>& a, const bitset<BS>& b) const {
        for (int i = 0; i < BS; i++) {
            if (a[i] < b[i])
                return true;
            else if (a[i] > b[i])
                return false;
        }
        return false;
    }
};


int count_unique_diff(map< bitset<BS>, int, cmpBitset >& countingBox){
    int count = 0 ;
    auto it2 = countingBox.begin();
    while (it2 != countingBox.end()) {
        count++;
        it2++;
    }
    return count ;

}

/*
 * GIFT-64 implementation adapted from
 * https://github.com/giftcipher/gift
 */
const unsigned char GIFT_S[16] = {1, 10,4,12, 6,15, 3, 9, 2,13,11, 7, 5, 0, 8,14};
const unsigned char GIFT_P64[]={
        0, 17, 34, 51, 48,  1, 18, 35, 32, 49,  2, 19, 16, 33, 50,  3,
        4, 21, 38, 55, 52,  5, 22, 39, 36, 53,  6, 23, 20, 37, 54,  7,
        8, 25, 42, 59, 56,  9, 26, 43, 40, 57, 10, 27, 24, 41, 58, 11,
        12, 29, 46, 63, 60, 13, 30, 47, 44, 61, 14, 31, 28, 45, 62, 15
};

void gift_round_64(unsigned char* input, unsigned char *K){

    unsigned char bits[64], perm_bits[64];
    int i, j ;

    // S-box operation
    for (i=0; i<16; i++){
        input[i] = GIFT_S[input[i]];
    }

    //PermBits
    //input to bits

    for (i=0; i<16; i++){
        for (j=0; j<4; j++){
            bits[4*i+j] = (input[i] >> j) & 0x1;
        }
    }
    //permute the bits
    for (i=0; i<64; i++){
        perm_bits[GIFT_P64[i]] = bits[i];
    }

    //perm_bits to input
    for (i=0; i<16; i++){
        input[i]=0;
        for (j=0; j<4; j++){
            input[i] ^= perm_bits[4*i+j] << j;
        }
    }
    // Key XOR
    for(i = 0 ; i<16; i++){
        input[i] ^= K[i] ;
    }
}


void convert_to_bits(int *X, unsigned char*S, unsigned char*Z){
    unsigned char temp[16] ;

    for(int i = 0; i< 16; i++){
        temp[i] = (S[i] ^ Z[i]) ;
    }

    for(int i = 0; i < 16; i++){
        for(int j = 0; j<4; j++){
            X[4*i +j] = ( temp[i] >> (j) ) & 1 ;
        }
    }
}

void dimension(unsigned char *A, int rounds, int steps, int *temp, int option, int *temp_uniq_diff){

    srand(time(NULL) + clock());
    int n = 64 ;
    int N2 = n + (n*(n-1))/2 ;
    int N3 = n + (n*(n-1))/2 + ((n-1)*(n-2))/2 ;

    Mat<GF2> M1, M2, M3;
    M1.SetDims(steps, n);
    M2.SetDims(steps, N2);
    M3.SetDims(steps, N3);

    int i, j, k, r, l;

    int X[64];
    unsigned char S[16], Z[16], K[16];

    for(r = 1 ; r<rounds; r++){
        map< bitset<BS>, int, cmpBitset> countingBox;
        for(i = 0; i<steps; i++){
            for(j = 0 ; j< 16; j++){
                S[j] = rand() % 16;
                Z[j] = S[j] ^ A[j];
            }
            for(j = 0 ; j<r; j++){
                for(k = 0 ; k< 16; k++){
                    K[k] = rand() % 16 ;
                }
                gift_round_64(S, K);
                gift_round_64(Z, K);
            }
            convert_to_bits(X, S, Z) ;

            //Add difference value to hash map
            bitset<BS> tmp;
            for (j = 0; j < BS; j++) {
                if (X[j] == 1) tmp[j] = 1;
                else tmp[j] = 0;
            }
            countingBox[tmp]++;

            if(option == 1){
                for(j = 0 ; j< n; j++){
                    M1[i][j] = X[j];
                }
            }
            else if(option == 2){
                for(j = 0 ; j< n; j++){
                    M2[i][j] = X[j];
                }
                l = n ;
                for(j = 0 ; j < n; j++){
                    for(k = j+1; k<n; k++){
                        M2[i][l] = X[j] & X[k] ;
                        l = l + 1 ;
                    }
                }
            }
            else if(option == 3){
                for(j = 0 ; j< n; j++){
                    M3[i][j] = X[j];
                }
                l = n ;
                for(j = 0 ; j < n; j++){
                    for(k = j+1; k<n; k++){
                        M3[i][l] = X[j] & X[k] ;
                        l = l + 1 ;
                    }
                }
                l = N2 ;
                for(j = 0 ; j<1; j++){
                    for(k = j+1; k<n; k++){
                        for(int q = k+1; q<n; q++){
                            M3[i][l] = X[j] & X[k] & X[q] ;
                            l = l + 1 ;
                        }
                    }
                }
            }
        }
        if(option == 1){
            temp[r] = gauss(M1);
        }
        else if(option == 2){
            temp[r] = gauss(M2);
        }
        else if(option == 3){
            temp[r] = gauss(M3);
        }
        temp_uniq_diff[r] = count_unique_diff(countingBox);
    }
}


void avg_rank(int s, unsigned char subspace, int exp, int r, int samples, int option){

    int dim[exp][r], uniq_diff[exp][r], temp[r], sum[r], sum_temp[r], temp_uniq_diff[r];
    int i, j, k;
    int n = 64;
    unsigned char A[16] = {0};

    for(i = 0; i<16; i++){
        A[i] = 0 ;
    }

    for(i = 0 ; i<r;i++){
        temp[i] = 0 ;
        sum[i] = 0 ;
        sum_temp[i] = 0 ;
        temp_uniq_diff[i] = 0 ;
    }

    for(i = 0 ; i<exp; i++){
        for(j = 0 ; j<r; j++){
            dim[i][j] = 0 ;
            uniq_diff[i][j] = 0 ;
        }
    }

    A[s] = subspace & 0xF ;

    for(i = 0; i<exp; i++){
        dimension(A, r, samples, temp, option, temp_uniq_diff);
        for(k = 0 ; k<r; k++){
            dim[i][k] = temp[k];
            uniq_diff[i][k] = temp_uniq_diff[k];
        }
    }

    for(i = 0 ; i<exp; i++){
        for(k = 0 ; k<r; k++){
            sum[k] += dim[i][k] ;
            sum_temp[k] += uniq_diff[i][k] ;
        }
    }

    int N;
    if(option == 1){
        N = 64;
    }
    else if(option == 2){
        N = n + (n*(n-1))/2 ;
    }
    else if(option == 3){
        N = n + (n*(n-1))/2 + ((n-1)*(n-2))/2 ;
    }

    printf("Round \t Avg. Dimension \t Final Dimension \t Avg. Unique Differences \n");
    for(i = 0 ; i<r; i++){
        printf("%d \t %9.2f \t\t %d \t\t\t%9.2f\n", i, sum[i]/(1.0*exp), N, sum_temp[i]/(1.0*exp));
    }

}


int main(){
    printf("#### GIFT-64 ####\n");
    printf("#### First-order ####\n");
    avg_rank(0, 0xE & 0xF, 100, 4, 84, 1);
    printf("---------------------------------------------------------------------------\n");

    printf("#### Second-order ####\n");
    avg_rank(0, 0xC & 0xF, 100, 4, 2226, 2);
    printf("---------------------------------------------------------------------------\n");

    printf("#### Third-order ####\n");
    avg_rank(6, 0xE & 0xF, 100, 6, 4302, 3);
    printf("---------------------------------------------------------------------------\n");
}



