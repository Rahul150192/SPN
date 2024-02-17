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

#define BS 128
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
 * GIFT-128 implementation adapted from
 * https://github.com/giftcipher/gift
 */
const unsigned char GIFT_S[16] = {1, 10,4,12, 6,15, 3, 9, 2,13,11, 7, 5, 0, 8,14};
const unsigned char GIFT_P128[]= {
        0, 33, 66, 99, 96,  1, 34, 67, 64, 97,  2, 35, 32, 65, 98,  3,
        4, 37, 70,103,100,  5, 38, 71, 68,101,  6, 39, 36, 69,102,  7,
        8, 41, 74,107,104,  9, 42, 75, 72,105, 10, 43, 40, 73,106, 11,
        12, 45, 78,111,108, 13, 46, 79, 76,109, 14, 47, 44, 77,110, 15,
        16, 49, 82,115,112, 17, 50, 83, 80,113, 18, 51, 48, 81,114, 19,
        20, 53, 86,119,116, 21, 54, 87, 84,117, 22, 55, 52, 85,118, 23,
        24, 57, 90,123,120, 25, 58, 91, 88,121, 26, 59, 56, 89,122, 27,
        28, 61, 94,127,124, 29, 62, 95, 92,125, 30, 63, 60, 93,126, 31
};

void gift_round_128(unsigned char* input, unsigned char *K){
    unsigned char bits[128], perm_bits[128];
    int i, j ;

    // S-box operation
    for (i=0; i<32; i++){
        input[i] = GIFT_S[input[i]];
    }

    //PermBits
    //input to bits

    for (i=0; i<32; i++){
        for (j=0; j<4; j++){
            bits[4*i+j] = (input[i] >> j) & 0x1;
        }
    }
    //permute the bits
    for (i=0; i<128; i++){
        perm_bits[GIFT_P128[i]] = bits[i];
    }
    //perm_bits to input
    for (i=0; i<32; i++){
        input[i]=0;
        for (j=0; j<4; j++){
            input[i] ^= perm_bits[4*i+j] << j;
        }
    }

    // Key XOR
    for(i = 0 ; i<32; i++){
        input[i] ^= K[i] ;
    }
}

void convert_to_bits(int *X, unsigned char*S, unsigned char*Z){
    unsigned char temp[32] ;
    for(int i = 0; i< 32; i++){
        temp[i] = (S[i] ^ Z[i]) ;
    }
    for(int i = 0; i < 32; i++){
        for(int j = 0; j<4; j++){
            X[4*i +j] = ( temp[i] >> (j) ) & 1 ;
        }
    }
}

void dimension(unsigned char *A, int rounds, int steps, int *temp, int option, int *temp_uniq_diff){
    srand(time(NULL) + clock());
    int n = 128 ;
    int N2 = n + (n*(n-1))/2 ;
    int N3 = n + (n*(n-1))/2 + ((n-1)*(n-2))/2 ;

    Mat<GF2> M1, M2, M3;
    M1.SetDims(steps, n);
    M2.SetDims(steps, N2);
    M3.SetDims(steps, N3);

    int i, j, k, r, l;

    int X[128];
    unsigned char S[32], Z[32], K[32];

    for(r = 1 ; r<rounds; r++){
        map< bitset<BS>, int, cmpBitset> countingBox;
        for(i = 0; i<steps; i++){
            for(j = 0 ; j< 32; j++){
                S[j] = rand() % 16;
                Z[j] = S[j] ^ A[j];
            }
            for(j = 0 ; j<r; j++){
                for(k = 0 ; k< 32; k++){
                    K[k] = rand() % 16 ;
                }
                gift_round_128(S, K);
                gift_round_128(Z, K);
            }
            convert_to_bits(X, S, Z) ;

            // Add difference value to hash map
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
    int n = 128;
    unsigned char A[32] = {0};

    for(i = 0; i<32; i++){
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
        N = 128;
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

    printf("#### GIFT-128 deterministic distinguishers ####\n");
    printf("#### First-order ####\n");
    avg_rank(0, 0xC & 0xF, 100, 5, 148, 1);
    printf("---------------------------------------------------------------------------\n");

    printf("#### Second-order ####\n");
    avg_rank(0, 0x2 & 0xF, 100, 5, 8784, 2);
    printf("---------------------------------------------------------------------------\n");

    printf("#### Third-order ####\n");
    avg_rank(0, 0xC & 0xF, 100, 5, 17286, 3);
    printf("---------------------------------------------------------------------------\n");




    printf("#### GIFT-128 probabilistic distinguishers ####\n");
    printf("#### First-order ####\n");
    avg_rank(0, 0x6 & 0xF, 100, 5, 148, 1);
    printf("---------------------------------------------------------------------------\n");

    printf("#### Second-order ####\n");
    avg_rank(0, 0x6 & 0xF, 100, 5, 8784, 2);
    printf("---------------------------------------------------------------------------\n");

    printf("#### Third-order ####\n");
    avg_rank(0, 0xE & 0xF, 100, 7, 17286, 3);
    printf("---------------------------------------------------------------------------\n");



}





