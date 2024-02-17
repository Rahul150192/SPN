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
 * SKINNY-128 code adapted from:
 * https://sites.google.com/site/skinnycipher/downloads
 */

const unsigned char P[16] = {0,1,2,3,7,4,5,6,10,11,8,9,13,14,15,12};
const unsigned char sbox_8[256] ={
        0x65 , 0x4c , 0x6a , 0x42 , 0x4b , 0x63 , 0x43 , 0x6b , 0x55 , 0x75 , 0x5a , 0x7a ,
        0x53 , 0x73 , 0x5b , 0x7b ,0x35 , 0x8c , 0x3a , 0x81 , 0x89 , 0x33 , 0x80 , 0x3b ,
        0x95 , 0x25 , 0x98 , 0x2a , 0x90 , 0x23 , 0x99 , 0x2b ,0xe5 , 0xcc , 0xe8 , 0xc1 ,
        0xc9 , 0xe0 , 0xc0 , 0xe9 , 0xd5 , 0xf5 , 0xd8 , 0xf8 , 0xd0 , 0xf0 , 0xd9 , 0xf9 ,
        0xa5 , 0x1c , 0xa8 , 0x12 , 0x1b , 0xa0 , 0x13 , 0xa9 , 0x05 , 0xb5 , 0x0a , 0xb8 ,
        0x03 , 0xb0 , 0x0b , 0xb9 ,0x32 , 0x88 , 0x3c , 0x85 , 0x8d , 0x34 , 0x84 , 0x3d ,
        0x91 , 0x22 , 0x9c , 0x2c , 0x94 , 0x24 , 0x9d , 0x2d ,0x62 , 0x4a , 0x6c , 0x45 ,
        0x4d , 0x64 , 0x44 , 0x6d , 0x52 , 0x72 , 0x5c , 0x7c , 0x54 , 0x74 , 0x5d , 0x7d ,
        0xa1 , 0x1a , 0xac , 0x15 , 0x1d , 0xa4 , 0x14 , 0xad , 0x02 , 0xb1 , 0x0c , 0xbc ,
        0x04 , 0xb4 , 0x0d , 0xbd ,0xe1 , 0xc8 , 0xec , 0xc5 , 0xcd , 0xe4 , 0xc4 , 0xed ,
        0xd1 , 0xf1 , 0xdc , 0xfc , 0xd4 , 0xf4 , 0xdd , 0xfd ,0x36 , 0x8e , 0x38 , 0x82 ,
        0x8b , 0x30 , 0x83 , 0x39 , 0x96 , 0x26 , 0x9a , 0x28 , 0x93 , 0x20 , 0x9b , 0x29 ,
        0x66 , 0x4e , 0x68 , 0x41 , 0x49 , 0x60 , 0x40 , 0x69 , 0x56 , 0x76 , 0x58 , 0x78 ,
        0x50 , 0x70 , 0x59 , 0x79 ,0xa6 , 0x1e , 0xaa , 0x11 , 0x19 , 0xa3 , 0x10 , 0xab ,
        0x06 , 0xb6 , 0x08 , 0xba , 0x00 , 0xb3 , 0x09 , 0xbb ,0xe6 , 0xce , 0xea , 0xc2 ,
        0xcb , 0xe3 , 0xc3 , 0xeb , 0xd6 , 0xf6 , 0xda , 0xfa , 0xd3 , 0xf3 , 0xdb , 0xfb ,
        0x31 , 0x8a , 0x3e , 0x86 , 0x8f , 0x37 , 0x87 , 0x3f , 0x92 , 0x21 , 0x9e , 0x2e ,
        0x97 , 0x27 , 0x9f , 0x2f ,0x61 , 0x48 , 0x6e , 0x46 , 0x4f , 0x67 , 0x47 , 0x6f ,
        0x51 , 0x71 , 0x5e , 0x7e , 0x57 , 0x77 , 0x5f , 0x7f ,0xa2 , 0x18 , 0xae , 0x16 ,
        0x1f , 0xa7 , 0x17 , 0xaf , 0x01 , 0xb2 , 0x0e , 0xbe , 0x07 , 0xb7 , 0x0f , 0xbf ,
        0xe2 , 0xca , 0xee , 0xc6 , 0xcf , 0xe7 , 0xc7 , 0xef , 0xd2 , 0xf2 , 0xde , 0xfe ,
        0xd7 , 0xf7 , 0xdf , 0xff};

void SubCell8(unsigned char state[4][4])
{
    int i,j;
    for(i = 0; i < 4; i++)
        for(j = 0; j <  4; j++)
            state[i][j] = sbox_8[state[i][j]];
}

void ShiftRows(unsigned char state[4][4])
{
    int i, j, pos;

    unsigned char state_tmp[4][4];
    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            //application of the ShiftRows permutation
            pos=P[j+4*i];
            state_tmp[i][j]=state[pos>>2][pos&0x3];
        }
    }

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            state[i][j]=state_tmp[i][j];
        }
    }
}

void MixColumn(unsigned char state[4][4])
{
    int j;
    unsigned char temp;

    for(j = 0; j < 4; j++){
        state[1][j]^=state[2][j];
        state[2][j]^=state[0][j];
        state[3][j]^=state[2][j];

        temp=state[3][j];
        state[3][j]=state[2][j];
        state[2][j]=state[1][j];
        state[1][j]=state[0][j];
        state[0][j]=temp;
    }
}


void skinny_round_128(unsigned char* input, const unsigned char* key)
{
    unsigned char state[4][4];
    unsigned char keyCells[4][4];

    int l = 0 ;
    for(int i= 0; i<4; i++){
        for(int j = 0 ; j<4; j++){
            state[i][j] = input[l];
            keyCells[i][j] = key[l];
            l++;
        }
    }

    SubCell8(state);
    for(int i= 0; i<4; i++){
        for(int j = 0 ; j<4; j++){
            state[i][j] ^= keyCells[i][j];
        }
    }
    ShiftRows(state);
    MixColumn(state);

    l = 0 ;
    for(int i = 0; i<4; i++){
        for(int j = 0 ; j<4; j++){
            input[l] = state[i][j];
            l++;
        }
    }
}



void convert_to_bits(int *X, unsigned char*S, unsigned char*Z){
    unsigned char temp[16] ;

    for(int i = 0; i< 16; i++){
        temp[i] = (S[i] ^ Z[i]) ;
    }

    for(int i = 0; i < 16; i++){
        for(int j = 0; j<8; j++){
            X[8*i +j] = ( temp[i] >> (j) ) & 1 ;
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
    unsigned char S[16], Z[16], K[16];

    for(r = 1 ; r<rounds; r++){
        map< bitset<BS>, int, cmpBitset> countingBox;
        for(i = 0; i<steps; i++){
            for(j = 0 ; j< 16; j++){
                S[j] = rand() % 256;
                Z[j] = S[j] ^ A[j];
            }
            for(j = 0 ; j<r; j++){
                for(k = 0 ; k< 16; k++){
                    K[k] = rand() % 256 ;
                }
                skinny_round_128(S, K);
                skinny_round_128(Z, K);
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

    A[s] = subspace & 0xFF ;


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

    printf("#### SKINNY-128 ####\n");
    printf("#### First-order ####\n");
    avg_rank(12, 0x02 & 0xFF, 100, 7, 148, 1);
    printf("---------------------------------------------------------------------------\n");

    printf("#### Second-order ####\n");
    avg_rank(12, 0x02 & 0xFF, 100, 7, 8784, 2);
    printf("---------------------------------------------------------------------------\n");

    printf("#### Third-order ####\n");
    avg_rank(12, 0x02 & 0xFF, 100, 7, 17286, 3);
    printf("---------------------------------------------------------------------------\n");

}







