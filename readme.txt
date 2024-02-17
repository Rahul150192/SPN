**** Dependencies ****
1. GMP [https://gmplib.org/]
2. NTL [https://libntl.org/]

**** Running Codes ****
1. g++ -g -O2 -std=c++11 -pthread -march=native gift64.cpp -o foo -lntl -lgmp -lm
2. ./foo

The current run will produce output as follows:

#### GIFT-64 ####
#### First-order ####
Round    Avg. Dimension          Final Dimension         Avg. Unique Differences
0             0.00               64                          0.00
1             3.00               64                          6.00
2            11.00               64                         42.96
3            44.88               64                         81.48
---------------------------------------------------------------------------
#### Second-order ####
Round    Avg. Dimension          Final Dimension         Avg. Unique Differences
0             0.00               2080                        0.00
1             6.00               2080                        6.00
2            55.00               2080                      157.58
3           953.55               2080                     1567.78
---------------------------------------------------------------------------
#### Third-order ####
Round    Avg. Dimension          Final Dimension         Avg. Unique Differences
0             0.00               4033                        0.00
1             6.00               4033                        6.00
2            57.00               4033                      216.65
3          1246.50               4033                     2733.80
4          3354.20               4033                     4346.59
5          4016.74               4033                     4366.96
---------------------------------------------------------------------------

**** Changing ciphers ****
In the command: g++ -g -O2 -std=c++11 -pthread -march=native gift64.cpp -o foo -lntl -lgmp -lm
change gift64.cpp to
1. gift128.cpp
2. present.cpp
3. skinny64.cpp
4. skinny128.cpp




