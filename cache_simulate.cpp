#include <string>
#include <iostream>
#include <sstream>

using namespace std;

int hexToDec(string hexStr) {
    int decNum;
    stringstream ss;
    ss << hex << hexStr;
    ss >> decNum;
    return decNum;
}

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <bitset>
#include <math.h>
#include <cstring>
#include <algorithm>
#include <tuple>
#include <vector>


using namespace std;
//access state:
#define NA 0 // no action
#define RH 1 // read hit
#define RM 2 // read miss
#define WH 3 // Write hit
#define WM 4 // write miss

struct config{
       int L1blocksize;
       int L1assoc;       // associat * block
       int L1size;
       int L2blocksize;
       int L2assoc;
       int L2size;
       };


class cache {
    public:
    int n_sets;
    int n_index_bits;
    int n_offset_bits;
    int n_tag_bits;
    
    int tag_match;
    bool valid_bit;

    // vector<int> offset;
    // vector<int> tags;
    // vector<int> index;
    // vector<int> dirty;

    unsigned long tag_int;
    unsigned long index_int;
    
    bool read_miss;
    
    int num_reads = 0;
    int num_read_misses = 0;
    int num_writes = 0;
    int num_write_misses = 0;
    int num_writebacks = 0;
    float miss_rate = 0;
    // vector < vector < unsigned long > > Mylru;

    // int lru[][];
    std::vector<std::vector<int>>lru;
    
    void computeN_Bits(int c_size,int associativity, int b_size){
        if(associativity>0){
            n_sets=(c_size)/(b_size*associativity);
            n_index_bits=log2(n_sets);
        }

        n_offset_bits = log2(b_size);
        n_tag_bits=64 - n_offset_bits - n_index_bits;
    
        lru = std::vector<std::vector<int>>(n_sets, std::vector<int>(associativity, 0));
    }  
};


std::tuple< int , int > compute_tag_index( int address , int no_of_sets,int block_size ){
    int addr = address / block_size;
    int index_int = addr % no_of_sets;
    int tag_int = addr / no_of_sets;
    return make_tuple(tag_int,index_int);
}

int address_from_tag_index(int tag , int no_of_sets , int index, int block_size){
    int address = ((no_of_sets * tag) + index ) * block_size   ;      //  ((no_of sets * tag) + index ) * block_size
    return address;
}

int main(int argc, char** argv) {

    
    if (argc != 7) {
        cerr << "Usage: " << argv[0] << " BLOCKSIZE L1_SIZE L1_ASSOC L2_SIZE L2_ASSOC TRACE_FILE" << endl;
        exit(EXIT_FAILURE);
    }

    // Parse input parameters
    int block_size = stoi(argv[1]);
    int l1_size = stoi(argv[2]);
    int l1_assoc = stoi(argv[3]);
    int l2_size = stoi(argv[4]);
    int l2_assoc = stoi(argv[5]);
    string filename = argv[6];
    
    
    
    
//     int mem_inst[1024];
// 	int index = 0;

    // int block_size = 64;
    // int l1_size = 1024;
    // int l1_assoc = 2;
    // int l2_size = 65536;
    // int l2_assoc = 8;
    
    config block;
    block.L1blocksize = block_size;
    block.L1assoc = l1_assoc ;
    block.L1size = l1_size ;
    block.L2blocksize = block_size;
    block.L2assoc = l2_assoc;
    block.L2size = l2_size;


    cache l1,l2 ;
    l1.computeN_Bits(block.L1size, block.L1assoc , block.L1blocksize);
    l2.computeN_Bits(block.L2size, block.L2assoc , block.L2blocksize);

    int L1block[l1.n_sets][3*block.L1assoc];  // tag bits, valid_bit, dirty bit
    int L2block[l2.n_sets][3*block.L2assoc];  
    

    for (int i = 0; i < l1.n_sets; i++) {
        for (int j = 0; j < 3*block.L1assoc; j++) {
            L1block[i][j] = 0;
        }
    }
    for (int i = 0; i < l2.n_sets; i++) {
        for (int j = 0; j < 3*block.L2assoc; j++) {
            L2block[i][j] = 0;
        }
    }

    int L1AcceState =0; // L1 access state variable, can be one of NA, RH, RM, WH, WM;
    int L2AcceState =0; // L2 access state variable, can be one of NA, RH, RM, WH, WM;

    ifstream traces;
    ofstream tracesout;
    // string outname;
    // outname = string("trace_small.txt") + ".out";

    int c1, c2;
    string access;
    
    // string filename = "trace1.txt" ;
    

    ifstream infile(filename);

    vector<string> operations;
    vector<string> addresses;

    string op;
    string addr;

    // int i = 0;
    

    while (infile >> op >> addr) {
        
        l1.tag_match = 0;
        l2.tag_match = 0;
        l1.valid_bit = 0;
        l2.valid_bit = 0;
        
        operations.push_back(op);
        addresses.push_back(addr);
        
        access = op;

        int addresses_int = hexToDec(addr);
 
        tie(l1.tag_int,l1.index_int) = compute_tag_index( addresses_int , l1.n_sets, block.L1blocksize);
        tie(l2.tag_int,l2.index_int) = compute_tag_index( addresses_int , l2.n_sets, block.L2blocksize);
        
        
        //////////////////////          L1 Read      ///////////////////////////////

        if (access == "r"){
            for(int i=0;i<3*block.L1assoc;i+=3){ //checking the tag bits    
                if(L1block[l1.index_int][i]==l1.tag_int){
                    l1.tag_match=1;
                    c1=i;
                    break;
                    }
                };
            if(l1.tag_match==1){                   //checking for a tag match
                if(L1block[l1.index_int][c1+1]==1){ // second index = valid bit
                    l1.valid_bit=1;
                }
            }    
            
            if(l1.valid_bit==1 && l1.tag_match==1 ){
                // cout << "Read Hit\n";
                L1AcceState=RH;
                l1.num_reads += 1;
                
                auto max_it = max_element(l1.lru[l1.index_int].begin(), l1.lru[l1.index_int].end());
                int max_val = *max_it;
                l1.lru[l1.index_int][c1/3] = max_val + 1; 
            
            }else {
                // cout << "Read Miss\n";
                L1AcceState=RM;
                l1.num_read_misses += 1;
                auto min_it = std::min_element(l1.lru[l1.index_int].begin(), l1.lru[l1.index_int].end());
                int min_index = std::distance(l1.lru[l1.index_int].begin(), min_it);
                if (L1block[l1.index_int][3*min_index + 1] == 1) {  // if evicting a block from l1 
                    if (L1block[l1.index_int][3*min_index + 2] == 1) {  // eviction from L1 and write bit is 1                        
                        
                        int evict_address = address_from_tag_index(L1block[l1.index_int][3*min_index] , l1.n_sets , l1.index_int, block.L1blocksize);
                        int tag_evict,index_evict;
                        tie(tag_evict,index_evict) = compute_tag_index(evict_address,l2.n_sets, block.L2blocksize);
                        
                        for(int j=0;j<3*block.L2assoc;j+=3){    //checking the tag bits
                            if(L2block[index_evict][j]==tag_evict) {
                                L2block[index_evict][j+1] = 1;     // valid
                                L2block[index_evict][j+2] = 1;     // dirty
                                // l2.num_writes += 1;
                                break;
                                }
                            };                   
                        l1.num_writebacks += 1;
                        l2.num_writes += 1;

                    }
                }
                L1block[l1.index_int][3*min_index] = l1.tag_int;   // changing tag bit
                L1block[l1.index_int][3*min_index + 1] = 1;        // setting valid bit
                L1block[l1.index_int][3*min_index + 2] = 0;        // setting dirty bit

                auto max_it = max_element(l1.lru[l1.index_int].begin(), l1.lru[l1.index_int].end());
                int max_val = *max_it;
                l1.lru[l1.index_int][min_index] = max_val + 1;
                l1.num_reads += 1;
                
            };


//////////////////////          L2 Read      ///////////////////////////////
            
            if (L1AcceState==RH){
                L2AcceState=NA;
            } else if (L1AcceState==RM) {
                for(int i=0;i<3*block.L2assoc;i+=3){    //checking the tag bits
                    if(L2block[l2.index_int][i]==l2.tag_int) {
                        l2.tag_match=1;
                        c2=i;
                        break;
                        }
                    };
        
                if(l2.tag_match==1){
                    if(L2block[l2.index_int][c2+1]==1){
                        l2.valid_bit=1;
                        }
                    }
        
                if(l2.valid_bit==1 && l2.tag_match==1){
                    // cout << "L2 Read Hit\n";
                    L2AcceState=RH;
                    l2.num_reads += 1;

                    auto max_it = max_element(l2.lru[l2.index_int].begin(), l2.lru[l2.index_int].end());
                    int max_val = *max_it;
                    l2.lru[l2.index_int][c2/3] = max_val + 1;
                
                } else if(l1.valid_bit!=1 || l1.tag_match!=1 ) {
                    // cout << "L2 Read Miss\n";
                    L2AcceState=RM;
                    l2.num_read_misses += 1;
                    
                    auto min_it = min_element(l2.lru[l2.index_int].begin(), l2.lru[l2.index_int].end());
                    int min_index = distance(l2.lru[l2.index_int].begin(), min_it);
                    

                    if (L2block[l2.index_int][3*min_index + 1] == 1) {  // if evicting a block from l2
                        if (L2block[l2.index_int][3*min_index + 2] == 1) {  // eviction from L1 and write bit is 1
                            l2.num_writebacks += 1;
                        }

                        int evict_address = address_from_tag_index(L2block[l2.index_int][3*min_index] , l2.n_sets , l2.index_int, block.L2blocksize);
                        int tag_evict,index_evict;
                        tie(tag_evict,index_evict) = compute_tag_index(evict_address,l1.n_sets, block.L1blocksize);
                        
                        for(int j=0;j<3*block.L1assoc;j+=3){    //checking the tag bits
                            if(L1block[index_evict][j]==tag_evict) {

                                if (L1block[index_evict][j+1]==1){
                                    L1block[index_evict][j+1] = 0;     // unvalid
                                    l1.lru[index_evict][j/3] = 0;      // taaki ab isko use kre for eviction
                                    if (L1block[index_evict][j+2]==1){
                                        l1.num_writebacks += 1;
                                        l1.num_writes += 1;
                                        l2.num_writebacks += 1;
                                    }
                                }
                                break;

                                };
                            } 

                    }
                    
                    
                    L2block[l2.index_int][3*min_index] = l2.tag_int;   // changing tag bit
                    L2block[l2.index_int][3*min_index + 1] = 1;        // setting valid bit
                    L2block[l2.index_int][3*min_index + 2] = 0;        // setting dirty bit
                    l2.num_reads += 1;
                    
                    auto max_it = max_element(l2.lru[l2.index_int].begin(), l2.lru[l2.index_int].end());
                    int max_val = *max_it;
                    l2.lru[l2.index_int][min_index] = max_val + 1;
                    
                    
                }
            }
        }
        else if (access == "w") {
    ///////////////////////     L1 Write     ////////////
            // cout << "WRITE" << endl; 
               
            for(int i=0;i<3*block.L1assoc;i+=3){ //checking the tag bits
                if(L1block[l1.index_int][i]==l1.tag_int){
                    l1.tag_match=1;
                    c1=i;
                    break;
                    }
                };

            if(l1.tag_match==1){
                if(L1block[l1.index_int][c1+1]==1){
                    l1.valid_bit=1;
                }
            };
                
            if(l1.valid_bit==1 && l1.tag_match==1 ) {
                // cout << "Write Hit\n";
                L1AcceState=WH;
                l1.num_writes += 1;

                L1block[l1.index_int][c1 + 2] = 1;     //change dirty bit
                  
                auto max_it = max_element(l1.lru[l1.index_int].begin(), l1.lru[l1.index_int].end());
                int max_val = *max_it;
                l1.lru[l1.index_int][c1/3] = max_val + 1;
                
            } else{
                // cout << "Write Miss\n";
                
                L1AcceState = WM;
                l1.num_write_misses += 1;

                auto min_it = std::min_element(l1.lru[l1.index_int].begin(), l1.lru[l1.index_int].end());
                int min_index = std::distance(l1.lru[l1.index_int].begin(), min_it);
                
                if (L1block[l1.index_int][3*min_index + 1] == 1) {  // if evicting a block from l1                                        
                    if (L1block[l1.index_int][3*min_index + 2] == 1) {  // eviction from L1 and write bit is 1                        
                        
                        int evict_address = address_from_tag_index(L1block[l1.index_int][3*min_index] , l1.n_sets , l1.index_int, block.L1blocksize);
                        int tag_evict,index_evict;
                        tie(tag_evict,index_evict) = compute_tag_index(evict_address,l2.n_sets, block.L2blocksize);
                        
                        for(int j=0;j<3*block.L2assoc;j+=3) {    //checking the tag bits
                            if(L2block[index_evict][j]==tag_evict) {
                                L2block[index_evict][j+1] = 1;     // valid
                                L2block[index_evict][j+2] = 1;     // dirty
                                // l2.num_writes += 1;
                                // l1.num_writebacks += 1;
                                break;
                                }
                            };                    
                        l1.num_writebacks += 1;
                        l2.num_writes += 1;

                    }
                }
                    
                L1block[l1.index_int][3*min_index] = l1.tag_int;   // changing tag bit
                L1block[l1.index_int][3*min_index + 1] = 1;        // setting valid bit
                L1block[l1.index_int][3*min_index + 2] = 1;        // setting dirty bit
                l1.num_writes += 1;
                
                auto max_it = max_element(l1.lru[l1.index_int].begin(), l1.lru[l1.index_int].end());
                int max_val = *max_it;
                l1.lru[l1.index_int][min_index] = max_val + 1;
                
                

            };
            //////////////////// L2 Write  //////////////////////////////////////////////////    
            if (L1AcceState == WH) {
                L2AcceState = NA;
            } else if (L1AcceState == WM) {
                
                for(int i=0;i<3*block.L2assoc;i+=3) { //checking the tag bits
                    if(L2block[l2.index_int][i]==l2.tag_int) {
                        l2.tag_match=1;
                        c2=i;
                        break;
                    }
                };
    
                if(l2.tag_match==1){
                    if(L2block[l2.index_int][c2+1]==1){
                        l2.valid_bit=1;
                    }
                }
                
                
                if(l2.valid_bit==1 && l2.tag_match==1 ) {
                    
                    // cout << "L2 Write Hit\n";
                    L2AcceState=WH;
                    l2.num_writes += 1;

                    auto max_it = max_element(l2.lru[l2.index_int].begin(), l2.lru[l2.index_int].end());
                    int max_val = *max_it;
                    l2.lru[l2.index_int][c2/3] = max_val + 1;


                } else {
                    
                    // cout << "L2 Write Miss\n";
                    L2AcceState=WM;
                    l2.num_write_misses += 1;
                    
                    auto min_it = min_element(l2.lru[l2.index_int].begin(), l2.lru[l2.index_int].end());
                    int min_index = distance(l2.lru[l2.index_int].begin(), min_it);
                    
                    if (L2block[l2.index_int][3*min_index + 1] == 1) {  // if evicting a block from l2
                        if (L2block[l2.index_int][3*min_index + 2] == 1) {  // eviction from L1 and write bit is 1
                            l2.num_writebacks += 1;
                        }

                        int evict_address = address_from_tag_index(L2block[l2.index_int][3*min_index] , l2.n_sets , l2.index_int, block.L1blocksize);
                        int tag_evict,index_evict;
                        tie(tag_evict,index_evict) = compute_tag_index(evict_address,l1.n_sets, block.L1blocksize);
                        
                        for(int j=0;j<3*block.L1assoc;j+=3){    //checking the tag bits
                            if(L1block[index_evict][j]==tag_evict) {

                                if (L1block[index_evict][j+1]==1){
                                    L1block[index_evict][j+1] = 0;     // unvalid
                                    l1.lru[index_evict][j/3] = 0;      // taaki ab isko use kre for eviction
                                    if (L1block[index_evict][j+2]==1){
                                        l1.num_writebacks += 1;
                                        l2.num_writes += 1;
                                        l2.num_writebacks += 1;
                                    }
                                }
                                break;

                                };
                            }  
                    }

                    L2block[l2.index_int][3*min_index] = l2.tag_int;   // changing tag bit
                    L2block[l2.index_int][3*min_index + 1] = 1;        // setting valid bit
                    L2block[l2.index_int][3*min_index + 2] = 0;        // setting dirty bit
                    l2.num_writes += 1;

                    auto max_it = max_element(l2.lru[l2.index_int].begin(), l2.lru[l2.index_int].end());
                    int max_val = *max_it;
                    l2.lru[l2.index_int][min_index] = max_val + 1;
                    
                    
                }
                
                
            }


        } // end access = "w"
    } // while end
    
    // cout << "L1 reads: " << l1.num_reads << endl;
    // cout << "L1 read misses: " << l1.num_read_misses << endl;
    // cout << "L1 writes: " << l1.num_writes << endl;
    // cout << "L1 write misses: " << l1.num_write_misses << endl;
    // cout << "L1 writebacks: " << l1.num_writebacks << endl;
    // cout << "L2 reads: " << l2.num_reads << endl;
    // cout << "L2 read misses: " << l2.num_read_misses << endl;
    // cout << "L2 writes: " << l2.num_writes << endl;
    // cout << "L2 write misses: " << l2.num_write_misses << endl;
    // cout << "MEM writebacks: " << l2.num_writebacks << endl;    
    
    l1.miss_rate = (float)(l1.num_read_misses + l1.num_write_misses) / (l1.num_reads + l1.num_writes);
    l2.miss_rate = (float)(l2.num_write_misses + l2.num_read_misses) / (l2.num_reads + l2.num_writes);
    
    cout << "\n\n===== Simulation Results =====";
    cout << endl;

    cout << "\ni. number of L1 reads:\t\t\t\t" << dec << l1.num_reads;
    cout << "\nii. number of L1 read misses:\t\t\t" << dec << l1.num_read_misses;
    cout << "\niii. number of L1 writes:\t\t\t" << dec << l1.num_writes;
    cout << "\niv. number of L1 write misses:\t\t\t" << dec << l1.num_write_misses;
    cout << "\nv. L1 miss rate:\t\t\t\t" << fixed << setprecision(4) << l1.miss_rate;
    cout << "\nvi. number of writebacks from L1 memory:\t" << dec << l1.num_writebacks;
    cout << endl;
    if (block.L2blocksize != 0)
    cout << "\nvii. number of L2 reads:\t\t\t" << dec << l2.num_reads;
    cout << "\nviii. number of L2 read misses:\t\t\t" << dec << l2.num_read_misses;
    cout << "\nix. number of L2 writes:\t\t\t" << dec << l2.num_writes;
    cout << "\nx. number of L2 write misses:\t\t\t" << dec << l2.num_write_misses;
    cout << "\nxi. L2 miss rate:\t\t\t\t" << fixed << setprecision(4) << l2.miss_rate;
    cout << "\nxii. number of writebacks from L2 memory:\t" << dec << l2.num_writebacks << "\n";
    cout << endl;
};
