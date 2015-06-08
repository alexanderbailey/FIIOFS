//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  join_ideal.cpp is part of 'FIIOFS'                                      //
//  Copyright (C) 2014 Alex Bailey                                          //
//                                                                          //
//  Licensing information can be found in the README file                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <dirent.h>

using namespace std;

// Function to encode a word in the free semigroup to its shortlex position

inline unsigned long long shortlex_encode(const char *V, int r) {
  unsigned long long s=0;
  if (r==1) { return(*V-1); }
  else {
    for (int i=1;i<=*V;i++) { s=s+*(V+i)*pow(r,*V-i); }
    return(s+(pow(r,*V)-r)/(r-1)); }
}

// Operator to compare two blocks of memory representing semigroups

inline bool shortlex_less_block (const char *u, const char *v, int n, int bit) {
  unsigned long long encoded_word1=0;
  unsigned long long encoded_word2=0;
  int l=0;
  while (encoded_word1==encoded_word2 && l<n-1)
  {
    encoded_word1=0;
    encoded_word2=0;
    memcpy(&encoded_word1,u+(bit/8)*l,bit/8);
    memcpy(&encoded_word2,v+(bit/8)*l,bit/8);
    l++;
  }
  if (encoded_word1 < encoded_word2) { return(1); } else { return(0); }
}

// Defining a structure for the block of memory that represents a semigroup

struct data_block {
  vector<char> data;
  int index;
  int bit;
  bool operator<(const data_block &rhs) const {
    if (shortlex_less_block(&data[0],&rhs.data[0],index,bit)) return (1); else return(0);
  }
};

// Main function begins here

int main(int c, char *a[]) {
  
  // Retrieving command line parameters
  // n = index, k = support, i = orbit size

  static int n=atoi(a[1]); static int k=atoi(a[2]);

  stringstream stream_folder_name;
  stream_folder_name << "ideals/n-" << n << "/";
  string folder_name=stream_folder_name.str();

  // Calculating bitrate of files

  int ibit;
  char largest_frob[n+1]; largest_frob[0]=n;
  {
    for (int j=1;j<=n;j++) largest_frob[j]=k-1;
    ibit=8*ceil((floor(log(shortlex_encode(&largest_frob[0],k))/log(2))+1)/8);  
  }

  // Finding relevant files in folder

  set<data_block> all_data_blocks;
  data_block temp_data_block; temp_data_block.index=n, temp_data_block.bit=ibit; temp_data_block.data.resize((ibit/8)*(n-1));
  cout << "Reading from folder: " << folder_name << endl;
  cout << "Index = " << n << ", support = " << k << endl;

  stringstream stream_intial_required_prefix;
  stream_intial_required_prefix << "n-" << n << "-k-" << k;
  string intial_required_prefix=stream_intial_required_prefix.str();

  set<long> orbit_sizes;

  struct dirent **all_files;
  int x=scandir(folder_name.c_str(),&all_files,0,alphasort);
  if (x<0) perror("scandir");
  else
  {
    for (int z=0;z<x;z++)
    {
      char file_prefix[50]={0};
      memcpy(&file_prefix[0],all_files[z]->d_name,floor(log10(n))+floor(log10(k))+7);
      if (strcmp(intial_required_prefix.c_str(),file_prefix)==0)
      {
        string file_name = folder_name + all_files[z]->d_name;
        ifstream file(file_name.c_str(),ios::in|ios::binary|ios::ate); streampos file_size=file.tellg();
        file.close();
        int pos_i=file_name.find("-i-");
        int pos_bit=file_name.find(".bit");
        // Checking all available orbit sizes from filenames
        if (strcmp(file_name.substr(pos_bit+floor(log10(ibit))+5,1).c_str(),"-")==0)
          { orbit_sizes.insert(atol(file_name.substr(pos_i+3,pos_bit-pos_i-3).c_str())); }
      }
    }
  }

  if (strcmp(a[3],"A")!=0) { orbit_sizes.clear(); orbit_sizes.insert(atol(a[3])); }

  set<long>::iterator i_iter_begin=orbit_sizes.begin();
  set<long>::iterator i_iter_end=orbit_sizes.end();
  set<long>::iterator i_iter=i_iter_begin;
  // For each file of a given orbit size...
  while (i_iter!=i_iter_end)
  {
    all_data_blocks.clear();
    cout << "Orbit size = " << *i_iter << endl;
    cout << "Joining following files:" << endl;

    stringstream stream_required_prefix;
    stream_required_prefix << "n-" << n << "-k-" << k << "-i-" << *i_iter << ".bit" << ibit << "-";
    string required_prefix=stream_required_prefix.str();

    for (int z=0;z<x;z++) {
      char file_prefix[50]={0};
      memcpy(&file_prefix[0],all_files[z]->d_name,floor(log10(n))+floor(log10(k))+floor(log10(*i_iter))+floor(log10(ibit))+17);
      if (strcmp(required_prefix.c_str(),file_prefix)==0)
      {
        string file_name=folder_name+all_files[z]->d_name;
        ifstream file(file_name.c_str(),ios::in|ios::binary|ios::ate); streampos file_size=file.tellg();
        unsigned long long number_of_blocks=file_size/((n-1)*(ibit/8));
        cout << "File name = " << all_files[z]->d_name << ", ";
        cout << "file size = " << file_size << ", ";
        cout << "number of semigroups = " << number_of_blocks << endl;
        char *data_stream=new char[file_size];
        file.seekg(0);
        file.read(data_stream,file_size);
        file.close();
        for (int h=0;h<number_of_blocks;h++)
        {
          memcpy(&temp_data_block.data[0],&data_stream[(ibit/8)*(n-1)*h],(ibit/8)*(n-1));
          // Add semigroups from files to set of all semigroups
          all_data_blocks.insert(temp_data_block);
        }
        delete[] data_stream;
      }
    }

    // Output joined file

    unsigned long long num_of_semigroups=all_data_blocks.size();
    char *buff=new char[(ibit/8)*(n-1)*(num_of_semigroups)];
    set<data_block>::iterator iter_begin=all_data_blocks.begin();
    set<data_block>::iterator iter_end=all_data_blocks.end(); 
    set<data_block>::iterator iter=iter_begin;
    size_t d=0;
    while (iter!=iter_end)
    {
      memcpy(&buff[(ibit/8)*(n-1)*d],&iter->data[0],(ibit/8)*(n-1));
      iter++; d++;
    }
    stringstream stream_output_file;
    stream_output_file << "ideals/n-" << n << "/n-" << n << "-k-" << k << "-i-" << *i_iter << ".bit" << ibit;
    string output_file=stream_output_file.str();
    cout << "Output file = " << output_file << ", number of semigroups = " << num_of_semigroups << endl;
    ofstream ofile(output_file.c_str(),ios::out|ios::binary);
    ofile.write(&buff[0],(ibit/8)*(n-1)*num_of_semigroups);
    ofile.close();
    delete[] buff;

    // Loop through all requested orbit sizes
    i_iter++;
  }

  free(all_files);

return(0);
}
