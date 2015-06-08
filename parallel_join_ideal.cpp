//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  parallel_join_ideal.cpp is part of 'FIIOFS'                             //
//  Copyright (C) 2014 Alex Bailey                                          //
//                                                                          //
//  Licensing information can be found in the README file                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <set>
#include <stdio.h>
#include <sstream>
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
  // n = index, k = support, i = orbit size, q = which job to run now

  static int n=atoi(a[1]); static int k=atoi(a[2]); long i=atol(a[3]); long q=atol(a[4]);

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

  // Finding relevant files in input folder

  set<data_block> all_data_blocks;
  data_block temp_data_block; temp_data_block.index=n, temp_data_block.bit=ibit; temp_data_block.data.resize((ibit/8)*(n-1));
  cout << "Reading from folder: " << folder_name << endl;
  cout << "Joining following files:" << endl;
  
  stringstream stream_required_prefix;
  stream_required_prefix << "n-" << n << "-k-" << k << "-i-" << i << ".bit" << ibit << "-";
  string required_prefix=stream_required_prefix.str();

  struct dirent **all_files;
  int x=scandir(folder_name.c_str(),&all_files,0,alphasort);
  size_t array_of_filesizes[x];
  vector<string> vector_of_filenames;
  int number_of_matching_files=0;
  if (x<0) perror("scandir");
  else {
    for (int z=0;z<x;z++) {
      char file_prefix[50]={0};
      memcpy(&file_prefix[0],all_files[z]->d_name,floor(log10(n))+floor(log10(k))+floor(log10(i))+floor(log10(ibit))+17);
      if (strcmp(required_prefix.c_str(),file_prefix)==0)
      {
        string file_name=folder_name+all_files[z]->d_name;
        ifstream file(file_name.c_str(),ios::in|ios::binary|ios::ate); streampos file_size=file.tellg();
        vector_of_filenames.push_back(file_name);
        array_of_filesizes[number_of_matching_files]=file_size;
        cout << "File " << number_of_matching_files << ": ";
        cout << "file name = " << all_files[z]->d_name << ", ";
        cout << "file size = " << file_size << ", ";
        cout << "number of semigroups = " << file_size/((n-1)*(ibit/8)) << endl;
        number_of_matching_files++;
        file.close();
      }
    }
  }
  free(all_files);

  // Reading in precomputed configuration file

  stringstream stream_config_file;
  stream_config_file << "join-config-n-" << n << "-k-" << k << "-i-" << i << ".bit" << ibit;
  string config_file=stream_config_file.str();

  ifstream file(config_file.c_str(),ios::in|ios::binary|ios::ate);
  streampos file_size=file.tellg();
  char *data_stream=new char[file_size];
  file.seekg(0,ios::beg);
  file.read(data_stream,file_size);
  file.close();

  unsigned long long p=file_size/(number_of_matching_files*sizeof(unsigned long long))-1;
  unsigned long long array_of_positions[(p+1)*number_of_matching_files];
  unsigned long long temp_position;
  for (int step=0;step<=p;step++) {
    for (int f=0;f<number_of_matching_files;f++) {
      temp_position=0;
      memcpy(&temp_position,&data_stream[(step*number_of_matching_files+f)*sizeof(unsigned long long)],sizeof(unsigned long long));
      array_of_positions[step*number_of_matching_files+f]=temp_position;
    }
  }
  delete[] data_stream;

  unsigned long long bStart;
  unsigned long long bEnd;
  unsigned long long number_of_blocks;
  char temp_semigroup[(n+1)*n];
  temp_semigroup[0]=1;
  temp_semigroup[1]=0;

  // Reading in relevant bytes specified by configuration file for this job

  cout << "Job " << q << endl << "---------------" << endl;
  for (int f=0;f<number_of_matching_files;f++) {
    bStart=array_of_positions[(q-1)*number_of_matching_files+f];
    bEnd=array_of_positions[q*number_of_matching_files+f];
    number_of_blocks=(bEnd-bStart)/((ibit/8)*(n-1));
    cout << "File " << f << ": reading " << number_of_blocks << " semigroups from bytes " << bStart << " to " << bEnd << endl;

    char *data_stream=new char[bEnd-bStart];
    ifstream file(vector_of_filenames[f].c_str(),ios::binary);
    file.seekg(bStart,ios::beg);
    file.read(data_stream,bEnd-bStart);
    file.close();

    for (int h=0;h<number_of_blocks;h++)
    {
      memcpy(&temp_data_block.data[0],&data_stream[(ibit/8)*(n-1)*h],(ibit/8)*(n-1));
      // Add semigroups from each file to set of all semigroups
      all_data_blocks.insert(temp_data_block);
    }
    delete[] data_stream;
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
  stream_output_file << "ideals/n-" << n << "/joined-n-" << n << "-k-" << k << "-i-" << i << ".bit" << ibit << "-";
  for (int j=floor(log10(q));j<floor(log10(p));j++) stream_output_file << "0";
  stream_output_file << q;
  string output_file=stream_output_file.str();
  cout << "Output file = " << output_file.c_str() << ", number of semigroups = " << num_of_semigroups << endl;
  ofstream ofile(output_file.c_str(),ios::out|ios::binary);
  ofile.write(&buff[0],(ibit/8)*(n-1)*num_of_semigroups);
  ofile.close();
  delete[] buff; 
return(0);
}
