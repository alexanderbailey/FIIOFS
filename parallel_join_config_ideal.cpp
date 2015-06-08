//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  parallel_join_config_ideal.cpp is part of 'FIIOFS'                      //
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

inline unsigned long long shortlex_encode(const char *V, char r) {
  unsigned long long s=0;
  if (r==1) { return(*V-1); }
  else {
    for (int i=1;i<=*V;i++) { s=s+*(V+i)*pow(r,*V-i); }
    return(s+(pow(r,*V)-r)/(r-1)); }
}

// Operator to compare two blocks of memory representing semigroups

inline bool shortlex_less_block (const char *u, const char *v, char n, int bit) {
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
  // n = index, k = support, i = orbit size, p = how many jobs

  static int n=atoi(a[1]); static int k=atoi(a[2]); long i=atol(a[3]); int p=atol(a[4]);

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
  char temp_data_block[(ibit/8)*(n-1)];

  // Finding relevant files in input folder

  cout << "Reading from folder: " << folder_name << endl;
  cout << "Joining following files:" << endl;

  stringstream stream_required_prefix;
  stream_required_prefix << "n-" << n << "-k-" << k << "-i-" << i << ".bit" << ibit << "-";
  string required_prefix=stream_required_prefix.str();
  struct dirent **all_files;
  int x=scandir(folder_name.c_str(),&all_files,0,alphasort);
  unsigned long long array_of_filesizes[x];
  vector<string> vector_of_filenames;
  data_block sample_data_block;
  sample_data_block.index=n;
  sample_data_block.bit=ibit;
  sample_data_block.data.resize((ibit/8)*(n-1));
  set<data_block> sample_of_data_blocks;
  unsigned long long m;
  unsigned long long bStart;
  int number_of_matching_files=0;
  char *sample_data_stream=new char[(ibit/8)*(n-1)];
  if (x<0) perror("scandir"); else {
       for (int z=0;z<x;z++) {
        char file_prefix[50]={0};
        memcpy(&file_prefix[0],all_files[z]->d_name,floor(log10(n))+floor(log10(k))+floor(log10(i))+floor(log10(ibit))+17);
        if (strcmp(required_prefix.c_str(),file_prefix)==0)
        {
          string file_name=folder_name+all_files[z]->d_name;
          ifstream file(file_name.c_str(),ios::in|ios::binary|ios::ate);
          streampos file_size=file.tellg();
          // Take a sample of semigroups by reading 5 uniformally distributed semigroups from each file
          for (int q=0;q<5;q++)
          {
            m=file_size/(ibit/8)/(n-1);
            bStart=ceil((q*m)/5)*(ibit/8)*(n-1);
            file.seekg(bStart,ios::beg);
            file.read(sample_data_stream,(ibit/8)*(n-1));
            memcpy(&sample_data_block.data[0],&sample_data_stream[0],(ibit/8)*(n-1));
            sample_of_data_blocks.insert(sample_data_block);
          }
          file.close();
          vector_of_filenames.push_back(file_name);
          array_of_filesizes[number_of_matching_files]=file_size;
	  cout << "File " << number_of_matching_files << ": ";
	  cout << "file name = " << all_files[z]->d_name << ", ";
          cout << "file size = " << file_size << ", ";
	  cout << "number of semigroups = " << file_size/((n-1)*(ibit/8)) << endl;
	  number_of_matching_files++;
        }
      }
    }
  delete[] sample_data_stream;
  free(all_files);

  // Begin determining where to split up each file...

  unsigned long long array_of_positions[(p+1)*number_of_matching_files];
  for (int h=0;h<number_of_matching_files;h++) { array_of_positions[h]=0; }
  for (int h=0;h<number_of_matching_files;h++) { array_of_positions[p*number_of_matching_files+h]=array_of_filesizes[h]; }

  set<data_block>::iterator iter_begin=sample_of_data_blocks.begin();
  set<data_block>::iterator iter_end=sample_of_data_blocks.end(); 
  set<data_block>::iterator iter=iter_begin;
  unsigned long long d=0;
  unsigned long long iterator_placeholder;

  m=sample_of_data_blocks.size();
  cout << "Calculating..." << endl;
  for (int step=1;step<p;step++) {
    // Create breaking points from the sample of semigroups based on the number of jobs
    cout << "Job " << step << endl;
    iterator_placeholder=ceil((step*m)/p);
    while (d!=iterator_placeholder)
    {
      iter++; d++;
    }

    memcpy(&temp_data_block[0],&iter->data[0],(ibit/8)*(n-1));

    for (int f=0;f<number_of_matching_files;f++) {
      if (array_of_positions[(step-1)*number_of_matching_files+f]==array_of_filesizes[f])
      {
        array_of_positions[step*number_of_matching_files+f]=array_of_filesizes[f];
      } else {
        char *data_stream=new char[(ibit/8)*(n-1)];
        unsigned long long begin_chunk=array_of_positions[(step-1)*number_of_matching_files+f];
        unsigned long long end_chunk=array_of_filesizes[f];
        unsigned long long middle_chunk;
        // Divide and conquer each file until we come within 10 semigroups of our breaking point
        do
        {
          middle_chunk=(ibit/8)*(n-1)*ceil(ceil((begin_chunk+end_chunk)/2)/((ibit/8)*(n-1)));
          ifstream file(vector_of_filenames[f].c_str(),ios::binary);
          file.seekg(middle_chunk,ios::beg);
          file.read(data_stream,(ibit/8)*(n-1));
          file.close();
          int z=shortlex_less_block(&temp_data_block[0],&data_stream[0],n,ibit);
          if (z==1) end_chunk=middle_chunk; else begin_chunk=middle_chunk;
        } while (end_chunk-begin_chunk>(ibit/8)*(n-1)*10);
	delete[] data_stream;
        unsigned long long number_of_blocks=(end_chunk-begin_chunk)/((ibit/8)*(n-1));
        // With last <10 semigroups find closest semigroup in file to breaking point
        if (number_of_blocks>0) {
	  char *data_stream2=new char[end_chunk-begin_chunk];
          ifstream file(vector_of_filenames[f].c_str(),ios::binary);
          file.seekg(begin_chunk,ios::beg);
          file.read(data_stream2,end_chunk-begin_chunk);
          file.close();
          int h=0;
          while (shortlex_less_block(&data_stream2[(ibit/8)*(n-1)*h],&temp_data_block[0],n,ibit)==1 and h<number_of_blocks-1)
          {
            h++;
          };
          // Add the byte position of this semigroup to an array of all positions
          if (shortlex_less_block(&temp_data_block[0],&data_stream2[(ibit/8)*(n-1)*h],n,ibit)==0)
          { array_of_positions[step*number_of_matching_files+f]=begin_chunk+(h+1)*(n-1)*(ibit/8); }
          else
          { array_of_positions[step*number_of_matching_files+f]=begin_chunk+h*(n-1)*(ibit/8); }
          delete[] data_stream2;
        }
        else
        {
          array_of_positions[step*number_of_matching_files+f]=begin_chunk;
        }
      }
    }
  }

  stringstream output1;
  stringstream output2;
  set<size_t> memory_required;
  size_t temp_memory;

  output1 << "Splitting files as follows:" << endl;
  output2 << endl << "Memory required:" << endl;
  for (int step=1;step<=p;step++) {
    output1 << endl << "Job " << step << endl;
    output2 << endl << "Job " << step << ": ";
    temp_memory=0;
    for (int f=0;f<number_of_matching_files;f++)
    {
      // Output which parts of which files will be read for each job
      temp_memory=temp_memory+array_of_positions[step*number_of_matching_files+f]-array_of_positions[(step-1)*number_of_matching_files+f];
      output1 << "File " << f << ": ";
      output1 << (array_of_positions[step*number_of_matching_files+f]-array_of_positions[(step-1)*number_of_matching_files+f])/((n-1)*(ibit/8)) << " semigroups";
      output1 << " from bytes " << array_of_positions[(step-1)*number_of_matching_files+f];
      output1 << " to " << array_of_positions[step*number_of_matching_files+f] << endl;
    }
    // Output how much memory is required for each job
    memory_required.insert(temp_memory);
    output2 << temp_memory;
  }

  cout << output1.str();
  cout << output2.str() << endl;
  // Output the maximum memory required out of all the jobs
  cout << "Maximum memory required:" << *memory_required.rbegin() << endl;

  char *buff=new char[(p+1)*number_of_matching_files*sizeof(unsigned long long)];
  memcpy(&buff[0],&array_of_positions[0],(p+1)*number_of_matching_files*sizeof(unsigned long long));

  // Save the configuration file with the byte positions

  stringstream stream_output_file;
  stream_output_file << "join-config-n-" << n << "-k-" << k << "-i-" << i << ".bit" << ibit;
  string output_file=stream_output_file.str();
  ofstream ofile(output_file.c_str(),ios::out|ios::binary);
  ofile.write(&buff[0],(p+1)*number_of_matching_files*sizeof(unsigned long long));
  ofile.close();
  delete[] buff; 
return(0);
}
