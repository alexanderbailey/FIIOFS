//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  next_index_ideal.cpp is part of 'FIIOFS'                                //
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
#include <algorithm>
#include <dirent.h>

using namespace std;

// Function to create array of all perumataions of Sym(k)

void perm(int k, char *A) {
  char id[k]; for (int i=0;i<k;i++) {id[i]=i; }
  unsigned int c=0;
  do {
    for (int i=0;i<k;i++) { *(A+c*k+i)=id[i]; } c++;
  } while (next_permutation(id,id+k));
}

// Factorial function

unsigned long long factorial(unsigned int n) { return (n==1||n==0)?1:factorial(n-1)*n; }

// Function to encode a word in the free semigroup to its shortlex position

inline unsigned long long shortlex_encode(const char *V, int r) {
  unsigned long long s=0;
  if (r==1) { return (*V-1); }
  else { for (int i=1;i<=*V;i++) { s=s+*(V+i)*pow(r,*V-i); }
    return(s+(pow(r,*V)-r)/(r-1)); }
}

// Function to decode a word from its shortlex position

inline void shortlex_decode(unsigned long long x, int r, char *V) {
  if (r==1) { *V=x+1; for (int i=1;i<=x+1;i++) *(V+i)=0; }
  else { int m=floor(log((r-1)*(x+2))/log(r)+0.00000000001);
    *V=m;
    if(m==1) { *(V+1)=x; }
    else
    {
      unsigned long long t=x+1-(pow(r,m)-1)/(r-1);
      for (int i=1;i<=m-1;i++)
      {
        *(V+i)=floor(t/pow(r,m-i));
        t=t-pow(r,m-i)*(*(V+i));
      }
      *(V+m)=t;
    }
  }
}

// Operator to compare words by their shortlex position

inline bool shortlex_less_word (const char *u, const char *v) {
  int i=0; while (*(u+i)==*(v+i) && i<*u) { i++; }
  if (*(u+i)<*(v+i)) { return(1); } else { return(0); }
}

// Defining a structure for a word in the free semigroup

struct word {
  vector<char> gaps;
  bool operator<(const word &rhs) const { return(shortlex_less_word(&gaps[0],&rhs.gaps[0])); }
}; 

// Defining a subsemigroup structure which contains the index and an array of the words in the gaps

struct semigroup {
  int index;
  vector<char> gaps;
  bool operator<(const semigroup &rhs) const {
    int k=0; while (k<index) {
      if (shortlex_less_word(&gaps[k*(index+1)],&rhs.gaps[k*(index+1)])) { return(1); } else {
        if (shortlex_less_word(&rhs.gaps[k*(index+1)],&gaps[k*(index+1)])) { return(0); } else { k++; } } };
    return(0);
  }
};

// Defining another subsemigroup structure that also includes the orbit size

struct orbit_semigroup {
  int index;
  long orbit;
  vector<char> gaps;
  bool operator<(const orbit_semigroup &rhs) const {
    if (orbit<rhs.orbit) { return (1); } else { if (rhs.orbit<orbit) { return(0); } else {
    int k=0; while (k<index) {
      if (shortlex_less_word(&gaps[k*(index+1)],&rhs.gaps[k*(index+1)])) { return(1); } else {
        if (shortlex_less_word(&rhs.gaps[k*(index+1)],&gaps[k*(index+1)])) { return(0); } else { k++; } } };
    return(0);
  } } }
};

// Function that orders the words in the gaps after a permutation has been applied

void order_semigroup(char *ordered_semigroup, char *semigroup, int n) {
  char *words[n];
  for (int k=0;k<=n-1;k++) { words[k]=semigroup+(n+1)*k; }
  sort(words,words+n,shortlex_less_word);
  for (int k=0;k<=n-1;k++) { memcpy(ordered_semigroup+(n+1)*k,words[k],n+1); }
}

// Function that reads semigroups from a file in to a set

void read_file(set<semigroup> & input_semigroups, string filename, int bitrate, int n, int k, int byteStart, int byteEnd) {
  size_t semigroupStart, semigroupEnd;
  unsigned long long encoded_word=0;
  semigroup temp_semigroup;
  temp_semigroup.index=n;
  temp_semigroup.gaps.resize((n+1)*n);
  temp_semigroup.gaps[0]=1;
  temp_semigroup.gaps[1]=0;
  ifstream file(filename.c_str(),ios::in|ios::binary);
  char *data_stream=new char[byteEnd-byteStart];
  file.seekg(byteStart,ios::beg);
  file.read(data_stream,byteEnd-byteStart);
  file.close();
  semigroupStart=byteStart/(bitrate/8)/(n-1);
  semigroupEnd=byteEnd/(bitrate/8)/(n-1);
  for (int j=0;j<=(semigroupEnd-semigroupStart)-1;j++)
  {
    for (int l=0;l<n-1;l++)
    {
      encoded_word=0;
      memcpy(&encoded_word,&data_stream[(bitrate/8)*((n-1)*j+l)],bitrate/8);
      shortlex_decode(encoded_word,k,&(temp_semigroup.gaps[(n+1)*(l+1)]));
    }
    input_semigroups.insert(temp_semigroup);
  }
  delete[] data_stream;
}

// Main function begins here

int main(int c, char *a[]) {

  // Retrieving command line parameters
  // n = index, k = support, p = number of parallel jobs, i = which job to run now

  static int n=atoi(a[1]); static int k=atoi(a[2]); long p=atol(a[3]); long i=atol(a[4]);

  stringstream stream_folder_name, stream_output_folder_name;
  stream_folder_name << "ideals/n-" << n-1 << "/";
  stream_output_folder_name << "ideals/n-" << n << "/";
  string folder_name=stream_folder_name.str();
  string output_folder_name=stream_output_folder_name.str();

  // Check if directory exists

  struct dirent **all_files;
  if (scandir(folder_name.c_str(),&all_files,0,alphasort)<0)
  {
    cout << "Input directory '" << folder_name << "' does not exist." << endl;
    exit(0);
  }
  if (scandir(output_folder_name.c_str(),&all_files,0,alphasort)<0)
  {
    cout << "Output directory '" << output_folder_name << "' does not exist. Create it first." << endl;
    exit(0);
  }

  // Calculating bitrate of input and output files

  int ibit;
  int obit;
  char largest_frob[n+1];
  largest_frob[0]=n-1;
  for (int j=1;j<=n-1;j++) largest_frob[j]=k-1;
  ibit=8*ceil((floor(log(shortlex_encode(&largest_frob[0],k))/log(2))+1)/8);
  largest_frob[0]=n;
  for (int j=1;j<=n;j++) largest_frob[j]=k-1;
  obit=8*ceil((floor(log(shortlex_encode(&largest_frob[0],k))/log(2))+1)/8);

  // Declaring variables

  set<semigroup> input_semigroups;
  set<semigroup> orbit;
  set<semigroup>::iterator iter_begin;
  set<semigroup>::iterator iter_end;
  set<orbit_semigroup> output_semigroups;

  set<word>::iterator iter_frob;
  set<word> new_words;
  set<word> gaps;
  
  unsigned long long count;

  semigroup temp_semigroup;
  temp_semigroup.index=n;
  temp_semigroup.gaps.resize((n+1)*n);
  temp_semigroup.gaps[0]=1;
  temp_semigroup.gaps[1]=0;

  word temp_word; temp_word.gaps.resize(n+1); 
  word temp_word2; temp_word2.gaps.resize(n+1);
  word right_word; right_word.gaps.resize(n);  

  bool fail;

  semigroup ordered_semigroup;
  ordered_semigroup.index=n;
  ordered_semigroup.gaps.resize((n+1)*n);

  orbit_semigroup temp_output_semigroup;
  temp_output_semigroup.index=n;
  temp_output_semigroup.gaps.resize((n+1)*n);

  // If k=1 then add the only possible new ideal to output_semigroups

  if (k==1)
  {
    for (int l=0;l<n;l++)
    {
      temp_output_semigroup.gaps[(n+1)*l]=l+1;
      for (int m=0;m<l;m++)
      {
        temp_output_semigroup.gaps[(n+1)*l+m+1]=0;
      }
    }
    temp_output_semigroup.orbit=1;
    output_semigroups.insert(temp_output_semigroup);
  }

  // If n=k then add the only possible new ideal to output_semigroups

  if (n==k)
  {
    for (int l=0;l<n;l++)
    {
      temp_output_semigroup.gaps[(n+1)*l]=1;
      temp_output_semigroup.gaps[(n+1)*l+1]=l;
    }
    temp_output_semigroup.orbit=1;
    output_semigroups.insert(temp_output_semigroup);
  }

  // Otherwise calculate new ideals from previous files

  if (k!=n and k!=1) 
  {
    unsigned long long f=factorial(k); char *A=new char[k*f]; perm(k,&A[0]);
    char temp_perm_semigroup[(n+1)*n];

    cout << "Running job " << i << " of " << p << endl;
    cout << "Reading from folder: " << folder_name << endl;
  
    // Finding relevant files in input folder

    cout << "All files with index " << n-1 << " and support " << k << endl;
    cout << "----------------------------------------------------------------------" << endl;
 
    size_t total_file_size=0;
      
    stringstream stream_required_prefix;
    stream_required_prefix << "n-" << n-1 << "-k-" << k << "-i-";
    string required_prefix=stream_required_prefix.str();
	  
    int x=scandir(folder_name.c_str(),&all_files,0,alphasort);
    size_t array_of_filesizes[x];
    vector<string> vector_of_filenames;
    int number_of_matching_files=0;
    if (x<0) perror("scandir"); else {
      char file_name[100];
      for (int z=0;z<x;z++) {
        char file_prefix[50]={0};
        memcpy(&file_prefix[0],all_files[z]->d_name,floor(log10(n-1))+floor(log10(k))+10);
        if (strcmp(required_prefix.c_str(),file_prefix)==0)
        {
          string file_name = folder_name + all_files[z]->d_name;
          ifstream file(file_name.c_str(),ios::in|ios::binary|ios::ate);
          vector_of_filenames.push_back(file_name);
          array_of_filesizes[number_of_matching_files]=file.tellg();
          file.close();
          total_file_size=total_file_size+array_of_filesizes[number_of_matching_files];
          cout << "File " << number_of_matching_files << ": ";
          cout << "file name = " << all_files[z]->d_name << ", ";
          cout << "file size = " << array_of_filesizes[number_of_matching_files] << ", ";
          cout << "number of semigroups = " << array_of_filesizes[number_of_matching_files]/(ibit/8)/(n-2) << "\n";
          number_of_matching_files++;
        }
        free(all_files[z]);
      }
      free(all_files);
    }
    cout << "----------------------------------------------------------------------" << endl;
    cout << "Total number of bytes = " << total_file_size << ", ";
    cout << "Total number of semigroups = " << total_file_size/(ibit/8)/(n-2) << "\n";

    // Calculating which parts of which files to read for this parallel job
    // and reading the semigroups in to the set input_semigroups

    unsigned long long m;
    m=total_file_size/(ibit/8)/(n-2);
    unsigned long long pStart=ceil(((i-1)*m)/p);
    unsigned long long pEnd=ceil((i*m)/p);
    if (pEnd-pStart!=0)
    {
      unsigned long long bStart;
      unsigned long long bEnd;
      bStart=pStart*(ibit/8)*(n-2);
      bEnd=pEnd*(ibit/8)*(n-2);
      cout << "Reading from semigroup " << pStart << " to " << pEnd << ", ";
      cout << "from byte " << bStart << " to byte " << bEnd << endl;
      cout << "----------------------------------------------------------------------" << endl;
      cout << "Reading in files:" << endl;
      cout << "----------------------------------------------------------------------" << endl;
      int array_index=0;
      unsigned long long running_total=array_of_filesizes[0];
      size_t byteStart;
      size_t byteEnd;
      unsigned long long encoded_word=0;
      while (running_total<=bStart) { array_index++; running_total=running_total+array_of_filesizes[array_index]; }
      if (running_total>=bEnd)
      {
        running_total=running_total-array_of_filesizes[array_index];
        cout << "File " << array_index << ": ";
        cout << "reading " << (bEnd-bStart)/(ibit/8)/(n-2) << " semigroups ";
        cout << "from bytes " << bStart-running_total << " ";
        cout << "to " << bEnd-running_total << ", ";
        cout << "bytes so far = " << bEnd << endl;
        read_file(input_semigroups,vector_of_filenames[array_index],ibit,n-1,k,bStart-running_total,bEnd-running_total);
      }
      else
      {
        cout << "File " << array_index << ": ";
        cout << "reading " << (running_total-bStart)/(ibit/8)/(n-2) << " semigroups ";
        cout << "from bytes " << bStart-running_total+array_of_filesizes[array_index] << " ";
        cout << "to " << array_of_filesizes[array_index] << ", ";
        cout << "bytes so far = " << running_total << endl;
        read_file(input_semigroups,vector_of_filenames[array_index],ibit,n-1,k,bStart-running_total+array_of_filesizes[array_index],array_of_filesizes[array_index]);
        array_index++; running_total=running_total+array_of_filesizes[array_index];
        while (running_total<bEnd)
        {
          cout << "File " << array_index << ": ";
          cout << "reading " << (array_of_filesizes[array_index])/(ibit/8)/(n-2) << " semigroups ";
          cout << "from bytes 0 ";
          cout << "to " << array_of_filesizes[array_index] << ", ";
          cout << "bytes so far = " << running_total << endl;
          read_file(input_semigroups,vector_of_filenames[array_index],ibit,n-1,k,0,array_of_filesizes[array_index]);
          array_index++; running_total=running_total+array_of_filesizes[array_index];
        }
        running_total=running_total-array_of_filesizes[array_index];
        cout << "File " << array_index << ": ";
        cout << "reading " << (bEnd-running_total)/(ibit/8)/(n-2) << " semigroups ";
        cout << "from bytes 0 ";
        cout << "to " << bEnd-running_total << ", ";
        cout << "bytes so far = " << bEnd << endl;
        read_file(input_semigroups,vector_of_filenames[array_index],ibit,n-1,k,0,bEnd-running_total);
      }
      cout << "----------------------------------------------------------------------" << endl;
      cout << "Number of input semigroups = " << input_semigroups.size() << endl;
      cout << "----------------------------------------------------------------------" << endl;
      cout << "Calculating..." << endl;
      cout << "----------------------------------------------------------------------" << endl;
      count=0;
      iter_begin=input_semigroups.begin();
      iter_end=input_semigroups.end();


      count=0;
      iter_begin=input_semigroups.begin();
      iter_end=input_semigroups.end();

      // We now construct all possible forms for new minimal generators...
      // For each semigroup S in input_semgroups
      for (set<semigroup>::iterator S=iter_begin;S!=iter_end;S++)
      {
        new_words.clear();
        // For each word w_i in the set of gaps of S
        for (int wi=0;wi<n-1;wi++)
        {
          // For each g in X_r
          for (int g=0;g<k;g++)
          {
            memcpy(&temp_word.gaps[0],&S->gaps[wi*n],S->gaps[wi*n]+1);
            temp_word.gaps[0]=temp_word.gaps[0]+1;
            temp_word.gaps[temp_word.gaps[0]]=g;
            //  Add [w_i*g] to set of possible new minimal generators
            new_words.insert(temp_word);
          }
        }
        memcpy(&temp_word.gaps[0],&S->gaps[n*(n-2)],S->gaps[n*(n-2)]+1);
        new_words.insert(temp_word);
        iter_frob=new_words.find(temp_word);
        iter_frob++;
        // Checking new word is bigger than Frobenius
        if (iter_frob!=new_words.end())
        {
          gaps.clear();
          for (int l=0;l<n-1;l++)
          {
            memcpy(&temp_semigroup.gaps[(n+1)*l],&S->gaps[n*l],S->gaps[n*l]+1);
            memcpy(&temp_word.gaps[0],&temp_semigroup.gaps[(n+1)*l],temp_semigroup.gaps[(n+1)*l]+1);
            gaps.insert(temp_word);
          }
          for (set<word>::iterator W=iter_frob;W!=new_words.end();W++)
          {
            memcpy(&temp_semigroup.gaps[(n+1)*(n-1)],&W->gaps[0],W->gaps[0]+1);
            // Checking new semigroup is an ideal
            right_word.gaps[0]=temp_semigroup.gaps[(n+1)*(n-1)]-1;
            memcpy(&right_word.gaps[1],&temp_semigroup.gaps[(n+1)*(n-1)+2],temp_semigroup.gaps[(n+1)*(n-1)]-1);
            // If new semigroup is an ideal we now find the orbit of the new semigroup
            if (gaps.find(right_word)!=gaps.end())
            {
              orbit.clear();
              orbit.insert(temp_semigroup);
              // For each permuatation in Sym(k)
              for (unsigned int z=1;z<f;z++)
              {
                for (int l=0;l<n;l++)
                {
                  // Apply permutation to new semigroup
                  temp_perm_semigroup[(n+1)*l]=temp_semigroup.gaps[(n+1)*l];
                  for (int o=1;o<=temp_semigroup.gaps[(n+1)*l];o++) {
                  temp_perm_semigroup[(n+1)*l+o]=A[k*z+temp_semigroup.gaps[(n+1)*l+o]]; }
                }
                // Order the permuted semigroup
                order_semigroup(&ordered_semigroup.gaps[0],&temp_perm_semigroup[0],n);
                // And the ordered semigroup to orbit set
                orbit.insert(ordered_semigroup);
              }
              memcpy(&temp_output_semigroup.gaps[0],&orbit.begin()->gaps[0],(n+1)*n);
              temp_output_semigroup.orbit=orbit.size();
              // Add minimal representative of orbit to set of output_semigroups
              output_semigroups.insert(temp_output_semigroup);
            }
          }
        }
        count++;
      }
      input_semigroups.clear();
      cout << "Number of output semigroups = " << output_semigroups.size() << endl;
      cout << "----------------------------------------------------------------------" << endl;
    }
    delete[] A;
  }

  // Output all the data from output_semigroups in to separate files for each orbit size

  unsigned long long temp_number;

  set<orbit_semigroup>::iterator o_iter_begin=output_semigroups.begin();
  set<orbit_semigroup>::iterator o_iter_end=output_semigroups.end();
  set<orbit_semigroup>::iterator iter1=o_iter_begin; size_t d_1=0;
  set<orbit_semigroup>::iterator iter2=iter1; iter2++;
  size_t size_of_buffer=(obit/8)*n*output_semigroups.size();
  char *buff=new char[size_of_buffer];
  for (int j=1;j<n;j++)
  { temp_number=shortlex_encode(&(iter1->gaps[(n+1)*j]),k);
    memcpy(&buff[(obit/8)*((d_1)*(n-1)+(j-1))],&temp_number,obit/8); }
  unsigned long long d=1;
  while (iter1!=o_iter_end)
  {
    if (iter1->orbit==iter2->orbit) { d++; } else
    {
      stringstream stream_output_file;
      stream_output_file << "ideals/n-" << n << "/n-" << n << "-k-" << k << "-i-" << iter1->orbit << ".bit" << obit;
      if (p!=1)
      {
        stream_output_file << "-";
        for (int j=floor(log10(i));j<floor(log10(p));j++) stream_output_file << "0";
        stream_output_file << i;
      }
      string output_file=stream_output_file.str();
      cout << "File: " << output_file << ", number of semigroups = " << d << "\n";
      ofstream ofile(output_file.c_str(),ios::out|ios::binary);
      ofile.write(&buff[(obit/8)*(((d_1)-d+1)*(n-1))],(obit/8)*(n-1)*d);
      ofile.close();
      d=1;
    }
    if (iter2==o_iter_end) { iter1++; d_1++; } else {
      iter1++; d_1++; iter2++;
      for (int j=1;j<n;j++)
      { temp_number=shortlex_encode(&(iter1->gaps[(n+1)*j]),k); 
        memcpy(&buff[(obit/8)*((d_1)*(n-1)+(j-1))],&temp_number,obit/8);
      }
    }
  }
  output_semigroups.clear();
  delete[] buff;
  cout << "Finished" << endl;
  return(0);
}

