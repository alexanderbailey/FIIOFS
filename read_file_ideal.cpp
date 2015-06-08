//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  read_file_ideal.cpp is part of 'FIIOFS'                                 //
//  Copyright (C) 2014 Alex Bailey                                          //
//                                                                          //
//  Licensing information can be found in the README file                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <math.h>

using namespace std;

// Function to print the set of gaps of a semigroup

void print_semigroup(const char *V, char n) {
  for (int i=0;i<=n-1;i++)
  {
    cout << "[ ";
    for (int j=1;j<=*(V+i*(n+1));j++) { cout << *(V+i*(n+1)+j)+1 << " "; }
    cout << "]";
  }
  cout << endl;
}

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

// Main function begins here

int main(int c, char *a[]) {

  // Retrieving command line parameters
  // n = index, k = support, i = size of orbit, p = number of chunks to split file, q = which chunk to read

  static int n=atoi(a[1]); static int k=atoi(a[2]); long i=atol(a[3]); long p=atol(a[4]); long q=atol(a[5]);

  stringstream stream_folder_name;
  stream_folder_name << "ideals/n-" << n << "/";
  string folder_name=stream_folder_name.str();

  // Calculating bitrate of file

  int bit;
  char largest_frob[n+1];
  largest_frob[0]=n;
  for (int j=1;j<=n;j++) largest_frob[j]=k-1;
  bit=8*ceil((floor(log(shortlex_encode(&largest_frob[0],k))/log(2))+1)/8);

  // Declaring variables

  char temp_semigroup[(n+1)*n];
  temp_semigroup[0]=1;
  temp_semigroup[1]=0;
  unsigned long long encoded_word;

  // Calcuating which part of the file to read

  size_t file_size;
  stringstream stream_file_name;
  stream_file_name << folder_name << "n-" << n << "-k-" << k << "-i-" << i << ".bit" << bit;
  string file_name=stream_file_name.str();
  ifstream file(file_name.c_str(),ios::in|ios::binary|ios::ate);
  file_size=file.tellg();

  unsigned long long m;
  m=file_size/(bit/8)/(n-1);
  unsigned long long pStart=ceil(((q-1)*m)/p);
  unsigned long long pEnd=ceil((q*m)/p);
  unsigned long long bStart=pStart*(bit/8)*(n-1);
  unsigned long long bEnd=pEnd*(bit/8)*(n-1);

  cout << "File: " << file_name << ", filesize = " << file_size << ", number of semigroups = " << file_size/(bit/8)/(n-1) << endl;
  cout << "Reading chunk " << q << " of " << p << endl;
  cout << "Reading semigroups "; if (pEnd==0) cout << "0"; else cout << pStart; cout << " to " << pEnd << ", ";
  cout << "from byte " << bStart << " to " << bEnd << endl;

  // Reading in q-th chunk of file

  char *data_stream=new char[bEnd-bStart];
  file.seekg(bStart,ios::beg);
  file.read(data_stream,bEnd-bStart);
  file.close();

  for (int j=0;j<pEnd-pStart;j++)
  {
    cout << "0 ";
    for (int l=0;l<n-1;l++)
    { 
      encoded_word=0;
      memcpy(&encoded_word,&data_stream[(bit/8)*((n-1)*j+l)],bit/8);
      cout << encoded_word << " ";
      shortlex_decode(encoded_word,k,&(temp_semigroup[(n+1)*(l+1)]));
    }
    // Printing gaps of semigroup with both shortlex encoding and without
    print_semigroup(&temp_semigroup[0],n);
  }
  delete[] data_stream;
  return(0);
}
