//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  get_poly_ideal.cpp is part of 'FIIOFS'                                  //
//  Copyright (C) 2014 Alex Bailey                                          //
//                                                                          //
//  Licensing information can be found in the README file                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include <sstream>
#include <set>
#include <math.h>
#include <limits>
#include <dirent.h>

using namespace std;

// Function for GCD

long double gcd(long double a, long double b)
{
  if ( a==0 ) return b;
  return gcd(fmod(b,a), a);
}

// Factorial function

unsigned long long factorial(int n) { return (n==1||n==0)?1:factorial(n-1)*n; }

// Function for Stirling numbers of the first kind

long double stirling(int n,int k) {
    if (n==k) return 1;
    if (n==0) return 0;
    return stirling(n-1,k-1)-(n-1)*stirling(n-1,k);
}

// Defining a structure for adding fractions

struct fraction {
  long double num;
  long double den;
  fraction() {}
  fraction(long double num, long double den) : num(num), den(den) {}
  fraction(const fraction &rhs) : num(rhs.num), den(rhs.den) {}
  fraction operator + (const fraction &rhs)  const {
    long double g=gcd(fabs(num*rhs.den+den*rhs.num),den*rhs.den);
    return fraction((num*rhs.den+den*rhs.num)/g, (den*rhs.den)/g);
  }
};

// Defining a structure for orbit data

struct orbit_data {
  int support;
  unsigned long long orbit_size;
  unsigned long long num_of_semigroups;
  bool operator<(const orbit_data &rhs) const {
    if (support<rhs.support) return (1); else {
      if (support>rhs.support) return (0); else {
        if (orbit_size<rhs.orbit_size) return(1); else {
          if (orbit_size>rhs.orbit_size) return (0); else {
            if (num_of_semigroups<rhs.num_of_semigroups) return(1); else { return(0); } } } } }
  }
};

// Main function begins here

int main(int c, char *a[]) {
  
  // Retrieving command line parameters
  // n = index

  static int n=atoi(a[1]);

  // Declaring variables

  set<orbit_data> all_data;
  orbit_data temp_data;
  
  // Finding relevant files in input folder

  stringstream stream_folder_name;
  stream_folder_name << "ideals/n-" << n << "/";
  string folder_name=stream_folder_name.str();

  cout << "Reading from folder: " << folder_name << endl;

  stringstream stream_required_prefix;
  stream_required_prefix << "n-" << n << "-k-";
  string required_prefix=stream_required_prefix.str();

  struct dirent **all_files;
  int x=scandir(folder_name.c_str(),&all_files,0,alphasort);
  if (x<0) perror("scandir");
  else
  {
    for (int z=0;z<x;z++)
    {
      char file_prefix[50]={0};
      memcpy(&file_prefix[0],all_files[z]->d_name,floor(log10(n))+6);
      if (strcmp(required_prefix.c_str(),file_prefix)==0)
      {
        string file_name = folder_name + all_files[z]->d_name;
        ifstream file(file_name.c_str(),ios::in|ios::binary|ios::ate); streampos file_size=file.tellg();
        file.close();
        int pos_k=file_name.find("-k-");
        int pos_i=file_name.find("-i-");
        int pos_bit=file_name.find(".bit");
        // Checking file name is of the correct form
        if (file_name.find("-",pos_bit+1)==string::npos)
        {
          // Retrieving support, orbit size and bitrate from filenames
          temp_data.support=atoi(file_name.substr(pos_k+3,pos_i-pos_k-3).c_str());
          temp_data.orbit_size=atol(file_name.substr(pos_i+3,pos_bit-pos_i-3).c_str());
          temp_data.num_of_semigroups=file_size/((n-1)*(atol(file_name.substr(pos_bit+4,string::npos).c_str())/8));
          cout << "File name = " << all_files[z]->d_name << ", ";
          cout << "file size = " << file_size << ", ";
          cout << "number of semigroups = " << temp_data.num_of_semigroups << endl;
          all_data.insert(temp_data);
        }
      }
      free(all_files[z]);
    }
  }
  free(all_files);

  // Calculating polynomial

  cout << "Calculating polynomial for n=" << n << endl;

  long double c_array[n]; for (int k=0;k<n;k++) c_array[k]=0;
  
  set<orbit_data>::iterator o_iter_begin=all_data.begin();
  set<orbit_data>::iterator o_iter_end=all_data.end();
  set<orbit_data>::iterator iter1=o_iter_begin;
  set<orbit_data>::iterator iter2=iter1; iter2++;
  stringstream output;
  output << "a_" << n << "(FS_r) = (1/" << iter1->support << "!)*( ";
  if (all_data.size()>1)
  {
    while (iter1!=o_iter_end)
    {
      output << iter1->orbit_size << "*" << iter1->num_of_semigroups;
      c_array[iter1->support-1]=c_array[iter1->support-1]+iter1->orbit_size*iter1->num_of_semigroups;
      if (iter1->support==iter2->support and iter2!=o_iter_end) { output << " + "; }
      else
      {
        output << " )*r";
        for (int k=1;k<iter1->support;k++) { output << "*(r-" << k << ")"; }
        if (iter2!=o_iter_end) { output << endl << "+ (1/" << iter2->support << "!)*( "; }
      }
      if (iter2!=o_iter_end) { iter1++; iter2++; } else { iter1++; }
    }
  }
  else
  {
    output << iter1->orbit_size << "*" << iter1->num_of_semigroups << " )*r";
  }

  // Simplifying polynomial

  output.precision(numeric_limits<long double>::digits10);

  output << endl << "Simplifying polynomial" << endl << "a_" << n << "(FS_r)=";

  fraction f[n], temp_fraction;
  for (int j=0;j<n;j++) { f[j].num=0; f[j].den=1; }
  
  if (c_array[n-1]!=0)
  {
    long double g=gcd(c_array[n-1],factorial(n));
    temp_fraction.num=c_array[n-1]/g;
    temp_fraction.den=factorial(n)/g;
    f[n-1]=f[n-1]+temp_fraction;
  }
  
  output << "(" << f[n-1].num;
  if (f[n-1].den!=1) output << "/" << f[n-1].den;
  output << ")*r^" << n;

  for (int j=1;j<n;j++)
  {
    for (int k=n-j;k<=n;k++)
    {
      if (c_array[k-1]!=0)
      {
        long double g=gcd(fabs(c_array[k-1]*stirling(k,n-j)),factorial(k));
        temp_fraction.num=c_array[k-1]*stirling(k,n-j)/g;
        temp_fraction.den=factorial(k)/g;
        f[n-j-1]=f[n-j-1]+temp_fraction;
      }
    }
    if (f[n-j-1].num<0) { output << "-"; } else { output << "+"; }
    output << "(" << fabs(f[n-j-1].num);
    if (f[n-j-1].den!=1) output << "/" << f[n-j-1].den;
    if (n-j==1) output << ")*r"; else output << ")*r^" << n-j;
  }

  cout << output.str() << endl;

return(0);
}
