#include <Rcpp.h>
using namespace Rcpp;
using namespace std; 
//Source for C++ perfect sum code problem 


// Prints all subsets of arr[0..n-1] with sum 0. 
// [[Rcpp::export]]
void printAllSubsets(NumericVector numbers,int target,NumericVector partial) 
{ 
  double s = sum(partial); 
  if (s== target) {
    Rcout << partial << "\n"; 
  }
  if(s >=target) {
    return; 
  }
  
  int n = numbers.length(); 
} 


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
arr = c(1,2,3,4,5,10) 
n = length(arr)
sum = 12
printAllSubsets(arr, n, sum) 
*/
