////////////////////////////////////////////////////////////////////////////////
//  isomap algorithm for R
//  
//  created by Kohta Ishikawa (quantumcorgi _at_ gmail.com)
//  date       Jan 3 2012
//  2012 Kohta Ishikawa all rights reserved.
//
//  There might be some bugs in this program.
//  You could tell me if you find it.
////////////////////////////////////////////////////////////////////////////////


#include "risomap.h"
#include "cpplapack.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>

//caluculate a distance matrix
Rcpp::NumericMatrix getDistMat(Rcpp::NumericMatrix data)
{
  int num = data.nrow(); //data number
  int dim = data.ncol(); //data dimension
  Rcpp::NumericMatrix dist(num,num);
  for(int i=0;i<num;i++){
    for(int j=i;j<num;j++){
      if(i==j){
	dist(i,j) = 0;
      }else{
         double d = 0;
         for(int k=0;k<dim;k++){
           d += (data(i,k) - data(j,k))*(data(i,k) - data(j,k));
         }
         dist(i,j) = sqrt(d);
         dist(j,i) = dist(i,j);
      }
    }
  }
  return dist;
}

//function objects for sorting
class PrInc{
public:
  bool operator()(const std::pair<int,double> left, const std::pair<int,double> right){
    return left.second < right.second;
  }
};

class PrDec{
public:
  bool operator()(const std::pair<int,double> left, const std::pair<int,double> right){
    return left.second > right.second;
  }
};

//calculate a k nearest neighbor matrix
CPPL::dsymatrix getKNNMat(Rcpp::NumericMatrix dist, int k)
{
  //knnGraph Matrix
  CPPL::dgematrix knnGraph = CPPL::dgematrix(dist.nrow(),dist.nrow());
  knnGraph.zero();
  //calculate kNN by data points
  int num = dist.nrow();
  for(int i=0;i<num;i++){
    //sort by distance
    std::vector< std::pair<int,double> > knn(num);
    for(int j=0;j<num;j++){
      knn[j] = std::pair<int,double>(j,dist(i,j));
    }
    std::sort(knn.begin(),knn.end(),PrInc());
    //create a matrix
    int cnt = 0;
    for(int j=0;j<num && cnt<k;j++){
      if(i!=knn[j].first){
        knnGraph(i,knn[j].first) = knn[j].second;
        cnt++;
      }
    }
  }

  //symmetrization of matrix
  CPPL::dsymatrix knnGraphSym = CPPL::dsymatrix(dist.nrow());
  knnGraphSym.zero();
  for(int i=0;i<num;i++){
    for(int j=i;j<num;j++){
      if(knnGraph(i,j)==0 && knnGraph(j,i)!=0){
        knnGraphSym(i,j) = knnGraph(j,i);
      }else if(knnGraph(j,i)==0 && knnGraph(i,j)!=0){
        knnGraphSym(j,i) = knnGraph(i,j);
      }else if(i==j){
	knnGraphSym(i,j) = knnGraph(i,j);
      }else if(knnGraph(i,j) == knnGraph(j,i)){
	knnGraphSym(i,j) = knnGraph(i,j);
      }
    }
  }

  return knnGraphSym;
}

//create a geodesic distance matrix from kNN graph matrix
CPPL::dsymatrix getGeodesicDistMat(CPPL::dsymatrix knnGraph)
{
  //calculate minimum distance between each points (Warshall-Floyd algorithm)
  int num = knnGraph.m;
  for(int i=0;i<num;i++){
    for(int j=0;j<num;j++){
      for(int k=0;k<num;k++){
	if(knnGraph(j,i)!=0 && knnGraph(i,k)!=0){
	  if( j!=k && ((knnGraph(j,k) == 0) || (knnGraph(j,i)+knnGraph(i,k) < knnGraph(j,k))) ){
	    knnGraph(j,k) = knnGraph(j,i) + knnGraph(i,k);
	  }
	}
      }
    }
  }
  return knnGraph;
}

//create a kernel matrix by doing centering transformation to geodesic distance matrix
CPPL::dsymatrix getKernelMatrix(CPPL::dsymatrix geodesicDist)
{
  int num = geodesicDist.m;
  //create a centering matrix
  CPPL::dsymatrix H(num);
  for(int i=0;i<num;i++){
    for(int j=i;j<num;j++){
      if(i==j){
        H(i,j) = 1 - 1.0/num;
      }else{
        H(i,j) = -1.0/num;
        H(j,i) = -1.0/num;
      }
    }
  }
  //squaring
  CPPL::dsymatrix sqrGeoDist(num);
  for(int i=0;i<num;i++){
    for(int j=0;j<num;j++){
      sqrGeoDist(i,j) = geodesicDist(i,j) * geodesicDist(i,j);
    }
  }
  //create a kernel matrix
  CPPL::dgematrix K(num,num);
  K = - 0.5 * H * sqrGeoDist * H;
  CPPL::dsymatrix Ksy(num);
  for(int i=0;i<num;i++){
    for(int j=i;j<num;j++){
      Ksy(i,j) = K(i,j);
    }
  }
  return Ksy;
}

//diagonalizing the kernel marix and letting negative eigenvalues zero.
//returning eigenvalues and eigenvectors sorted by eigenvalues in discending order
//(multidimensional scaling method)
std::pair< std::vector<double>, std::vector<CPPL::dcovector> > getCMDS(CPPL::dsymatrix kernelMat)
{
  std::vector<double> wr;
  std::vector<CPPL::dcovector> vrr;
  kernelMat.dsyev(wr,vrr);

  std::vector< std::pair<int,double> > wrVec(wr.size());
  for(int i=0;i<wr.size();i++){
    if(wr[i]<0){
      wr[i] = 0.0;
    }
    wrVec[i] = std::pair<int,double>(i,wr[i]);
  }
  std::sort(wrVec.begin(),wrVec.end(),PrDec());

  std::vector<double> wrSorted(wr.size());
  std::vector<CPPL::dcovector> vrrSorted(wr.size());
  for(int i=0;i<wr.size();i++){
    wrSorted[i] = wrVec[i].second;
    vrrSorted[i] = vrr[wrVec[i].first];
  }
  return std::pair< std::vector<double>, std::vector<CPPL::dcovector> >(wrSorted, vrrSorted);
}

//isomap function for R
RcppExport SEXP risomap(SEXP data, SEXP k)
{
  Rcpp::NumericMatrix rdata(data);
  //create a distance matrix
  Rcpp::NumericMatrix dist = getDistMat(rdata);
  //create a kNN graph matrix
  CPPL::dsymatrix kNNGraph = getKNNMat(dist,Rcpp::as<int>(k));
  //create a geodesic distance matrix
  CPPL::dsymatrix geodesicDist = getGeodesicDistMat(kNNGraph);
  //create a kernel matrix and apply MDS method
  CPPL::dsymatrix kernelMat = getKernelMatrix(geodesicDist);
  std::pair< std::vector<double>, std::vector<CPPL::dcovector> >
    eigenCMDS = getCMDS(kernelMat);

  int num = rdata.nrow();
  //calculate coordinates
  Rcpp::NumericMatrix projection(num,num);
  for(int i=0;i<num;i++){
    for(int j=0;j<num;j++){
      projection(i,j) = sqrt(eigenCMDS.first[j])*eigenCMDS.second[j](i); 
    }
  }
  return Rcpp::wrap(projection);
}


//for test
RcppExport SEXP rcmds(SEXP data)
{
  Rcpp::NumericMatrix dist(data);
  int num = dist.nrow();
  CPPL::dsymatrix distMat(num);
  for(int i=0;i<num;i++){
    for(int j=i;j<num;j++){
      distMat(i,j) = dist(i,j);
    }
  }
  CPPL::dsymatrix kernelMat = getKernelMatrix(distMat);
  std::pair< std::vector<double>, std::vector<CPPL::dcovector> >
    eigenCMDS = getCMDS(kernelMat);
  Rcpp::NumericMatrix result(num,num);
  for(int i=0;i<num;i++){
    for(int j=0;j<num;j++){
      result(i,j) = eigenCMDS.second[j](i);
    }
  }
  return Rcpp::wrap(result);
}
